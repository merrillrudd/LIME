#' Operating model
#'
#' \code{sim_pop} Age-converted-to-length-based operating model specifying true population dynamics
#'
#' @param lh list of life history information, from create_lh_list
#' @param Nyears number of years to simulate
#' @param Fdynamics Specify name of pattern of fishing mortality dynamics, Constant, Endogenous, Ramp, Increasing, or None
#' @param Rdynamics Specify name of pattern of recruitment dynamics, Constant, Pulsed, Pulsed_up, or BH
#' @param Nyears_comp number of years of length composition data
#' @param comp_sample vector with sample sizes of length composition data annually
#' @param init_depl initial depletion; if FALSE, will use F1 from lh list
#' @param nburn number of years of burn-in for operating model
#' @param seed set seed for generating stochastic time series
#' @param modname save model name for true dynamics in named list output

#' @return named list of attributes of true population/data
#' @export
sim_pop <- function(lh, Nyears, Fdynamics, Rdynamics, Nyears_comp, comp_sample, init_depl, nburn, seed, modname){

    ## SB_t = spawning biomass over time
    ## F_t = fishing mortality over time
    ## Cn_at = number of individuals that die from fishing mortality
    ## N_at = abundance by number at age over time

    with(lh, {
    ##########################
    ## Initial calcs
    ##########################

    tyears_only <- nburn+Nyears
    Nyears_real <- Nyears
    nburn <- nburn*nseasons
    Nyears <- Nyears*nseasons
    tyears <- tyears_only*nseasons


    ##########################
    ## Random variables
    ##########################
    set.seed(seed)
    ## recruitment deviations
    RecDev <- rnorm(tyears, mean=-(SigmaR^2)/2, sd=SigmaR)
    ## autocorrelated recruitment deviations
    RecDev_AR <- rep(NA, length(RecDev))
    RecDev_AR[1] <- RecDev[1]
    for(t in 2:length(RecDev)){
        RecDev_AR[t] <- RecDev_AR[t-1]*rho + sqrt(1-rho^2)*RecDev[t]
    }

    ## fishing mortality deviations
    FishDev <- rnorm(tyears, mean=-(SigmaF^2)/2, sd=SigmaF)

    ## abundance index observation error
    IndexDev <- rnorm(tyears, mean=-(SigmaI^2)/2, sd=SigmaI)

    ## catch observation error
    CatchDev <- rnorm(tyears, mean=-(SigmaC^2)/2, sd=SigmaC)

    ##########################
    ## Data objects
    ##########################
    TB_t <- VB_t <- SB_t <- F_t <- R_t <- D_t <- rep(NA, tyears)                               
    Cn_at <- N_at <- N_at0 <- matrix(NA, nrow=length(L_a), ncol=tyears)

    #####################################
    ## Fishing and recruitment dynamics
    #####################################   
    ## reference points
    F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=ages, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, ref=0.4)$root, error=function(e) NA)
    if(init_depl==FALSE) Finit <- F1
    if(init_depl!=FALSE) Finit <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=ages, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, ref=init_depl)$root, error=function(e) NA)
    if(is.na(Finit)) stop("F corresponding to initial depletion does not exist")
    Fmax <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=ages, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, ref=0.2)$root, error=function(e) NA)
    if(is.na(Fmax)) Fmax <- 3

    if(Fdynamics=="Ramp") Framp_t <- c(rep(Finit, nburn), "rampup"=seq(Finit, Fmax, length=floor(Nyears/2)), 
        "peak"=rep(Fmax, floor((Nyears-floor(Nyears/2))/2)), 
        "managed"=rep(Fmax/2, Nyears-floor(Nyears/2)-floor((Nyears-floor(Nyears/2))/2)))
    if(Fdynamics=="Constant") Fconstant_t <- rep(Finit, tyears)
    if(Fdynamics=="Increasing") Finc_t <- c(rep(Finit, nburn), seq(Finit, Fmax, length=Nyears))
    if(Fdynamics=="Decreasing") Fdec_t <- c(rep(Finit, nburn), seq(Finit, 0, length=Nyears))
    if(Fdynamics=="None") F_t <- rep(0, tyears)
    if(Fdynamics=="4010") F_t <- rep(NA, tyears)

    if(Rdynamics=="Pulsed") Rpulse_t <- c(rep(R0, nburn), "initial"=rep(R0, floor(Nyears/3)),
        "pulse_down"=rep(R0/3, floor(Nyears/3)), "pulse_up"=rep(R0, Nyears-(2*floor(Nyears/3))))
    if(Rdynamics=="Pulsed_up") Rpulse_t <- c(rep(R0, nburn), "initial"=rep(R0, floor(Nyears/3)), "pulse_up"=rep(R0*3, floor(Nyears/3)), "pulse_down"=rep(R0, Nyears-(2*floor(Nyears/3))))
    if(Rdynamics=="Constant" | Rdynamics=="AR") Rconstant_t <- rep(R0, tyears)

        if(Fdynamics=="Ramp"){
            F_t <- Framp_t * exp(FishDev)
        }
        if(Fdynamics=="Constant"){
            F_t <- Fconstant_t * exp(FishDev)
        }
        if(Fdynamics=="Increasing"){
            F_t <- Finc_t * exp(FishDev)
        }
        if(Fdynamics=="Decreasing"){
            F_t <- Fdec_t * exp(FishDev)
        }
        if(Fdynamics=="Endogenous"){
            F_t[1] <- Finit
        }
        if(Fdynamics=="4010"){
            F_t[1] <- Finit
        }
        if(Rdynamics=="Constant"){
            R_t <- Rconstant_t/nseasons * exp(RecDev)
        }
        if(Rdynamics=="AR"){
            R_t <- Rconstant_t/nseasons * exp(RecDev_AR)
        }
        if(Rdynamics=="Pulsed"){
            R_t <- Rpulse_t/nseasons * exp(RecDev)
        }
        if(Rdynamics=="Pulsed_up"){
            R_t <- Rpulse_t/nseasons * exp(RecDev)
        }
        if(Rdynamics=="BH"){
            R_t[1] <- R0/nseasons * exp(RecDev[1])
        }

    ## year 1
    for(a in 1:length(L_a)){
        if(a==1){
            N_at[a,1] <- R_t[1]
            N_at0[a,1] <- R_t[1]
        }
        if(a>1 & a<length(L_a)){
            N_at[a,1] <- N_at[a-1,1] * exp(-M - F_t[1] * S_a[a-1])
            N_at0[a,1] <- N_at0[a-1,1] * exp(-M)
        }
        if(a==length(L_a)){
            N_at[a,1] <- (N_at[a-1,1] * exp(-M - F_t[1] * S_a[a])) / (1-exp(-M - F_t[1] * S_a[a]))
            N_at0[a,1] <- (N_at0[a-1,1] * exp(-M))/(1-exp(-M))
        }

    }
    VB_t[1] <- sum(N_at[,1] * W_a * S_a)
    TB_t[1] <- sum(N_at[,1] * W_a)
    SB_t[1] <- sum(N_at[,1] * W_a * Mat_a)
    Cn_at[,1] <- N_at[,1] * (1 - exp(-M - F_t[1] * S_a)) * (F_t[1] *S_a)/(M + F_t[1] * S_a)

    ##########################
    ## Projection
    ##########################
    SB0 <- sum(N_at0[,1]*Mat_a*W_a)
    D_t[1] <- SB_t[1]/SB0

    # D_t <- seq(0,1,length=tyears)
    # F_t <- rep(NA, tyears)
    # F_t[1] <- 0
    # for(t in 2:tyears){
    #             if(D_t[t-1] < 0.10) F_t[t] <- 0
    #             if(D_t[t-1] >= 0.40) F_t[t] <- F40 * exp(FishDev[t])
    #             if(D_t[t-1] >= 0.10 & D_t[t-1] < 0.40) F_t[t] <- ((F40/0.3)*D_t[t-1] - ((0.10*F40)/0.30)) * exp(FishDev[t])
    # }

    for(y in 2:tyears){
        ## fishing effort and recruitment, not dependent on age structure
        if(Fdynamics=="Endogenous"){
            if(y <= nburn) F_t[y] <- Finit
            if(y > nburn) F_t[y] <- F_t[y-1]*(SB_t[y-1]/(Fequil*SB0))^Frate * exp(FishDev[y])
        }
        if(Fdynamics=="4010"){
            if(y <= nburn) F_t[y] <- Finit
            if(y > nburn){
                if(D_t[y-1] < 0.10) F_t[y] <- 0
                # if(D_t[y-1] >= 0.40) F_t[y] <- F40 * exp(FishDev[y])
                # if(D_t[y-1] >= 0.10 & D_t[y-1] < 0.40) F_t[y] <- ((F40/0.3)*D_t[y-1] - ((0.10*F40)/0.30)) * exp(FishDev[y])
                if(D_t[y-1] >= 0.10) F_t[y] <- ((F40/0.3)*D_t[y-1] - ((0.10*F40)/0.30)) * exp(FishDev[y])
            }
        }
        if(Rdynamics=="BH"){
            if(h==1) h_use <- 0.7
            if(h!=1) h_use <- h
            R_t[y] <- (4 * h_use * R0 * SB_t[y-1] / ( SB0*(1-h_use) + SB_t[y-1] * (5*h_use-1)))/nseasons * exp(RecDev[y])
        }

        ## age-structured dynamics
        for(a in 1:length(L_a)){

            if(a==1){
                N_at[a,y] <- R_t[y]
                N_at0[a,y] <- R_t[y]
            }
            if(a>1 & a<length(L_a)){
                N_at[a,y] <- N_at[a-1,y-1] * exp(-M - F_t[y-1] * S_a[a-1])
                N_at0[a,y] <- N_at0[a-1,y-1] * exp(-M)
            }
            if(a==length(L_a)){
                N_at[a,y] <- (N_at[a-1,y-1] * exp(-M - F_t[y-1] * S_a[a-1])) + (N_at[a,y-1] * exp(-M - F_t[y-1] * S_a[a]))
                N_at0[a,y] <- (N_at0[a-1,y-1] * exp(-M)) + (N_at0[a,y-1] * exp(-M))
            }

            ## spawning biomass
            SB_t[y] <- sum((N_at[,y] * W_a * Mat_a))
            VB_t[y] <- sum(N_at[,y] * W_a * S_a)
            TB_t[y] <- sum(N_at[,y] * W_a)

            ## catch
            Cn_at[,y] <- N_at[,y] * (1 - exp(-M - F_t[y] * S_a)) * (F_t[y] * S_a)/ (M + F_t[y] * S_a)
            D_t <- SB_t/SB0


        }
    }
    Cn_t <- colSums(Cn_at)
    Cw_t <- colSums(Cn_at * W_a)
    N_t <- colSums(N_at[-1,])

    I_t <- qcoef * TB_t #* exp(IndexDev - (SigmaI^2)/2)
    C_t <- Cn_t #* exp(CatchDev - (SigmaC^2)/2)

    ## age to length comp
    obs_per_year <- rep(comp_sample, tyears)
    LFinfo <- AgeToLengthComp(lh=lh, tyears=tyears, N_at=N_at, comp_sample=obs_per_year)
    LF0info <- AgeToLengthComp(lh=lh, tyears=tyears, N_at=N_at0, comp_sample=obs_per_year)

    plba <- LFinfo$plba
    plb <- LFinfo$plb
    page <- LFinfo$page
    LF <- LFinfo$LF
    LF0 <- LF0info$LF

    ########################################################
    ## Expected mean length in catch 
    ########################################################
    ML_t <- vector(length=tyears)
    for(y in 1:tyears){
        vul_pop <- sum(N_at[,y]*S_a)
        vul_lengths <- sum(vul_pop*plb[y,]*mids)
        ML_t[y] <- vul_lengths/vul_pop
    }

    ########################################################
    ## cut out burn-in
    ########################################################

    I_tout <- I_t[-c(1:nburn)]
    C_tout <- C_t[-c(1:nburn)]
    Cw_tout <- Cw_t[-c(1:nburn)]
            names(C_tout) <- names(Cw_tout) <- names(I_tout) <- 1:Nyears

    LFout <- LF[-c(1:nburn),]
        rownames(LFout) <- 1:Nyears
    LF0out <- LF0[-c(1:nburn),]
        rownames(LF0out) <- 1:Nyears

    R_tout <- R_t[-c(1:nburn)]
    N_tout <- N_t[-c(1:nburn)]
    SB_tout <- SB_t[-c(1:nburn)]
    TB_tout <- TB_t[-c(1:nburn)]
    VB_tout <- VB_t[-c(1:nburn)]
    D_tout <- D_t[-c(1:nburn)]
    F_tout <- F_t[-c(1:nburn)]
    ML_tout <- ML_t[-c(1:nburn)]
    N_atout <- N_at[,-c(1:nburn)]
    N_at0out <- N_at0[,-c(1:nburn)]

        # LFindex <- (Nyears-Nyears_comp+1):Nyears
        LFindex <- rev(seq(Nyears, by=-nseasons, length.out=Nyears_comp))
        # LFindex <- rev(seq(Nyears-nseasons+1, by=-nseasons, length.out=Nyears_comp))
        LFout <- LFout[LFindex,]
        LF0out <- LF0out[LFindex,]
        if(is.vector(LFout)==FALSE) colnames(LFout) <- highs
        if(is.vector(LF0out)==FALSE) colnames(LF0out) <- highs
        if(is.vector(LFout)){
            LFout <- t(as.matrix(LFout))
            rownames(LFout) <- LFindex
        }
        if(is.vector(LF0out)){
            LF0out <- t(as.matrix(LF0out))
            rownames(LF0out) <- LFindex
        }

    ## static SPR
    SPR_t <- sapply(1:length(F_tout), function(x) calc_ref(ages=ages, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, F=F_tout[x]))
    SPR <- SPR_t[length(SPR_t)]

    ## outputs
    lh$I_t <- I_tout
    lh$C_t <- C_tout
    lh$Cw_t <- Cw_tout
    lh$DataScenario <- modname
    lh$LF <- LFout
    lh$LF0 <- LF0out
    lh$R_t <- R_tout
    lh$N_t <- N_tout
    lh$SB_t <- SB_tout
    lh$D_t <- D_tout
    lh$F_t <- F_tout
    lh$ML_t <- ML_tout
    lh$N_at <- N_atout
    lh$N_at0 <- N_at0out
    lh$plb <- plb
    lh$plba <- plba
    lh$page <- page
    lh$SPR <- SPR
    lh$SPR_t <- SPR_t
    lh$VB_t <- VB_tout
    lh$TB_t <- TB_tout
    lh$nlbins <- length(mids)
    lh$Nyears <- Nyears
    lh$years <- 1:Nyears
    lh$obs_per_year <- obs_per_year
    if(Rdynamics!="AR") lh$RecDev <- RecDev
    if(Rdynamics=="AR") lh$RecDev <- RecDev_AR
    lh$FishDev <- FishDev

    return(lh)

}) ## end with function

}
