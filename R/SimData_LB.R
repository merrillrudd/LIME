#' Operating model
#'
#' \code{SimData_LB} Age-converted-to-length-based operating model specifying true population dynamics
#'
#' @param Nyears number of years to simulate
#' @param AgeMax maximum age in population
#' @param SigmaR Specified level of recruitment variation (as a standard deviation)
#' @param M natural mortality rate
#' @param F1 Fishing mortality level in first year
#' @param S_a vector of selectivity at age
#' @param h steepness
#' @param qcoef catchability coefficient for abundance index
#' @param Frate parameter in Endogenous fishing mortality, see Thorson et al. 2014 CJFAS effort dynamics
#' @param Fequil parameter in Endogenous fishing mortality, see Thorson et al. 2014 CJFAS effort dynamics
#' @param SigmaF specified level of variation in fishing mortality annually (as a standard deviation)
#' @param Fdynamics Specify name of pattern of fishing mortality dynamics, "Constant", "Endogenous", "Ramped", "Increasing", or "None"
#' @param Rdynamics Specify name of pattern of recruitment dynamics, "Constant", "Pulsed", "Pulsed_up", or "BH"
#' @param R0 equilibrium, initial recruitment
#' @param Fmax maximum possible fishing mortality
#' @param CVlen variation in growth at age as a coefficient of variation
#' @param mids midpoint of length bins
#' @param highs high end of length bins
#' @param lows low end of length bins
#' @param W_a weight at age
#' @param L_a length-at-age
#' @param linf asymptotic length
#' @param vbk Brody growth coefficient
#' @param Mat_a maturity at age
#' @param M50 Age at 50 percent maturity
#' @param comp_sample number of indiviuals sampled for lenght composition
#' @param Nyears_comp number of years of length composition data
#' @param alt_yrs only some years sampled for data? TRUE or FALSE (default FALSE)
#' @param sample only a sample of catch is reported? default FALSE
#' @param nburn number of years of burn-in for operating model
#' @param seed set seed for generating stochastic time series
#' @param modname save model name for true dynamics in named list output


#' @return named list of attributes of true population/data
#' @export
SimData_LB <- function(Nyears, AgeMax, SigmaR, M, F1, S_a, h, qcoef,
    Frate, Fequil, SigmaF, Fdynamics, Rdynamics, R0, Fmax, CVlen, mids, highs, lows,
    W_a, L_a, linf, vbk, Mat_a, M50, comp_sample, Nyears_comp, alt_yrs=FALSE, sample=FALSE,
    nburn, seed, modname){

    ## SB_t = spawning biomass over time
    ## F_t = fishing mortality over time
    ## Cn_at = number of individuals that die from fishing mortality
    ## N_at = abundance by number at age over time

    ##########################
    ## Initial calcs
    ##########################

    tyears <- nburn+Nyears


    ##########################
    ## Random variables
    ##########################
    set.seed(seed)
    RecDev <- rnorm(tyears, mean=-SigmaR^2/2, sd=SigmaR)
    FishDev <- rnorm(tyears, mean=-SigmaF^2/2, sd=SigmaF)
    EffDev <- rnorm(tyears, mean=-SigmaF^2/2, sd=SigmaF)

    ##########################
    ## Data objects
    ##########################
    TB_t <- VB_t <- SB_t <- F_t <- R_t <- rep(NA, tyears)                               
    Cn_at <- N_at <- matrix(NA, nrow=length(L_a), ncol=tyears)

    #####################################
    ## Fishing and recruitment dynamics
    #####################################   

    if(Fdynamics=="Ramp") Framp_t <- c(rep(F1, nburn), "rampup"=seq(F1, Fmax, length=floor(Nyears/2)), 
        "peak"=rep(Fmax, floor((Nyears-floor(Nyears/2))/2)), 
        "managed"=rep(Fmax/3, Nyears-floor(Nyears/2)-floor((Nyears-floor(Nyears/2))/2)))
    if(Fdynamics=="Constant") Fconstant_t <- c(rep(F1, nburn), rep(Fequil, Nyears))
    if(Fdynamics=="Increasing") Finc_t <- c(rep(F1, nburn), seq(F1, Fmax, length=Nyears))
    if(Fdynamics=="None") F_t <- rep(0, tyears)

    if(Rdynamics=="Pulsed") Rpulse_t <- c(rep(R0, nburn), "initial"=rep(R0, floor(Nyears/3)),
        "pulse_down"=rep(R0/3, floor(Nyears/3)), "pulse_up"=rep(R0, Nyears-floor(Nyears/3)))
    if(Rdynamics=="Pulsed_up") Rpulse_t <- c(rep(R0, nburn), "initial"=rep(R0, floor(Nyears/3)), "pulse_up"=rep(R0*3, floor(Nyears/3)), "pulse_down"=rep(R0, Nyears-floor(Nyears/3)))
    if(Rdynamics=="Constant") Rconstant_t <- rep(R0, tyears)

        if(Fdynamics=="Ramp"){
            F_t <- Framp_t * exp(FishDev)
        }
        if(Fdynamics=="Constant"){
            F_t <- Fconstant_t * exp(FishDev)
        }
        if(Fdynamics=="Increasing"){
            F_t <- Finc_t * exp(FishDev)
        }
        if(Fdynamics=="Endogenous"){
            F_t[1] <- F1
        }
        if(Rdynamics=="Constant"){
            R_t <- Rconstant_t * exp(RecDev)
        }
        if(Rdynamics=="Pulsed"){
            R_t <- Rpulse_t * exp(RecDev)
        }
        if(Rdynamics=="Pulsed_up"){
            R_t <- Rpulse_t * exp(RecDev)
        }
        if(Rdynamics=="BH"){
            R_t[1] <- R0
        }

    ## year 1
    for(a in 1:length(L_a)){
        if(a==1){
            N_at[a,1] <- R_t[1]
        }
        if(a>1 & a<length(L_a)){
            N_at[a,1] <- N_at[a-1,1]*exp(-M-F_t[1]*S_a[a-1])
        }
        if(a==length(L_a)){
            N_at[a,1] <- (N_at[a-1,1]*exp(-M-F_t[1]*S_a[a-1]))/(1-exp(-M-F_t[1]*S_a[a-1]))
        }

    }
    VB_t[1] <- sum(N_at[,1] * W_a * S_a)
    TB_t[1] <- sum(N_at[,1] * W_a)
    SB_t[1] <- sum(N_at[,1] * W_a * Mat_a)
    Cn_at[,1] <- N_at[,1] * (1-exp(-M - F_t[1]*S_a)) * (F_t[1]*S_a)/(M+F_t[1]*S_a)

    ##########################
    ## Projection
    ##########################
    Na0 <- rep(NA, length(W_a))
        if(Rdynamics=="Pulsed"){
            R0 <- median(Rpulse_t[-c(1:nburn)])
        }
    Na0[1] <- R0
    for(a in 2:length(W_a)){
        Na0[a] <- R0 * exp(-M*(a-1))
    }
    SB0 <- sum(Na0*Mat_a*W_a)

    for(y in 2:tyears){
        ## fishing effort and recruitment, not dependent on age structure
        if(Fdynamics=="Endogenous"){
            if(y <= nburn) F_t[y] <- F1
            if(y > nburn) F_t[y] <- F_t[y-1]*(SB_t[y-1]/(Fequil*SB0))^Frate * exp(FishDev[y])
        }

        ## age-structured dynamics
        for(a in 1:length(L_a)){
            if(a>1 & a<length(L_a)){
                N_at[a,y] <- N_at[a-1,y-1]*exp(-M-F_t[y-1]*S_a[a-1])
            }
            if(a==length(L_a)){
                N_at[a,y] <- (N_at[a-1,y-1] + N_at[a,y-1])*exp(-M-F_t[y-1]*S_a[a-1])
            }

            ## spawning biomass
            SB_t[y] <- sum((N_at[,y] * W_a * Mat_a))
            VB_t[y] <- sum(N_at[,y] * W_a * S_a)
            TB_t[y] <- sum(N_at[,y] * W_a)

            ## catch
            Cn_at[,y] <- N_at[,y] * (1-exp(-M-F_t[y]*S_a)) * (F_t[y]*S_a)/(M+F_t[y]*S_a)

            if(Rdynamics=="BH"){
                if(h==1) h_use <- 0.7
                if(h!=1) h_use <- h
                R_t[y] <- (4 * h_use * R0 * SB_t[y] / ( SB0*(1-h_use) + SB_t[y]*(5*h_use-1))) * exp(RecDev[y])
            }

            if(a==1){
                N_at[a,y] <- R_t[y]
            }
        }
    }
    Cn_t <- colSums(Cn_at)
    Cw_t <- colSums(Cn_at * W_a)
    N_t <- colSums(N_at[-1,])
    D_t <- SB_t/SB0

    if(sample!=FALSE){
        C_t <- Cn_t*sample
        Cw_t <- Cw_t*sample
    }
    if(sample==FALSE){
        C_t <- Cn_t
        Cw_t <- Cw_t
    }
    CPUE_t <- qcoef * TB_t * exp(EffDev)

    ## age to length comp
    LFinfo <- AgeToLengthComp(L_a=L_a, CVlen=CVlen, highs=highs, lows=lows, tyears=tyears, N_at=N_at, S_a=S_a, comp_sample=rep(comp_sample, tyears))

    plba <- LFinfo$plba
    plb <- LFinfo$plb
    page <- LFinfo$page
    LF <- LFinfo$LF

    ########################################################
    ## True mean length in vulnerable population each year 
    ########################################################
    L_t <- vector(length=tyears)
    for(y in 1:tyears){
        vul_pop <- sum(N_at[,y]*S_a)
        vul_lengths <- sum(vul_pop*plb[y,]*mids)
        L_t[y] <- vul_lengths/vul_pop
    }

    ########################################################
    ## cut out burn-in
    ########################################################

    CPUE_tout <- CPUE_t[-c(1:nburn)]
    C_tout <- C_t[-c(1:nburn)]
    Cw_tout <- Cw_t[-c(1:nburn)]
            names(C_tout) <- names(Cw_tout) <- names(CPUE_tout) <- 1:Nyears

    LFout <- LF[-c(1:nburn),]
        rownames(LFout) <- 1:Nyears

    R_tout <- R_t[-c(1:nburn)]
    N_tout <- N_t[-c(1:nburn)]
    SB_tout <- SB_t[-c(1:nburn)]
    TB_tout <- TB_t[-c(1:nburn)]
    VB_tout <- VB_t[-c(1:nburn)]
    D_tout <- D_t[-c(1:nburn)]
    F_tout <- F_t[-c(1:nburn)]
    L_tout <- L_t[-c(1:nburn)]
    N_atout <- N_at[,-c(1:nburn)]

    if(alt_yrs==TRUE){
        yrs <- Nyears:1
        index <- rep(c(1,1,0,0), Nyears/4)
        yr_vec <- rev(yrs[which(index==1)])

        C_tout <- C_tout[which(names(C_tout) %in% yr_vec)]
        Cw_tout <- Cw_tout[which(names(Cw_tout) %in% yr_vec)]
        CPUE_tout <- CPUE_tout[which(names(CPUE_tout) %in% yr_vec)]
        LFout <- LFout[which(rownames(LFout) %in% yr_vec),]
    }
        LFindex <- (Nyears-Nyears_comp+1):Nyears
        LFout <- LFout[LFindex,]
        if(is.vector(LFout)){
            LFout <- t(as.matrix(LFout))
            rownames(LFout) <- (Nyears-Nyears_comp+1):Nyears
        }

    ## static SPR
    SPR_t <- sapply(1:length(F_tout), function(x) calc_ref(Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, F=F_tout[x]))
    SPR <- SPR_t[length(SPR_t)]
    
    lbins <- lows
    if(Nyears_comp>1) ML_t <- sapply(1:nrow(LFout), function(x) sum(LFout[x,]*lbins)/sum(LFout[x,]))
    if(Nyears_comp==1) ML_t <- sum(LFout*lbins)/sum(LFout)

    DataList <- list("I_t"=CPUE_tout, "C_t"=C_tout, "Cw_t"=Cw_tout, "DataScenario"=modname,
        "LF"=LFout, "SigmaR"=SigmaR, "R_t"=R_tout, "N_t"=N_tout, "SB_t"=SB_tout, "D_t"=D_tout, "F_t"=F_tout, 
        "L_t"=L_tout, "N_at"=N_atout, "M50"=M50, "Mat_a"=Mat_a, "SB0"=SB0, "Nyears"=Nyears, "L_a"=L_a,
        "W_a"=W_a, "AgeMax"=AgeMax, "M"=M, "S_a"=S_a, "plb"=plb, "plba"=plba, "page"=page, "R0"=R0, 
        "SPR"=SPR, "SPR_t"=SPR_t, "VB_t"=VB_tout, "TB_t"=TB_tout,
        "ML_t"=ML_t, "nlbins"=length(mids))

    return(DataList)
}
