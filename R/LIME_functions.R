#' Age-converted-to-length structure
#'
#' \code{AgeToLengthComp} Converts vulnerable numbers at age to length composition of catch
#'
#' @param L_a vector, growth curve: predicted length at each age
#' @param CVlen coefficient of variation for growth curve
#' @param highs vector of upper end of length bins
#' @param lows vector of lower end of length bins
#' @param tyears number of years of data
#' @param N_at matrix of numbers in the population at each age over time
#' @param S_a vector of selectivity at age
#' @param comp_sample vector of number of individuals sampled each year (set as 1 for proportions)
#'
#' @return data frame - number of individuals in each length bin in each year
#' @export
AgeToLengthComp <- function(L_a, CVlen, highs, lows, tyears, N_at, S_a, comp_sample){

	################################################
	## Probability being in a length bin given age
	################################################
    lbprobs <- function(mnl,sdl) return(pnorm(highs,mnl,sdl)-pnorm(lows,mnl,sdl))
    vlprobs <- Vectorize(lbprobs,vectorize.args=c("mnl","sdl"))
    plba <- t(vlprobs(L_a,L_a*CVlen))
    plba <- plba/rowSums(plba)

    ################################################
	## Probability being in harvested at an age
	################################################
    page <- matrix(ncol=dim(plba)[1], nrow=tyears)
    for (y in 1:tyears) page[y,] <- N_at[,y] * S_a
    page <- page/rowSums(page)    

    ################################################
	## Probability of sampling a given length bin
	################################################
    plb <- matrix(ncol=length(highs), nrow=tyears)
    for (y in 1:tyears) plb[y,] <- page[y,] %*% plba
    plb <- plb/rowSums(plb)    

    #######################
	## Length frequencies 
	#######################
    LF <- array(0, dim=dim(plb))
      rownames(LF) <- 1:tyears
    for(y in 1:tyears){
    	LF[y,] <- rmultinom(n=1, size=comp_sample[y], prob=plb[y,])
    }

    Outs <- NULL
    Outs$plba <- plba
    Outs$plb <- plb
    Outs$page <- page
    Outs$LF <- LF

    return(Outs)
}


#' Calculate biological reference points
#'
#' \code{Calc_derived_quants} calculates derived quanties for status or productivity
#'
#' @param Obj, the fitted TMB object

#' @return List, a tagged list of potentially useful benchmarks
#' @export
Calc_derived_quants = function( Obj ){
  # Extract elements
  Data = Obj$env$data
  ParHat = Obj$env$parList()
  Report = Obj$report()

  SPR <- with(Report, calc_ref(Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, F=F_t[length(F_t)], ref=FALSE))
  F30 <- tryCatch(with(Report, uniroot(calc_ref, lower=0, upper=50, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, ref=0.3)$root), error=function(e) NA)
  F40 <- tryCatch(with(Report, uniroot(calc_ref, lower=0, upper=50, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, ref=0.4)$root), error=function(e) NA)
  FF30 <- FF40 <- NULL
  if(is.na(F30)==FALSE) FF30 <- Report$F_t[length(Report$F_t)]/F30
  if(is.na(F40)==FALSE) FF40 <- Report$F_t[length(Report$F_t)]/F40

  # Total biomass
  TB_t = as.vector( Report$W_a %*% t(Report$N_ta) )
  Cw_t <- as.vector(Report$W_a %*% t(Report$Cn_ta));

  # MSY calculations
  Yield_Fn = function( Fmean, Return_type="Yield" ){
    # Modify data
    Data_new = Data
    Data_new[["n_t"]] = 1000
    Data_new[["n_c"]] = 1000
    Data_new[["n_i"]] = 1000
    Data_new[["n_lc"]] = 1000
    Data_new[["T_yrs"]] = 1:1000
    Data_new[["C_yrs"]] = 1:1000
    Data_new[["I_yrs"]] = 1:1000
    Data_new[["LC_yrs"]] = 1:1000
    Data_new[["obs_per_yr"]] = rep(Data[["obs_per_yr"]][1], 1000)
    Data_new[["RecDev_biasadj"]] = rep(0, Data_new[["n_t"]])
    Data_new[["C_t"]] = rep(1, Data_new[["n_t"]])
      names(Data_new[["C_t"]]) <- Data_new[["C_yrs"]]
    Data_new[["LF"]] = matrix(0, nrow=Data_new[["n_t"]], ncol=Data_new[["n_lb"]])
      rownames(Data_new[["LF"]]) <- Data_new[["LC_yrs"]]
    Data_new[["I_t"]] = cbind(1, rep(1,Data_new[["n_t"]]))
      names(Data_new[["I_t"]]) <- Data_new[["I_yrs"]]
    # Modify parameters
    ParHat_new = ParHat
    ParHat_new[["log_F_t_input"]] = rep( log(Fmean+1e-10), Data_new[["n_t"]])
    ParHat_new[["Nu_input"]] = rep(0, Data_new[["n_t"]])
    Obj_new = MakeADFun(data=Data_new[1:length(Data_new)], parameters=ParHat_new[1:length(ParHat_new)], inner.control=list(maxit=1e3), DLL="LIME" )
    # Extract and return stuff
    Report_new = Obj_new$report()
    if(Return_type=="Yield") Return = rev(Report_new$Cw_t_hat)[1]
    if(Return_type=="Report") Return = Report_new
    return(Return)
  }

  # Calculate Fmsy
  Fmsy = optimize( f=Yield_Fn, interval=c(0,3), maximum=TRUE)$maximum
  Report_msy = Yield_Fn( Fmean=Fmsy, Return_type="Report" )
  Report_0 = Yield_Fn( Fmean=0, Return_type="Report" )
  TBmsy = rev(as.vector(Report$W_a %*% t(Report_msy$N_ta)))[1]
  SBmsy = rev(Report_msy$SB_t)[1]
  MSY = rev(Report_msy$Cw_t_hat)[1]
  TB0 = rev(as.vector(Report$W_a %*% t(Report_0$N_ta)))[1]
  SB0 = rev(Report_0$SB_t)[1]

  # Return
  Return <- list("SPR"=SPR, "F30"=F30, "F40"=F40, "FF30"=FF30, "FF40"=FF40, "Fmsy"=Fmsy, "SB0"=SB0, "TB0"=TB0, "TB_t"=TB_t, "SB_t"=Report$SB_t, "MSY"=MSY, "TBmsy"=TBmsy, "SBmsy"=SBmsy)
  return( Return )
}


#' Per recruit calculation and spawning potential ratio
#'
#' \code{calc_ref} Calculates derived spawning potential ratio: lifetime total egg production in fished:unfished states
#'
#' @param Mat_a from report file / true file for simulation
#' @param W_a from report file / true file for simulation
#' @param M from report file / true file for simulation
#' @param S_a from report file / true file for simulation
#' @param F typically terminal estimated/true F, can be any year
#' @param ref FALSE outputs SPR, ref= a value between 0 and 1 can be used with uniroot to find the F at which SPR=ref

#' @return List, a tagged list of potentially useful benchmarks
#' @details Use this function with uniroot to find the value of F that results in SPR matching the specified reference value (e.g. 0.30 to find F30)
#' @export
calc_ref <- function(Mat_a, W_a, M, S_a, F, ref=FALSE){

    ## calculate spawning biomass per recruit in fished and unfished condition
    ## a function of specified level of fishing mortality and ability to estimate selectivity parameters
        Na0 <- Naf <- rep(NA, length(W_a))
        Na0[1] <- Naf[1] <- 1
        for(a in 2:length(W_a)){
            if(a<length(W_a)){
                Na0[a] <- Na0[a-1]*exp(-M)
                Naf[a] <- Naf[a-1]*exp(-M-S_a[a-1]*F)
            }
            if(a==length(W_a)){
                Na0[a] <- (Na0[a-1]*exp(-M))/(1-exp(-M))
                Naf[a] <- (Naf[a-1]*exp(-M-F*S_a[a-1]))/(1-exp(-M-F*S_a[a-1]))
            }
        }

        ## ignore recruits
        SB0 <- sum(Na0[-1]*Mat_a[-1]*W_a[-1])
        SBf <- sum(Naf[-1]*Mat_a[-1]*W_a[-1])

        ## automatically returns SPR
        SPR <- SBf/SB0
        if(ref==FALSE) return(SPR)
            
        ## can use uniroot function on call to calc_ref to calculate the fishing mortality rate that results in a specified SPR, then compare current fishing mortality with this reference point
        if(ref!=FALSE){
            diff <- ref - SPR
            return(diff)
        }
}


#' Calculate relative error
#'
#' \code{calcRE} Calculates the relative error between true and estimated values of SPR and F/Fref
#'
#' @param modpath_vec vector of directories to search for saved true and estimated population parameters
#' @param itervec vector of iterations to check for results
#' @param value options: "SPR" or "FFref" (FFref is defined as F/F40)
#' @param yr if timeseries==FALSE, specify the year for which to take the relative error
#' @param timeseries default: timeseries=FALSE, if timeseries=TRUE, takes the relative error for all years up to the specified year using the 'yr' argument

#' @return list with Relative Error, Squared Error, Estimation Error, and vector flagging convergence issues
#' @details Only set up to run with a simulation study where there are multiple iterations of the model runs. 
#' @export
    calcRE <- function(modpath_vec, itervec, value, yr, timeseries=FALSE){
        if(timeseries==FALSE){
            RE <- SQerr <- EE <- matrix(NA, nrow=length(itervec), ncol=length(modpath_vec))
            converge <- matrix(0, nrow=length(itervec), ncol=length(modpath_vec))
        }
        if(timeseries==TRUE){
            RE <- SQerr <- EE <- array(NA, dim=c(yr, length(modpath_vec), length(itervec)))
            converge <- matrix(0, nrow=length(itervec), ncol=length(modpath_vec))
        }
        for(m in 1:length(modpath_vec)){
            for(iter in itervec){
                ## report file
                if(grepl("LBSPR", modpath_vec[m])==FALSE){
                    if(file.exists(file.path(modpath_vec[m], iter, "Report.rds"))) Rep <- readRDS(file.path(modpath_vec[m], iter, "Report.rds"))
                    if(file.exists(file.path(modpath_vec[m], iter, "NAs_final_gradient.txt")) | file.exists(file.path(modpath_vec[m], iter, "high_final_gradient.txt"))){
                        converge[iter,m] <- 1
                        next
                    } 
                    if(file.exists(file.path(modpath_vec[m], iter, "Report.rds"))==FALSE) next
                }
                if(grepl("LBSPR", modpath_vec[m])){
                    if(length(which(grepl("LBSPR", list.files(file.path(modpath_vec[m], iter)))))==1) Rep <- readRDS(file.path(modpath_vec[m], iter, "LBSPR_results.rds"))
                    if(file.exists(file.path(modpath_vec[m], iter, "non_convergence.txt"))){
                        converge[iter,m] <- 1
                        next
                    }
                }   

                ## estimates of SPR
                if(value=="SPR"){   

                    ## estimated values
                    if(grepl("LBSPR", modpath_vec[m])==FALSE){
                        if(timeseries==FALSE) Est <- with(Rep, calc_ref(Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, F=F_t[yr], ref=FALSE))
                        if(timeseries==TRUE) Est <- with(Rep, sapply(1:yr, function(x) calc_ref(Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, F=F_t[x], ref=FALSE)))
                    }
                    if(grepl("LBSPR", modpath_vec[m])){
                        Est <- Rep$SPR
                    }   

                
                    ## true values
                    if(file.exists(file.path(modpath_vec[m], iter, "SPR_site.rds"))){
                        True_file <- readRDS(file.path(modpath_vec[m], iter, "SPR_site.rds"))
                        if(timeseries==FALSE) True <- mean(True_file[yr,])
                        if(timeseries==TRUE) True <- rowMeans(True_file)
                    }
                    if(file.exists(file.path(modpath_vec[m], iter, "SPR_site.rds"))==FALSE){
                        True_file <-  readRDS(file.path(modpath_vec[m], iter, "True.rds"))
                        if(timeseries==FALSE) True <- True_file$SPR
                        if(timeseries==TRUE) True <- True_file$SPR_t
                    }
                }   

                ## estimates of F/Fref
                if(value=="FFref"){ 

                    if(grepl("LBSPR", modpath_vec[m])==FALSE){
                        Est <- with(Rep, F_t[yr]/(uniroot(calc_ref, lower=0, upper=50, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, ref=0.4)$root))
                    }
                    if(grepl("LBSPR", modpath_vec[m])){
                        stop("Cannot calculate F-based reference point for LBSPR")
                    }
        
                    
                    ## true
                    True_file <- readRDS(file.path(modpath_vec[m], iter, "True.rds"))
                    True <- with(True_file, F_t[yr]/(uniroot(calc_ref, lower=0, upper=50, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, ref=0.4)$root))
                }   
    

                if(timeseries==FALSE){
                    RE[iter,m] <- (Est[length(Est)] - True[length(True)])/True[length(True)]
                    SQerr[iter,m] <- (Est[length(Est)] - True[length(True)])^2
                    EE[iter,m] <- log(Est[length(Est)]) - log(True[length(True)])
                }
                if(timeseries==TRUE){
                    RE[iter,m,] <- (Est - True)/True
                    SQerr[iter,m,] <- (Est - True)^2
                    EE[iter,m,] <- log(Est) - log(True)
                } 

                rm(Est)
                rm(True)
                rm(Rep)               
            }
            
        }

        Outs <- NULL
        Outs$RelErr <- RE
        Outs$SqErr <- SQerr
        Outs$EstErr <- EE
        Outs$nonconvergence <- converge

        return(Outs)
    }

#' Check for identifiability of fixed effects -- frm TMBhelpers (slightly adjusted by MBR)
#'
#' \code{Check_Identifiable} calculates the matrix of second-derivatives of the marginal likelihood w.r.t. fixed effects, to see if any linear combinations are unidentifiable
#'
#' @param obj, The compiled TMB object
#'
#' @return A tagged list of the hessian and the message
#' @details Slight adjustment made in the RowMax calculation - when the eigenvector was a vector, the apply function was not working 

#' @export
Check_Identifiable2 = function( obj ){
  # Finite-different hessian
  ParHat = TMBhelper:::extract_fixed( obj )
  List = NULL
  List[["Hess"]] = optimHess( par=ParHat, fn=obj$fn, gr=obj$gr )

  # Check eigendecomposition
  List[["Eigen"]] = eigen( List[["Hess"]] )
  List[["WhichBad"]] = which( List[["Eigen"]]$values < sqrt(.Machine$double.eps) )

  # Check for parameters
  RowMax = tryCatch(apply( List[["Eigen"]]$vectors[,List[["WhichBad"]]], MARGIN=1, FUN=function(vec){max(abs(vec))} ), error=function(e) NA)
  if(all(is.na(RowMax))) RowMax <- sapply(1:length(List[["Eigen"]]$vectors[,List[["WhichBad"]]]), function(x) max(abs(vec)))
  List[["BadParams"]] = data.frame("Param"=names(obj$par), "MLE"=ParHat, ifelse(RowMax>0.1, "Bad","OK"))

  # Message
  if( length(List[["WhichBad"]])==0 ){
    message( "All parameters are identifiable" )
  }else{
    print( List[["BadParams"]] )
  }

  # Return
  return( invisible(List) )
}

#' Saved life history values
#'
#' \code{choose_lh_list} Fixed lists of life history/starting value information for pre-loaded species
#'
#' @param species species code name
#' @param selex "asymptotic"= assume asymptotic selectivity, "dome"= assume dome-shaped selectivity
#' 
#' @return List, a tagged list of life history/starting value information
#' @export
choose_lh_list <- function(species, selex, param_adjust=FALSE, val=FALSE){

    ## Costa Rican spotted rose snapper
    if(species=="CRSNAP"){

        ## growth - from Bystrom thesis
        vbk <- 0.21
        linf <- 64.58
        t0 <- -0.01
        CVlen <- 0.2
        lwa <- 0.0245
        lwb <- 2.790
            
        ## mortality
        M <- 0.43 ## based on vbk
        AgeMax <- 23

        ## recruitment
        R0 <- 1
        h <- 1

        ## index
        qcoef <- 1e-2 #1e-5

        ## fishing mortality
        Fequil <- 0.25
        Frate <- 0.2
        F1 <- 0.2
        Fmax <- 1

        ## variation terms
        SigmaF <- 0.3
        SigmaC <- 0.2
        SigmaI <- 0.2
        SigmaR <- 0.6

        ## selectivity
        S50 <- 3
        S95 <- 5 ## 4.43 fully recruited (39.1 cm TL) (from Bystrom thesis)
        
        ## maturity
        ML50 <- 34 ## Rojas 2006 Gulf of Nicoya

        ## length bins
        binwidth <- 1

        ## fishing mortality
        F1 <- 0.34  

        ## sensitivities
        if("linf" %in% param_adjust) linf <- val[which(param_adjust=="linf")]
        if("vbk" %in% param_adjust) vbk <- val[which(param_adjust=="vbk")]
        if("M" %in% param_adjust) M <- val[which(param_adjust=="M")]
        if("CVlen" %in% param_adjust) CVlen <- val[which(param_adjust=="CVlen")]
        if("SigmaR" %in% param_adjust) SigmaR <- val[which(param_adjust=="SigmaR")]
        if("ML50" %in% param_adjust) ML50 <- val[which(param_adjust=="ML50")]
        if("R0" %in% param_adjust) R0 <- val[which(param_adjust=="R0")]
        if("binwidth" %in% param_adjust) binwidth <- val[which(param_adjust=="binwidth")]


        ## derived
        mids <- seq((binwidth/2), linf*1.5, by=binwidth) # from 120 cm
        highs <- mids + (binwidth/2)
        lows <- mids - (binwidth)/2
        
        ages <- 0:AgeMax
        Amat <- round(t0-log(1-(ML50/linf))/vbk)
        A95 <- Amat+1
        ML95 <- round(linf*(1-exp(-vbk*(A95-t0))))
        SL50 <- round(linf*(1-exp(-vbk*(S50-t0))))
        SL95 <- round(linf*(1-exp(-vbk*(S95-t0))))

        ## growth at age
        L_a <- linf*(1-exp(-vbk*(ages - t0)))
        W_a <- lwa*L_a^lwb  

        ## maturity
        Mat_l <- 1 / (1 + exp(ML50 - mids))
        Mat_ages <- round(t0-log(1-(mids/linf))/vbk)
        names(Mat_l) <- Mat_ages
        Mat_a <- rep(NA, (AgeMax+1))
        for(a in 1:(AgeMax+1)){
            if(a==1) Mat_a[a] <- 1e-20
            if(a>1){
                fill <- Mat_l[which(names(Mat_l)==(a-1))][length(Mat_l[which(names(Mat_l)==(a-1))])]
                if(length(fill)==1) Mat_a[a] <- fill
                if(length(fill)==0) Mat_a[a] <- Mat_a[a-1]
            }
        }

        ## selectivity 
        if(selex=="asymptotic"){
            S_a <- c(1e-20, 1 / (1 + exp(S50 - ages[-1]))) # Selectivity at age
            S_a[1] <- 1e-20
            Syoung <- NA
            Sold <- NA
        }
        if(selex=="dome"){
            S_a_calc <- rep(NA, length(ages))
            Syoung <- 3
            Sold <- 15
            A <- sqrt(2/pi)/(Syoung + Sold)
            for(a in 1:length(ages)){
                if(a <= S95) S_a_calc[a] <- A*exp(-((ages[a] - S95)^2)/(2*Syoung^2))
                if(a > S95) S_a_calc[a] <- A*exp(-((ages[a] - S95)^2)/(2*Sold^2))
            }
            S_a_calc[1] <- 1e-20
            S_a <- S_a_calc/max(S_a_calc)
        }
    }

    ## Kenyan Siganus sutor (rabbitfish) ** not yet tested
    if(species=="SIGSUT"){
        ## from FishBase, Hicks and McClanahan 2012
        ## growth
        vbk <- 0.87
        linf <- 36.2
        t0 <- -0.24
        CVlen <- 0.1
        lwa <- 0.05970
        lwb <- 2.754

        ## recruitment
        R0 <- 1
        h <- 1

        ## index
        qcoef <- 1e-2

        ## fishing mortality
        Fequil <- 0.25
        Frate <- 0.2
        F1 <- 0.2
        Fmax <- 1

        ## variation terms
        SigmaF <- 0.3
        SigmaC <- 0.2
        SigmaI <- 0.2
        SigmaR <- 0.6

        ## bins
        binwidth <- 1

        ## selectivity
        ML50 <- 20.2
        SL50 <- 11.3

        ## fishing mortality
        F1 <- 1

        ## derived
        M <- 1.49

        ## sensitivities
        if("linf" %in% param_adjust) linf <- val[which(param_adjust=="linf")]
        if("vbk" %in% param_adjust) vbk <- val[which(param_adjust=="vbk")]
        if("M" %in% param_adjust) M <- val[which(param_adjust=="M")]
        if("CVlen" %in% param_adjust) CVlen <- val[which(param_adjust=="CVlen")]
        if("SigmaR" %in% param_adjust) SigmaR <- val[which(param_adjust=="SigmaR")]
        if("ML50" %in% param_adjust) ML50 <- val[which(param_adjust=="ML50")]
        if("R0" %in% param_adjust) R0 <- val[which(param_adjust=="R0")]
        if("binwidth" %in% param_adjust) binwidth <- val[which(param_adjust=="binwidth")]

        AgeMax <- round(-log(0.01)/M)
        ages <- 0:AgeMax
        Amat <- round(t0-log(1-(ML50/linf))/vbk)
        A95 <- Amat+1
        ML95 <- round(linf*(1-exp(-vbk*(A95-t0))))
        S50 <- ceiling(t0-log(1-(SL50/linf))/vbk)
        S95 <- S50+1
        SL95 <- round(linf*(1-exp(-vbk*(S95-t0))))

        mids <- seq((binwidth/2), by=binwidth, length=76) ## max 76 from the data, but set to something at least 1.2* expected Linf
        highs <- mids + (binwidth/2)
        lows <- mids - (binwidth)/2

        ## growth at age
        L_a <- linf*(1-exp(-vbk*(ages - t0)))
        W_a <- lwa*L_a^lwb  

        ## maturity
        Mat_l <- 1 / (1 + exp(ML50 - mids))
        Mat_ages <- round(t0-log(1-(mids/linf))/vbk)
        names(Mat_l) <- Mat_ages
        Mat_a <- rep(NA, (AgeMax+1))
        for(a in 1:(AgeMax+1)){
            if(a==1) Mat_a[a] <- 1e-20
            if(a>1){
                fill <- Mat_l[which(names(Mat_l)==(a-1))][length(Mat_l[which(names(Mat_l)==(a-1))])]
                if(length(fill)==1) Mat_a[a] <- fill
                if(length(fill)==0) Mat_a[a] <- Mat_a[a-1]
            }
        }

        ## selectivity 
        if(selex=="asymptotic"){
            S_l <- 1/(1+exp(-log(19)*(mids-SL50)/(SL95-SL50)))
            S_ages <- round(t0-log(1-(mids/linf))/vbk)
            names(S_l) <- S_ages
            S_a <- rep(NA, (AgeMax+1))
            for(a in 1:(AgeMax+1)){
                if(a==1) S_a[a] <- 1e-20
                if(a>1){
                    fill <- S_l[which(names(S_l)==(a-1))][length(S_l[which(names(S_l)==(a-1))])]
                    if(length(fill)==1) S_a[a] <- fill
                    if(length(fill)==0) S_a[a] <- S_a[a-1]
                }
            }
            # S_a <- c(1e-20, 1/(1+exp(-log(19)*(ages[-1]-S50)/(S95-S50)))) # Selectivity at age
            Syoung <- NA
            Sold <- NA
        }
        if(selex=="dome"){
            S_a_calc <- rep(NA, length(ages))
            Syoung <- 0.8
            Sold <- 8
            A <- sqrt(2/pi)/(Syoung + Sold)
            for(a in 1:length(ages)){
                if(a <= S95) S_a_calc[a] <- A*exp(-((ages[a] - S95)^2)/(2*Syoung^2))
                if(a > S95) S_a_calc[a] <- A*exp(-((ages[a] - S95)^2)/(2*Sold^2))
            }
            S_a_calc[1] <- 1e-20
            S_a <- S_a_calc/max(S_a_calc)
        }

    }

    ## Namibian hake
    if(species=="HAKE"){
        ## growth
        vbk <- 0.2
        linf <- 100
        t0 <- -0.01
        CVlen <- 0.2
        lwa <- 5e-06
        lwb <- 3

        ## recruitment
        R0 <- 1
        h <- 1

        ## index
        qcoef <- 1e-2

        ## fishing mortality
        Fequil <- 0.25
        Frate <- 0.2
        F1 <- 0.2
        Fmax <- 1

        ## variation terms
        SigmaF <- 0.3
        SigmaC <- 0.2
        SigmaI <- 0.2
        SigmaR <- 0.6

        ## bins
        binwidth <- 1

        ## selectivity
        ML50 <- 20.2
        S50 <- 2
        SL50 <- round(t0-log(1-(2/linf))/vbk)

        ## fishing mortality
        F1 <- 0.25

        ## derived
        M <- 0.15

        ## sensitivities
        if("linf" %in% param_adjust) linf <- val[which(param_adjust=="linf")]
        if("vbk" %in% param_adjust) vbk <- val[which(param_adjust=="vbk")]
        if("M" %in% param_adjust) M <- val[which(param_adjust=="M")]
        if("CVlen" %in% param_adjust) CVlen <- val[which(param_adjust=="CVlen")]
        if("SigmaR" %in% param_adjust) SigmaR <- val[which(param_adjust=="SigmaR")]
        if("ML50" %in% param_adjust) ML50 <- val[which(param_adjust=="ML50")]
        if("R0" %in% param_adjust) R0 <- val[which(param_adjust=="R0")]
        if("binwidth" %in% param_adjust) binwidth <- val[which(param_adjust=="binwidth")]

        AgeMax <- round(-log(0.01)/M)
        ages <- 0:AgeMax
        Amat <- round(t0-log(1-(ML50/linf))/vbk)
        A95 <- Amat+1
        ML95 <- round(linf*(1-exp(-vbk*(A95-t0))))
        S95 <- S50+1
        SL95 <- round(linf*(1-exp(-vbk*(S95-t0))))

        mids <- seq((binwidth/2), by=binwidth, length=76) ## max 76 from the data, but set to something at least 1.2* expected Linf
        highs <- mids + (binwidth/2)
        lows <- mids - (binwidth)/2

        ## growth at age
        L_a <- linf*(1-exp(-vbk*(ages - t0)))
        W_a <- lwa*L_a^lwb  

        ## maturity
        Mat_a <- c(1e-20, 1 / (1 + exp(Amat - ages[-1])))

        ## selectivity 
        if(selex=="asymptotic"){
            S_a <- c(1e-20, 1 / (1 + exp(S50 - ages[-1])) ) # Selectivity at age
            Syoung <- NA
            Sold <- NA
        }
        if(selex=="dome"){
            S_a_calc <- rep(NA, length(ages))
            Syoung <- 0.8
            Sold <- 8
            A <- sqrt(2/pi)/(Syoung + Sold)
            for(a in 1:length(ages)){
                if(a <= S95) S_a_calc[a] <- A*exp(-((ages[a] - S95)^2)/(2*Syoung^2))
                if(a > S95) S_a_calc[a] <- A*exp(-((ages[a] - S95)^2)/(2*Sold^2))
            }
            S_a_calc[1] <- 1e-20
            S_a <- S_a_calc/max(S_a_calc)
        }

    }


    ## output list
    Outs <- NULL
    Outs$vbk <- vbk
    Outs$linf <- linf
    Outs$t0 <- t0
    Outs$binwidth <- binwidth
    Outs$CVlen <- CVlen
    Outs$SigmaC <- SigmaC
    Outs$SigmaI <- SigmaI
    Outs$SigmaR <- SigmaR ## starting value only, will be estimated
    Outs$R0 <- R0
    Outs$lwa <- lwa
    Outs$lwb <- lwb
    Outs$S50 <- S50
    Outs$S95 <- S95
    Outs$SL50 <- SL50
    Outs$SL95 <- SL95
    Outs$Syoung <- Syoung
    Outs$Sold <- Sold
    Outs$h <- h
    Outs$qcoef <- qcoef
    Outs$SigmaF <- SigmaF
    Outs$M <- M
    Outs$F1 <- F1
    Outs$AgeMax <- AgeMax
    Outs$mids <- mids
    Outs$highs <- highs
    Outs$lows <- lows
    Outs$S_a <- S_a
    Outs$L_a <- L_a
    Outs$W_a <- W_a
    Outs$Amat <- Amat
    Outs$ML50 <- ML50
    Outs$ML95 <- ML95
    Outs$Mat_a <- Mat_a
    Outs$Fequil <- Fequil
    Outs$Frate <- Frate
    Outs$F1 <- F1
    Outs$Fmax <- Fmax

    return(Outs)

}

#' Create input parameters for TMB model
#'
#' \code{create_inputs} Gets list of parameter inputs into the proper format
#'
#' @param param parameter name to adjust (sensitivity analysis)
#' @param val value of parameter name to adjust (sensitivity analysis)
#' @param lh_list tagged list of life history/starting value information
#' @param data_avail_list artifact from sensitivity analysis, adjusts some information on sample size, years, etc. for different data-availability scenarios

#' @return List, a tagged list of potentially useful benchmarks
#' @details Specifically used to merge life history information with other model settings; flexibility to change parameter inputs for sensitivity analysis without changing the baseline life history information that was used to generate data in a simulation study, or carefully compiled for real=life application
#' @export
create_inputs <- function(param, val, lh_list, data_avail_list){
    
        ## copy life history
        dat_input <- c(lh_list, data_avail_list)

        ## have the log ready in the input file for some variance parameterss
        dat_input$log_sigma_C <- log(lh_list$SigmaC)
        dat_input$log_sigma_I <- log(lh_list$SigmaI)

        ## change input values for sensitivity analysis
        if(param!=FALSE){
            for(pp in 1:length(param)){
                dat_input[[param[pp]]] <- val[which(param==param[pp])]
            }
        }
        
        dat_input$log_CV_L <- log(dat_input$CVlen)

    return(dat_input)
}

#' Create hypothetical life histories for simulation model
#'
#' \code{create_lh_list} Specifies 4 distinct life history types for simulation testing
#'
#' @param lh choose from different life histories, 1-4
#' @param param_adjust possibility of adjusting true parameter values - names of parameters, default=FALSE
#' @param val possibility of adjusting true parameter values - values for parameter names in "param_adjust" (must specify val if param_adjust is specified), default=FALSE
#' @param selex "asymptotic" selectivity or "dome"-shaped selectivity as the true selectivity in the operating model
#' @param nlbins specified number of length bins so that all life histories match (useful in simulation study)

#' @return List, a tagged list of life history traits
#' @details Life histories used from Hordyk et al. 2015 ICES Journal of Marine Science, simulation test of LB-SPR method
#' @export
create_lh_list <- function(lh, param_adjust=FALSE, val=FALSE, selex, nlbins=50){

    ## sand sole P. melanostictus
    if(lh==1){
        ## growth
        vbk <- 0.79
        linf <- 37.6
        t0 <- -0.01
        CVlen <- 0.1
        lwa <- 0.00912
        lwb <- 3.09

        ## natural mortality
        M <- 0.42

        ## recruitment
        R0 <- 1
        h <- 1

        ## fishing mortality
        Fequil <- 0.25
        Frate <- 0.2
        F1 <- 0.2
        Fmax <- 1

        ## index
        qcoef <- 1e-2

        ## deviations
        SigmaF <- 0.3
        SigmaC <- 0.2
        SigmaI <- 0.2
        SigmaR <- 0.6
        

        ## bins
        binwidth <- 1
        
        ## selectivity
        SL50 <- 24.0
        SL95 <- 26.0

        ## maturity
        ML50 <- 29.0
        ML95 <- 32.0
    }

    ## Puget Sound rockfish (S. emphaeus)
    if(lh==2){
        ## growth
        vbk <- 0.535
        linf <- 17.0
        t0 <- -0.01
        CVlen <- 0.1
        lwa <- 0.01259
        lwb <- 3.08

        ## natural mortality
        M <- 0.44

        ## recruitment
        R0 <- 1
        h <- 1

        ## fishing mortality
        Fequil <- 0.25
        Frate <- 0.2
        F1 <- 0.2
        Fmax <- 1

        ## index
        qcoef <- 1e-2

        ## deviations
        SigmaF <- 0.3
        SigmaC <- 0.2
        SigmaI <- 0.2
        SigmaR <- 0.6
        

        ## bins
        binwidth <- 1
        
        ## selectivity
        SL50 <- 9.4
        SL95 <- 10.8

        ## maturity
        ML50 <- 12.1
        ML95 <- 17.0
    }

    ## yellowtail flathead (P. endrachtensis)
    if(lh==3){
        ## growth
        vbk <- 0.41
        linf <- 53.0
        t0 <- -0.01
        CVlen <- 0.1
        lwa <- 0.00490
        lwb <- 3.06

        ## natural mortality
        M <- 0.63

        ## recruitment
        R0 <- 1
        h <- 1

        ## fishing mortality
        Fequil <- 0.25
        Frate <- 0.2
        F1 <- 0.2
        Fmax <- 1

        ## index
        qcoef <- 1e-2

        ## deviations
        SigmaF <- 0.3
        SigmaC <- 0.2
        SigmaI <- 0.2
        SigmaR <- 0.6
        

        ## bins
        binwidth <- 1
        
        ## selectivity
        SL50 <- 22.0
        SL95 <- 26.0

        ## maturity
        ML50 <- 25.9
        ML95 <- 34.4
    }

    ## Pacific saury (C. saira)
    if(lh==4){
        ## growth
        vbk <- 0.41
        linf <- 34.2
        t0 <- -0.01
        CVlen <- 0.1
        lwa <- 0.00240
        lwb <- 3.15

        ## natural mortality
        M <- 1.25

        ## recruitment
        R0 <- 1
        h <- 1

        ## fishing mortality
        Fequil <- 0.25
        Frate <- 0.2
        F1 <- 0.2
        Fmax <- 1

        ## index
        qcoef <- 1e-2

        ## deviations
        SigmaF <- 0.3
        SigmaC <- 0.2
        SigmaI <- 0.2
        SigmaR <- 0.6
        

        ## bins
        binwidth <- 1
        
        ## selectivity
        SL50 <- 13.0
        SL95 <- 14.5

        ## maturity
        ML50 <- 19.4
        ML95 <- 20.4
    }

    if(lh==5){
      ## growth - from Bystrom thesis
        vbk <- 0.21
        linf <- 64.58
        t0 <- -0.01
        CVlen <- 0.2
        lwa <- 0.0245
        lwb <- 2.790
            
        ## mortality
        M <- 0.43 ## based on vbk

        ## recruitment
        R0 <- 1
        h <- 1

        ## fishing mortality
        Fequil <- 0.25
        Frate <- 0.2
        F1 <- 0.2
        Fmax <- 1

        ## index
        qcoef <- 1e-2

        ## variation terms
        SigmaF <- 0.3
        SigmaC <- 0.2
        SigmaI <- 0.2
        SigmaR <- 0.6

        ## bins
        binwidth <- 1

        ## selectivity
        SL50 <- 30.0
        SL95 <- 39.1
        
        ## maturity
        ML50 <- 34 ## Rojas 2006 Gulf of Nicoya
        ML95 <- 40
    }
    
    ## sensitivities
    if("linf" %in% param_adjust) linf <- val[which(param_adjust=="linf")]
    if("vbk" %in% param_adjust) vbk <- val[which(param_adjust=="vbk")]
    if("M" %in% param_adjust) M <- val[which(param_adjust=="M")]
    if("CVlen" %in% param_adjust) CVlen <- val[which(param_adjust=="CVlen")]
    if("SigmaR" %in% param_adjust) SigmaR <- val[which(param_adjust=="SigmaR")]
    if("ML50" %in% param_adjust) ML50 <- val[which(param_adjust=="ML50")]

    ## derived
    AgeMax <- round(-log(0.01)/M)
    ages <- 0:AgeMax
    Amat <- round(t0-log(1-(ML50/linf))/vbk)
    A95 <- Amat+1
    S50 <- round(t0-log(1-(SL50/linf))/vbk)
    S95 <- round(t0-log(1-(SL95/linf))/vbk)
    if(S50==S95) S95 <- S50 + 1

    mids <- seq((binwidth/2), nlbins, by=binwidth)
    highs <- mids + (binwidth/2)
    lows <- mids - (binwidth)/2

    ## growth at age
    L_a <- linf*(1-exp(-vbk*(ages - t0)))
    W_a <- lwa*L_a^lwb

    ## maturity
    Mat_l <- 1 / (1 + exp(ML50 - mids))
    Mat_ages <- round(t0-log(1-(mids/linf))/vbk)
    names(Mat_l) <- Mat_ages
    Mat_a <- rep(NA, (AgeMax+1))
    for(a in 1:(AgeMax+1)){
        if(a==1) Mat_a[a] <- 1e-20
        if(a>1){
            fill <- Mat_l[which(names(Mat_l)==(a-1))][length(Mat_l[which(names(Mat_l)==(a-1))])]
            if(length(fill)==1) Mat_a[a] <- fill
            if(length(fill)==0) Mat_a[a] <- Mat_a[a-1]
        }
    }
    # Mat_a <- c(1e-20, 1 / (1 + exp(Amat - ages[-1])))

    ## selectivity 
    if(selex=="asymptotic"){
        S_l <- 1 / (1 + exp(SL50 - mids))
        S_ages <- round(t0-log(1-(mids/linf))/vbk)
        names(S_l) <- S_ages
        S_a <- rep(NA, (AgeMax+1))
        for(a in 1:(AgeMax+1)){
            if(a==1) S_a[a] <- 1e-20
            if(a>1){
                fill <- S_l[which(names(S_l)==(a-1))][length(S_l[which(names(S_l)==(a-1))])]
                if(length(fill)==1) S_a[a] <- fill
                if(length(fill)==0) S_a[a] <- S_a[a-1]
            }
        }
        # S_a <- c(1e-20, 1/(1+exp(-log(19)*(ages[-1]-S50)/(S95-S50)))) # Selectivity at age
        Syoung <- NA
        Sold <- NA
    }
    if(selex=="dome"){
        S_a_calc <- rep(NA, length(ages))
        Syoung <- S50
        Sold <- max(c(S50+1),AgeMax*0.75)
        A <- sqrt(2/pi)/(Syoung + Sold)
        for(a in 1:length(ages)){
            if(a==1) S_a_calc[a] <- 1e-20
            if(a <= S95) S_a_calc[a] <- A*exp(-((ages[a] - S95)^2)/(2*Syoung^2))
            if(a > S95) S_a_calc[a] <- A*exp(-((ages[a] - S95)^2)/(2*Sold^2))
        }
        S_a <- S_a_calc/max(S_a_calc)
    }

    ## output list
    Outs <- NULL
    Outs$vbk <- vbk
    Outs$linf <- linf
    Outs$t0 <- t0
    Outs$binwidth <- binwidth
    Outs$CVlen <- CVlen
    Outs$SigmaC <- SigmaC
    Outs$SigmaI <- SigmaI
    Outs$SigmaR <- SigmaR
    Outs$R0 <- R0
    Outs$lwa <- lwa
    Outs$lwb <- lwb
    Outs$S50 <- S50
    Outs$S95 <- S95
    Outs$SL50 <- SL50
    Outs$SL95 <- SL95
    Outs$Syoung <- Syoung
    Outs$Sold <- Sold
    Outs$h <- h
    Outs$Fequil <- Fequil
    Outs$Frate <- Frate
    Outs$F1 <- F1
    Outs$qcoef <- qcoef
    Outs$Fmax <- Fmax
    Outs$SigmaF <- SigmaF
    Outs$M <- M
    Outs$AgeMax <- AgeMax
    Outs$mids <- mids
    Outs$highs <- highs
    Outs$lows <- lows
    Outs$S_a <- S_a
    Outs$L_a <- L_a
    Outs$W_a <- W_a
    Outs$Amat <- Amat
    Outs$ML50 <- ML50
    Outs$ML95 <- ML95
    Outs$Mat_a <- Mat_a

    return(Outs)

}

#' Settings for simulation testing of data availability scenarios
#'
#' \code{data_avail_settings} Specifies 4 distinct life history types for simulation testing
#'
#' @param avail_set names of data availability scenarios, must match names in function, see details
#' @param ESS Effective sample size, default 1000
#' @param simulation default=TRUE, not set up for simulation=FALSE 

#' @return List, a tagged list of model settings
#' @details avail_set: "Rich_LC"= 20 years length comp, index, and catch, high ESS; "Moderate_LC"=20 years length comp, index, and catch, lower ESS;  "Sample_LC"=20 years length comp, index, and catch, but only sampled every few years, and catch is reported at a rate of 20 percent; "Index_LC1"=20 years of abundance index, no catch data, 1 year of Length comp (can also specify "Index_LC10" for 10 years of length comp); "Catch_LC10"=20 years of catch data, 1 year of length comp (can also specify "Catch_LC10" for 10 years of length comp); "LC1" is 1 year of length comp only; can also specify LC2, 5, 10, and 20.
#' @export
data_avail_settings <- function(avail_set, ESS, simulation=TRUE){
        
    ### simulation
    if(simulation==TRUE){
        settings <- list()
        Nyears <- 20
        if("Rich_LC" %in% avail_set){
            settings$Rich_LC <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=Nyears, alt_yrs=FALSE, sample=FALSE)
        }
        if("Moderate_LC" %in% avail_set){
            settings$Moderate_LC <- list(Nyears=Nyears, comp_sample=50, obs_per_yr=rep(ESS,Nyears), Nyears_comp=Nyears, alt_yrs=FALSE, sample=FALSE)
        }
        if("Sample_LC" %in% avail_set){
            settings$Sample_LC <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=Nyears, alt_yrs=TRUE, sample=0.2)
        }
        if("Index_LC1" %in% avail_set){
            settings$Index_LC1 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=1, alt_yrs=FALSE, sample=FALSE)
        }
        if("Catch_LC1" %in% avail_set){
            settings$Catch_LC1 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=1, alt_yrs=FALSE, sample=FALSE)
        }
        if("Index_LC10" %in% avail_set){
            settings$Index_LC10 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=10, alt_yrs=FALSE, sample=FALSE)
        }
        if("Catch_LC10" %in% avail_set){
            settings$Catch_LC10 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=10, alt_yrs=FALSE, sample=FALSE)
        }

        ### must have at least 1 year of composition data specified to test method with 0 years length comp - just don't include it in data supplied to model
        if("Index_LC0" %in% avail_set){
            settings$Index_LC0 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=1, alt_yrs=FALSE, sample=FALSE)
        }
        if("Catch_LC0" %in% avail_set){
            settings$Catch_LC0 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=1, alt_yrs=FALSE, sample=FALSE)
        }

        if("LC1" %in% avail_set){
            settings$LC1 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=1, alt_yrs=FALSE, sample=FALSE)
        }
        if("LC2" %in% avail_set){
            settings$LC2 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=2, alt_yrs=FALSE, sample=FALSE)
        }
        if("LC5" %in% avail_set){
            settings$LC5 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=5, alt_yrs=FALSE, sample=FALSE)
        }
        if("LC10" %in% avail_set){
            settings$LC10 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=10, alt_yrs=FALSE, sample=FALSE)
        }
        if("LC20" %in% avail_set){
            settings$LC20 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=20, alt_yrs=FALSE, sample=FALSE)
        }
    }

    return(settings)

}

#' Formatting real data for specific species for TMB model
#'
#' \code{formatDataList} Reads data files from directory for specific species and gets the data into a form to input to the TMB model
#'
#' @param species species code, "SIGSUT" or "CRSNAP"
#' @param data_dir Directory to look for data files 

#' @return List, a tagged list of data to input into TMB model
#' @details required output: I_t: index with each element named by year 1-x, C_t: catch with each element named 1-x, LF: length frequency with years 1-x labeled on the rows and length bin labeled on the columns, LFprop: proportions in each length bin, same dimensions as LF, years: actual years of data, years_i: index years 1-x, lbins: length bins, ML_t: mean length with each element named year 1-x, Nyears: number of years, Nyears_comp: number of years o f length composition data, obs_per_yr: effective sample size of length composition annually
#' @export
formatDataList <- function(species, data_dir){

    if(species=="SIGSUT"){
        sp <- "Siganus_sutor"
        data <- read.csv(file.path(data_dir, paste0(sp, "_data.csv")), header=TRUE)
        meltDat <- melt(data, id.vars=c("Date", "Year", "Landing", "Management.type", "Group", 
            "Fishgear","X.Fishers", "Catch.category", "Family", "Name", "Catch.category",
            "Commercial.noncomm."))
        meltDat$value <- as.numeric(meltDat$value)
        dat_new <- meltDat[-which(is.na(meltDat$value)),-which(colnames(meltDat)=="variable")]      
    
        years <- unique(dat_new$Year)[order(unique(dat_new$Year))]
        tyears <- min(years):max(years)
        tyears_i <- 1:length(tyears)
        years_i <- tyears_i[which(tyears %in% years)]
        lbins <- 1:(max(dat_new$value)*1.5)     

        lcomp1 <- matrix(0, nrow=length(years), ncol=length(lbins))
        colnames(lcomp1) <- lbins
        rownames(lcomp1) <- years_i   
        lcomp_p <- lcomp1   

        catch <- cpue <- vector(length=length(years))
        cpue_day_out <- list()      

        obs_high <- obs_low <- rep(0, length(years))
        for(y in 1:length(years)){
            sub <- dat_new[which(dat_new$Year==years[y]),]   
            dates <- unique(as.character(sub$Date))  
            obs_high[y] <- length(dates)
            months <- unique(as.numeric(sapply(1:length(dates), function(x) strsplit(as.character(dates[x]), "/")[[1]][2])))
            obs_low[y] <- length(months)

            ## annual catch
            catch[y] <- nrow(sub)      

            ## annual cpue - average of daily/group cpue
            cpue_sub <- rep(NA, length(unique(sub$Date)))
            ngroups <- rep(NA, length(unique(sub$Date)))
            for(d in 1:length(unique(sub$Date))){
                subd <- sub[which(sub$Date==unique(sub$Date)[d]),]    
                ngroups <- length(unique(subd$Group)[which(is.na(unique(subd$Group))==FALSE)])
                if(ngroups==0) cpue_sub[d] <- NA
                if(ngroups!=0) cpue_sub[d] <- nrow(subd)/ngroups
                # rm(ngroups)
                # rm(subd)
            }
            if(all(is.na(cpue_sub))) cpue[y] <- NA
            if(all(is.na(cpue_sub)==FALSE)) cpue[y] <- mean(cpue_sub, na.rm=TRUE)
            cpue_day_out[[y]] <- cpue_sub
            rm(cpue_sub)       

            for(l in 1:length(lbins)){
                lcomp1[y,l] <- length(which(floor(sub$value)==lbins[l]))
            }
        }       

        names(cpue) <- names(catch) <- years_i
        cpue <- cpue[-c(which(is.na(cpue)),which(cpue==0))]     

        meanlen <- rep(NA, nrow(lcomp1))
        names(meanlen) <- rownames(lcomp1)
        for(i in 1:nrow(lcomp1)){
            lcomp_p[i,] <- lcomp1[i,]/sum(lcomp1[i,])
            meanlen[i] <- sum(sapply(1:ncol(lcomp1), function(x) lcomp1[i,x]*as.numeric(colnames(lcomp1)[x])))/sum(lcomp1[i,])
        }  

        DataList <- NULL
        DataList$I_t <- cpue
        DataList$C_t <- catch
        DataList$LF <- lcomp1
        DataList$LFprop <- lcomp_p
        DataList$years <- years
        DataList$years_i <- years_i
        DataList$lbins <- lbins
        DataList$meanlen <- meanlen
        DataList$Nyears <- length(min(years):max(years))
        DataList$Nyears_comp <- nrow(LFprop)
        DataList$obs_per_year <- obs_high
    }

    if(species=="CRSNAP"){

        lg <<- read.csv(file.path(data_dir, "cr_snapper_filtered.csv"))

        ## subset by gears
        lg_bl <- lg[which(lg$Gear=="Bottom Longline"),]
        lg_g <- lg[which(lg$Gear=="Gillnet"),]      

        ## annual mean length
        ml <- mean_length(data=lg, plot=FALSE)      

        ## life history info
        cr_lh <- choose_lh_list(species="CRSNAP", selex="asymptotic")       

        ## length frequency data
        lf <- length_frequency(binwidth=1, linf=cr_lh$linf, lmat=cr_lh$L50, data=lg, plot=FALSE, weight=TRUE)       

        ## catch and effort data
        fishery_data <- catch_effort(data=lg, sep_index=TRUE)
        catch <- fishery_data$catch
        cpue_bl <- fishery_data$cpue_bl
        cpue_g <- fishery_data$cpue_g
        tyears <- c((fishery_data$years[1]-10):(fishery_data$years[1]-1),fishery_data$years)
        tyears_i <- 1:length(tyears)
        years <- fishery_data$years
        years_i <- tyears_i[which(tyears %in% years)]
        obs_high <- obs_low <- rep(0, length(tyears_i))
        obs_high[which(tyears_i %in% years_i)] <- fishery_data$obs_high ## number of days fished per year 
        obs_low[which(tyears_i %in% years_i)] <- fishery_data$obs_low ## number of months fished per year

        ## bottom longline cpue index
        cpue_input <- cpue_bl[which(is.na(cpue_bl)==FALSE)] ## choose bottom longline cpue only
        cpue_yrs <- which(tyears %in% names(cpue_bl)[which(is.na(cpue_bl)==FALSE)])
        names(cpue_input) <- cpue_yrs       

        ## length frequency
        lf_input <- lf[which(rowSums(lf)>0),]
        lf_yrs_i <- tyears_i[which(tyears %in% names(which(rowSums(lf)>0)))]
        rownames(lf_input) <- lf_yrs_i      

        lcomp_p <- t(apply(lf_input, 1, function(x) x/sum(x)))
        rownames(lf_input) <- lf_yrs_i      

        lbins <- 1:ncol(lcomp_p)        

        ## mean length
        meanlen_input <- ml$all_gears
        meanlen_yrs_i <- tyears_i[which(tyears %in% names(meanlen_input))]
        names(meanlen_input) <- meanlen_yrs_i

                DataList <- NULL
                DataList$I_t <- cpue_input
                DataList$C_t <- NULL
                DataList$LF <- lf_input
                DataList$LFprop <- lcomp_p
                DataList$years <- tyears
                DataList$years_i <- tyears_i
                DataList$lbins <- lbins
                DataList$ML_t <- meanlen_input
                DataList$Nyears <- length(tyears_i)
                DataList$Nyears_comp <- nrow(lcomp_p)
                DataList$obs_per_year <- obs_high


    }
    return(DataList)

}


#' TMB input formatting
#'
#' \code{FormatInput_LB} Formats data, parameters, random effects, and mapped parameters for TMB input
#'
#' @param Nyears total number of years to model
#' @param DataList list of data required for model
#' @param linf fixed life history parameter from input list
#' @param vbk fixed life history parameter from input list
#' @param t0 fixed life history parameter from input list
#' @param M fixed life history parameter from input list
#' @param AgeMax fixed life history parameter from input list
#' @param lbhighs fixed length bins from input list
#' @param lbmids fixed length bins from input list
#' @param Mat_a fixed life history parameter from input list
#' @param lwa fixed life history parameter from input list
#' @param lwb fixed life history parameter from input list
#' @param log_sigma_C starting value from input list
#' @param log_sigma_I starting value from input list
#' @param log_CV_L starting value from input list
#' @param F1 starting value from input list
#' @param SigmaR starting value from input list
#' @param qcoef starting value from input list
#' @param R0 starting value from input list
#' @param S50 starting value from input list
#' @param model data type availability
#' @param RecDev_biasadj starting values for rec devs
#' @param site artifact of development: how many sites to estimate values?
#' @param Fpen penalty on fishing mortality 0= off, 1=on
#' @param Dpen penalty on terminal depletion 0= off, 1=on
#' @param Dprior prior info for depletion penalty
#' @param SigRpen penalty on sigmaR, 0=off, 1=on
#' @param SigRprior prior info for sigmaR penalty
#' @param obs_per_yr effective sample sizes for length comp (for each of total years)
#' @param SigmaF starting value from input list
#' @param RecType artifact of development: 0=random effects on recruitment
#' @param FType artifact of development: 0=fishing mortality fixed effects
#' @param LType artifact of development: 1= pooled length comp between sites
#' @param h starting value from input list
#' @param SelexTypeDesc artifact of development - asymptotic or dome-shaped selectivity
#' @param est_sigma list of variance parameters to estimate, must match parameter names
#' @param REML default FALSE
#' @param estimate_same TRUE=estimate least-common-denominator parameters for all data availability scenarios, FALSE=estimate parameters specific to data availability
#' @param start_f year (in numbers, not actual year) to start estimating fishing mortality (e.g. year 11 out of 20 to get estimates for last 10 years); the value of F in this year will be used as the estimate and SE for all previous years. 0=estimate all years.

#' @return List, a tagged list of Data, Parameters, Random, Map
#' @export
FormatInput_LB <- function(Nyears, DataList, linf, vbk, t0, M, AgeMax,
    lbhighs, lbmids, Mat_a, lwa, lwb, log_sigma_C, log_sigma_I, log_CV_L, F1, SigmaR, 
    qcoef, R0, S50, model, RecDev_biasadj, site,
    Fpen, Dpen, Dprior, SigRpen, SigRprior, obs_per_yr, SigmaF, RecType, FType, LType, h, SelexTypeDesc, est_sigma, REML, estimate_same, start_f){

        ## Rich, Moderate, and Sample names include catch, index, and length comp data - just vary in the sampling and ESS of composition data
        ## fit to length composition data
        if((grepl("Rich",model) | grepl("Moderate", model) | grepl("Sample", model)) & grepl("LC", model)){
            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=length(DataList$C_t),
                n_i=length(DataList$I_t), 
                n_lc=nrow(DataList$LF),
                n_ml=0, start_f=start_f,
                T_yrs=1:Nyears, C_yrs=as.numeric(names(DataList$C_t)),
                I_yrs=as.numeric(names(DataList$I_t)),
                LC_yrs=as.numeric(rownames(DataList$LF)),
                ML_yrs=as.vector(0),
                rel_c=0, rel_i=0,
                obs_per_yr=obs_per_yr, RecType=RecType,
                I_t=DataList$I_t, C_t=DataList$C_t, 
                ML_t=as.vector(0), LF=DataList$LF,
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen, SigRpen=SigRpen, SigRprior=SigRprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }

        ## same as above but fit to mean length data instead of length composition data
        if((grepl("Rich",model) | grepl("Moderate", model) | grepl("Sample", model)) & grepl("ML", model)){
            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=length(DataList$C_t),
                n_i=length(DataList$I_t), 
                n_lc=0,
                n_ml=nrow(DataList$LF), start_f=start_f,
                T_yrs=1:Nyears, C_yrs=as.numeric(names(DataList$C_t)),
                I_yrs=as.numeric(names(DataList$I_t)),
                LC_yrs=as.vector(0),
                ML_yrs=as.numeric(rownames(DataList$LF)),
                rel_c=0, rel_i=0,
                obs_per_yr=obs_per_yr, RecType=RecType,
                I_t=DataList$I_t, C_t=DataList$C_t, 
                ML_t=rowMeans(DataList$LF), LF=as.matrix(0),
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen, SigRpen=SigRpen, SigRprior=SigRprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }

        ## index and length composition data
        if(grepl("Index", model) & grepl("LC", model) & grepl("LC0", model)==FALSE){
            if(is.matrix(DataList$LF)){
                n_lc <- nrow(DataList$LF)
                LC_yrs <- as.numeric(rownames(DataList$LF))
                LF <- as.matrix(DataList$LF)
            }
            if(is.vector(DataList$LF)){
                n_lc <- 1
                LC_yrs <- Nyears
                LF <- t(as.matrix(DataList$LF))
            }
            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=0,
                n_i=length(DataList$I_t), 
                n_lc=n_lc,
                n_ml=0, start_f=start_f,
                T_yrs=1:Nyears, C_yrs=as.vector(0),
                I_yrs=as.numeric(names(DataList$I_t)),
                LC_yrs=LC_yrs,
                ML_yrs=as.vector(0),
                rel_c=0, rel_i=0, 
                obs_per_yr=obs_per_yr, RecType=RecType,
                I_t=DataList$I_t, C_t=as.vector(0), 
                ML_t=as.vector(0), LF=LF,
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen, SigRpen=SigRpen, SigRprior=SigRprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }
        ## index only
        if(grepl("Index", model) & grepl("LC", model) & grepl("LC0", model)){
            n_lc <- 0
            LC_yrs <- as.vector(0)
            LF <- as.matrix(0)
            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=0,
                n_i=length(DataList$I_t), 
                n_lc=n_lc,
                n_ml=0, start_f=start_f,
                T_yrs=1:Nyears, C_yrs=as.vector(0),
                I_yrs=as.numeric(names(DataList$I_t)),
                LC_yrs=LC_yrs,
                ML_yrs=as.vector(0),
                rel_c=0, rel_i=0, 
                obs_per_yr=obs_per_yr, RecType=RecType,
                I_t=DataList$I_t, C_t=as.vector(0), 
                ML_t=as.vector(0), LF=LF,
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen, SigRpen=SigRpen, SigRprior=SigRprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }

        ## index and mean length data
        if(grepl("Index", model) & grepl("ML", model)){
            if(is.matrix(DataList$LF)){
                n_ml <- nrow(DataList$LF)
                ML_yrs <- as.numeric(rownames(DataList$LF))
                LF <- as.matrix(DataList$LF)
            }
            if(is.vector(DataList$LF)){
                n_ml <- 1
                ML_yrs <- Nyears
                LF <- t(as.matrix(DataList$LF))
            }
            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=0,
                n_i=length(DataList$I_t), 
                n_lc=0,
                n_ml=n_ml, start_f=start_f,
                T_yrs=1:Nyears, C_yrs=as.vector(0),
                I_yrs=as.numeric(names(DataList$I_t)),
                LC_yrs=as.vector(0),
                ML_yrs=ML_yrs,
                rel_c=0, rel_i=0,
                obs_per_yr=obs_per_yr, RecType=RecType,
                I_t=DataList$I_t, C_t=as.vector(0), 
                ML_t=rowMeans(LF), LF=as.matrix(0),
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen,  SigRpen=SigRpen, SigRprior=SigRprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }

        ## catch and length composition data
        if(grepl("Catch", model) & grepl("LC", model) & grepl("LC0", model)==FALSE){
            if(is.matrix(DataList$LF)){
                n_lc <- nrow(DataList$LF)
                LC_yrs <- as.numeric(rownames(DataList$LF))
                LF <- as.matrix(DataList$LF)
            }
            if(is.vector(DataList$LF)){
                n_lc <- 1
                LC_yrs <- Nyears
                LF <- t(as.matrix(DataList$LF))
            }
            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=length(DataList$C_t),
                n_i=0, 
                n_lc=n_lc,
                n_ml=0, start_f=start_f,
                T_yrs=1:Nyears, C_yrs=as.numeric(names(DataList$C_t)),
                I_yrs=as.vector(0),
                LC_yrs=LC_yrs,
                ML_yrs=as.vector(0),
                rel_c=0, rel_i=0, 
                obs_per_yr=obs_per_yr, RecType=RecType,
                I_t=as.vector(0), C_t=DataList$C_t, 
                ML_t=as.vector(0), LF=LF,
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen, SigRpen=SigRpen, SigRprior=SigRprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }

        ## catch only
        if(grepl("Catch", model) & grepl("LC", model) & grepl("LC0", model)){
            n_lc <- 0
            LC_yrs <- as.vector(0)
            LF <- as.matrix(0)
            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=length(DataList$C_t),
                n_i=0, 
                n_lc=n_lc,
                n_ml=0, start_f=start_f,
                T_yrs=1:Nyears, C_yrs=as.numeric(names(DataList$C_t)),
                I_yrs=as.vector(0),
                LC_yrs=LC_yrs,
                ML_yrs=as.vector(0),
                rel_c=0, rel_i=0, 
                obs_per_yr=obs_per_yr, RecType=RecType,
                I_t=as.vector(0), C_t=DataList$C_t, 
                ML_t=as.vector(0), LF=LF,
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen, SigRpen=SigRpen, SigRprior=SigRprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }


        ## catch and mean length data
        if(grepl("Catch", model) & grepl("ML", model)){
            if(is.matrix(DataList$LF)){
                n_ml <- nrow(DataList$LF)
                ML_yrs <- as.numeric(rownames(DataList$LF))
                LF <- as.matrix(DataList$LF)
            }
            if(is.vector(DataList$LF)){
                n_ml <- 1
                ML_yrs <- Nyears
                LF <- t(as.matrix(DataList$LF))
            }
            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=length(DataList$C_t),
                n_i=0, 
                n_lc=0,
                n_ml=n_ml, start_f=start_f,
                T_yrs=1:Nyears, C_yrs=as.numeric(names(DataList$C_t)),
                I_yrs=s.vector(0),
                LC_yrs=as.vector(0),
                ML_yrs=ML_yrs,
                rel_c=0, rel_i=0,
                obs_per_yr=obs_per_yr, RecType=RecType,
                I_t=as.vector(0), C_t=DataList$C_t, 
                ML_t=rowMeans(LF), LF=as.matrix(0),
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen,  SigRpen=SigRpen, SigRprior=SigRprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }

        ## length composition data only 
        if(grepl("LC", model) & grepl("Index", model)==FALSE & grepl("Catch", model)==FALSE & grepl("Rich", model)==FALSE & grepl("Moderate", model)==FALSE & grepl("Sample", model)==FALSE){
            if(is.matrix(DataList$LF)){
                n_lc <- nrow(DataList$LF)
                LC_yrs <- as.numeric(rownames(DataList$LF))
                LF <- as.matrix(DataList$LF)
            }
            if(is.vector(DataList$LF)){
                n_lc <- 1
                LC_yrs <- Nyears
                LF <- t(as.matrix(DataList$LF))
            }
            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=0,
                n_i=0, 
                n_lc=n_lc,
                n_ml=0, start_f=start_f,
                T_yrs=1:Nyears, C_yrs=as.vector(0),
                I_yrs=as.vector(0),
                LC_yrs=LC_yrs,
                ML_yrs=as.vector(0),
                rel_c=0, rel_i=0,
                obs_per_yr=obs_per_yr, RecType=RecType,
                I_t=as.vector(0), C_t=as.vector(0), 
                ML_t=as.vector(0), LF=LF,
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen,  SigRpen=SigRpen, SigRprior=SigRprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }

        ## set input parameters - regardless of data availability        
        Parameters <- list(log_F_sd=log(SigmaF), log_F_t_input=log(rep(F1,Nyears)),log_q_I=log(qcoef), beta=log(R0), log_sigma_R=log(SigmaR), logS50=log(S50), log_sigma_C=log_sigma_C, log_sigma_I=log_sigma_I, log_CV_L=log_CV_L,Nu_input=rep(0,Nyears))

        ## turn off parameter estimation - depends on data availability
            Map = list()

           
            if(estimate_same==TRUE){
                Map[["beta"]] <- NA
                Map[["beta"]] <- factor(Map[["beta"]])
                Map[["log_q_I"]] <- NA
                Map[["log_q_I"]] <- factor(Map[["log_q_I"]])
            }

            if("log_F_sd" %in% est_sigma==FALSE){
                Map[["log_F_sd"]] <- NA
                Map[["log_F_sd"]] <- factor(Map[["log_F_sd"]])
            }

            if("log_sigma_R" %in% est_sigma==FALSE){
                Map[["log_sigma_R"]] <- NA
                Map[["log_sigma_R"]] <- factor(Map[["log_sigma_R"]])
            }

            if("log_sigma_C" %in% est_sigma==FALSE){
                Map[["log_sigma_C"]] <- NA
                Map[["log_sigma_C"]] <- factor(Map[["log_sigma_C"]])
            }
            if("log_sigma_I" %in% est_sigma==FALSE){
                Map[["log_sigma_I"]] <- NA
                Map[["log_sigma_I"]] <- factor(Map[["log_sigma_I"]])
            }
            if("log_CV_L" %in% est_sigma==FALSE){
                Map[["log_CV_L"]] <- NA
                Map[["log_CV_L"]] <- factor(Map[["log_CV_L"]])
            }
            if(all(est_sigma==FALSE)){
                Map[["log_CV_L"]] <- NA
                Map[["log_CV_L"]] <- factor(Map[["log_CV_L"]]) 
                Map[["log_sigma_C"]] <- NA
                Map[["log_sigma_C"]] <- factor(Map[["log_sigma_C"]])  
                Map[["log_sigma_I"]] <- NA
                Map[["log_sigma_I"]] <- factor(Map[["log_sigma_I"]])
                Map[["log_sigma_R"]] <- NA
                Map[["log_sigma_R"]] <- factor(Map[["log_sigma_R"]])

            }

            if(start_f>0){
                Map[["log_F_t_input"]] = 1:length(Parameters[["log_F_t_input"]])
                Map[["log_F_t_input"]][1:start_f] <- NA
                Map[["log_F_t_input"]] <- factor(Map[["log_F_t_input"]])   
            }

            if(estimate_same==FALSE){
                if(grepl("Index", model)){        
                    Map[["beta"]] <- NA
                    Map[["beta"]] <- factor(Map[["beta"]])
                }

                if(grepl("Catch", model)){
                    Map[["log_q_I"]] <- NA
                    Map[["log_q_I"]] <- factor(Map[["log_q_I"]])
                }

                if(grepl("LC", model) & grepl("Index", model)==FALSE & grepl("Catch", model)==FALSE & grepl("Rich", model)==FALSE & grepl("Moderate", model)==FALSE & grepl("Sample", model)==FALSE){                        
                   Map[["log_q_I"]] <- NA
                   Map[["log_q_I"]] <- factor(Map[["log_q_I"]])       
                   Map[["beta"]] <- NA
                   Map[["beta"]] <- factor(Map[["beta"]])
                }
            }

        if(length(Map)==0) Map <- NULL


        if(REML==FALSE) Random <- c("Nu_input")
        if(REML==TRUE){
            # Random_vec <- c("Nu_input", "log_F_t_input", "log_q_I", "beta", "logS50", "logS95") # 
            Random_vec <- c("Nu_input", "log_F_t_input", "log_q_I", "beta", "logS50") # 
            Random <- Random_vec[which(Random_vec %in% names(Map) == FALSE)]
        }
        if("log_F_sd" %in% est_sigma) Random <- c(Random, "log_F_t_input")




    Return <- list("Parameters"=Parameters, "Data"=Data, "Random"=Random, "Map"=Map)
    return(Return)

}

#' Generate data for simulation testing
#'
#' \code{generateData} Generates data from the operating model for use in simulation testing
#'
#' @param modpath directory to save generated data
#' @param modname name of model (to identify differences between different model runs)
#' @param itervec number of iterations of data to generate
#' @param spatial does asymptotic length vary by spatially? TRUE or FALSE
#' @param Fdynamics Specify name of pattern of fishing mortality dynamics, "Constant", "Endogenous", "Ramped", "Increasing", or "None"
#' @param Rdynamics Specify name of pattern of recruitment dynamics, "Constant", "Pulsed", "Pulsed_up", or "BH"
#' @param LType 1 (default): pool site-specific length composition into 1 dataset; 0: keep length composition collected at different sites separate by site
#' @param plotML plot of generated mean length data, default=FALSE
#' @param plotLF_compare plot length comp pooled vs not pooled, default=FALSE
#' @param plotLF plot length composition at chosen LType, default=FALSE
#' @param selex "asymptotic" for asympotitc selectivity in length composition generation (currently only option programmed)
#' @param write write generated dataset? default=TRUE. FALSE helpful for plotting.
#' @param lh_list list of life history inputs
#' @param data_avail_list list of other model settings
#' @param param_adjust vector of names of parameters to adjust true value, default FALSE to include no parameter
#' @param val vector of values aligning with names in param_adjust to adjust value to, default FALSE to include no parameter
#' @param rewrite TRUE will re-run OM and observation model. FALSE will skip if it's already written in directory.

#' @return print how many iterations were written into the model directory
#' @export
generateData <- function(modpath, modname, itervec, spatial, Fdynamics, Rdynamics, LType=1, plotML=FALSE, plotLF_compare=FALSE, plotLF=FALSE, selex="asymptotic", write=TRUE, lh_list, data_avail_list, param_adjust=FALSE, val=FALSE, rewrite){

    lh_num <- ifelse(grepl("LH1", modpath), 1, ifelse(grepl("LH2", modpath), 2, ifelse(grepl("LH3", modpath), 3, ifelse(grepl("LH4", modpath), 4, ifelse(grepl("LH5", modpath), 5, ifelse(grepl("CRSNAP", modpath), "CRSNAP", ifelse(grepl("SIGSUT", modpath), "SIGSUT", ifelse(grepl("HAKE", modpath), "HAKE", stop("No match to life history number")))))))))
  lh_choose <- lh_list[[lh_num]]
  if(param_adjust[1]!=FALSE) lh_choose[param_adjust] <- val
  
  Nyears_comp <- data_avail_list$Nyears_comp
  Nyears <- data_avail_list$Nyears

  for(iter in itervec){

    iterpath <- file.path(modpath, iter)
    if(write==TRUE) dir.create(iterpath, showWarnings=FALSE) 
    if(rewrite==FALSE){
      if(file.exists(file.path(iterpath, "True.rds"))) next
    }

    ## simulated data with no spatial structure in growth
    DataList <- with(c(lh_choose, data_avail_list), SimData_LB(Nyears=Nyears, AgeMax=AgeMax, 
      M=M, F1=F1, h=h, S_a=S_a, qcoef=qcoef, Frate=Frate, Fequil=Fequil, 
      SigmaF=SigmaF, Fdynamics=Fdynamics, Rdynamics=Rdynamics, 
      R0=R0, Fmax=Fmax, CVlen=CVlen, mids=mids, highs=highs, 
      lows=lows, W_a=W_a, L_a=L_a, Mat_a=Mat_a, Amat=Amat, 
      comp_sample=comp_sample, SigmaR=SigmaR, Nyears_comp=Nyears_comp, 
      alt_yrs=FALSE, sample=FALSE, nburn=20, seed=iter, modname=modname)) 

    # simulated data with spatial structure
    if(spatial==TRUE){

      set.seed(max(itervec)+iter)
      ## spatial process in Linf - varies with each iteration
      spatial_sim <- spatialgrowth_sim(n_i=20, linf=lh_choose$linf)
      if(write==TRUE) saveRDS(spatial_sim, file.path(iterpath, "spatial_sim.rds"))  

      ## life history - truth with spatial structure - varies with each iteration
      lh_spatial <- lapply(1:nrow(spatial_sim), function(x) choose_lh_list(species=lh_num, param_adjust=c("linf","ML50"), val=c(spatial_sim[x,"linf_i"], lh_choose$ML50*(spatial_sim[x,"linf_i"]/lh_choose$linf)), selex="asymptotic")) 

      ## simulated data with spatial structure in growth
      DataList_site <- lapply(1:length(lh_spatial), function(x) with(c(lh_spatial[[x]], data_avail_list), SimData_LB(Nyears=Nyears, AgeMax=AgeMax,
          M=M, F1=F1, h=h, S_a=S_a, qcoef=qcoef, Frate=Frate, Fequil=Fequil, 
          SigmaF=SigmaF, Fdynamics=Fdynamics, Rdynamics=Rdynamics,
          R0=R0, Fmax=Fmax, CVlen=CVlen, mids=mids, highs=highs,
          lows=lows, W_a=W_a, L_a=L_a, Mat_a=Mat_a, Amat=Amat,
          comp_sample=comp_sample, SigmaR=SigmaR, Nyears_comp=Nyears_comp,
          alt_yrs=FALSE, sample=FALSE, nburn=20, seed=iter, modname=modname)))  
      SPR_site <- sapply(1:length(DataList_site), function(x) DataList_site[[x]]$SPR_t)
      RelAbund_site <- sapply(1:length(DataList_site), function(x) DataList_site[[x]]$D_t[length(DataList_site[[x]]$D_t)])
      if(write==TRUE) saveRDS(SPR_site, file.path(iterpath, "SPR_site.rds"))

      ## length frequency data at each site
      if(Nyears_comp>1) LF_site <- lapply(1:length(DataList_site), function(x) DataList_site[[x]]$LF)
      if(Nyears_comp==1) LF_site <- lapply(1:length(DataList_site), function(x) as.matrix(DataList_site[[x]]$LF))
      ncols_site <- sapply(1:length(LF_site), function(x) ncol(LF_site[[x]]))
      for(i in 1:length(LF_site)){
        ncol <- ncol(LF_site[[i]])
        if(ncol < max(ncols_site)){
          add <- max(ncols_site) - ncol
          LF_site[[i]] <- cbind(LF_site[[i]], matrix(0, nrow=nrow(LF_site[[i]]), ncol=add))
        }
      } 

      ## length frequency by site
      LF_site_array <- array(NA, dim=c(dim(LF_site[[1]]), length(LF_site)))
      ML_t_site <- matrix(NA, nrow=nrow(LF_site[[1]]), ncol=length(LF_site))
      for(i in 1:length(LF_site)){
        LF_site_array[,,i] <- as.matrix(LF_site[[i]])
        if(Nyears_comp>1) ML_t_site[,i] <- sapply(1:nrow(LF_site_array[,,i]), function(x) sum(LF_site[[i]][x,]*1:ncol(LF_site[[i]]))/sum(LF_site[[i]][x,]))
        if(Nyears_comp==1) ML_t_site[,i] <- sum(LF_site[[i]]*1:length(LF_site[[i]]))/sum(LF_site[[i]])
      }
      rownames(LF_site_array) <- rownames(ML_t_site) <- (Nyears-Nyears_comp+1):Nyears #

      ## length frequency pooled across sites
      LF_pool <- matrix(NA, nrow=nrow(LF_site[[1]]), ncol=max(ncols_site))
      for(i in 1:nrow(LF_site[[1]])){
        for(j in 1:ncol(LF_site[[1]])){
          LF_pool[i,j] <- sum(sapply(1:length(LF_site), function(x) LF_site[[x]][i,j]))
        }
      }
      rownames(LF_pool) <- (Nyears-Nyears_comp+1):Nyears

    if(plotLF==TRUE){
      ## length frequency in the last year at each site
      par(mfrow=c(4,4), mar=c(0,0,0,0), omi=c(1,1,1,1))
      barplot(DataList$LF[nrow(DataList$LF),]/sum(DataList$LF[nrow(DataList$LF),]), axes=F, xlim=c(0,45), ylim=c(0,0.2), col="tomato")
      mtext(side=3, "no spatial process", font=2, line=-3, cex=2)
      axis(2, at=pretty(c(0,0.2)))      
      for(i in 1:length(LF_site)){
        barplot(LF_site[[i]][nrow(LF_site[[i]]),]/sum(LF_site[[i]][nrow(LF_site[[i]]),]), axes=F, xlim=c(0,45), ylim=c(0,0.2))
        mtext(paste0("site ", i), side=3, font=2, line=-3, cex=2)
        if(i %in% 12:15) axis(1, at=pretty(c(0,45)))
        if(i %% 4==0) axis(2, at=pretty(c(0,0.2)))
      }
      mtext(side=1, "Length bin (1 cm)", outer=TRUE, line=3)
      mtext(side=2, "Proportion", outer=TRUE, line=3)
    }
    if(plotLF_compare==TRUE){
      par(mfrow=c(2,1), mar=c(0,0,0,0), omi=c(1,1,1,1))
      barplot(DataList$LF[nrow(DataList$LF),]/sum(DataList$LF[nrow(DataList$LF),]), axes=F, xlim=c(0,45), ylim=c(0,0.2), col="tomato")
      mtext(side=3, "no spatial process", font=2, line=-3, cex=2)
      axis(2, at=pretty(c(0,0.2)))      

      barplot(LF_pool[nrow(LF_pool),]/sum(LF_pool[nrow(LF_pool),]), axes=F, xlim=c(0,45), ylim=c(0,0.2), col="black")
      mtext(side=3, "spatial process pooled", font=2, line=-3, cex=2)
      axis(1, at=pretty(c(0,45)))
      axis(2, at=pretty(c(0,0.2)))
      mtext(side=1, "Length bin (1 cm)", outer=TRUE, line=3)
      mtext(side=2, "Proportion", outer=TRUE, line=3)      
    }


      ## mean length over time pooled across sites
      ML_t_pool <- sapply(1:nrow(LF_pool), function(x) sum(LF_pool[x,]*1:ncol(LF_pool))/sum(LF_pool[x,])) 

      if(plotML==TRUE){
      # ## plot mean length at each site
      # # png("SIM_Mean_length_by_site_year.png", res=200, units="in", width=12, height=9)
      par(mfrow=c(4,4), mar=c(0,0,0,0), omi=c(1,1,1,1))
      plot(ML_t_pool, col="red", type="o", pch=17, lwd=2, ylim=c(0,linf), xaxt="n", yaxt="n")
      axis(2, cex=1.2, las=2)
      mtext(side=3, "pooled", font=2, line=-1.5)
      for(i in 1:length(DataList_site)){
        plot(DataList_site[[i]]$ML_t, col="black", pch=19, lwd=2, ylim=c(0,linf), xaxt="n", yaxt="n")
        lines(ML_t_pool, col="red", pch=17, type="o")
        if(i %in% c(12:15)) axis(1, cex=1.2)
        if(i %in% c(4,8,12)) axis(2, cex=1.2, las=2)
        mtext(paste0("site ", i), side=3, font=2, line=-1.5)
      }
      legend("bottomright", legend=c("site-specific", "pooled"), pch=c(19,17), col=c("black", "red"))
      mtext("Year", outer=TRUE, line=3, side=1)
      mtext("Mean length in catch (cm)", outer=TRUE, line=3,  side=2)
      # # dev.off()
    }
      if(LType==1){
        DataList$LF <- LF_pool
        DataList$ML_t <- ML_t_pool  
      }
      if(LType==0){
        DataList$LF <- LF_site_array
        DataList$ML_t <- ML_t_site
      }


      rm(spatial_sim)
      rm(lh_spatial)
    }
   
    DataList_out <- DataList
    DataList_out$LF <- DataList$LF[,1:length(lh_choose$mids)]

    Inputs <- FormatInput_LB(Nyears=Nyears, DataList=DataList_out, linf=lh_choose$linf, vbk=lh_choose$vbk, t0=lh_choose$t0, M=lh_choose$M, AgeMax=lh_choose$AgeMax, lbhighs=lh_choose$highs, lbmids=lh_choose$mids, Mat_a=lh_choose$Mat_a, lwa=lh_choose$lwa, lwb=lh_choose$lwb, log_sigma_C=log(lh_choose$SigmaC), log_sigma_I=log(0.001), log_CV_L=log(0.001), F1=DataList$F_t[1], SigmaR=lh_choose$SigmaR, qcoef=lh_choose$qcoef, R0=1, S50=lh_choose$S50, model="Rich_LC", RecDev_biasadj=rep(0,Nyears), Fpen=1, Dpen=0, Dprior=c(0,0), SigRpen=1, SigRprior=c(lh_choose$SigmaR, 0.2), obs_per_yr=rep(1000,Nyears), SigmaF=lh_choose$SigmaF, RecType=0, FType=0, LType=1, h=lh_choose$h, SelexTypeDesc="asymptotic", est_sigma="log_sigma_R", REML=FALSE, site=1, estimate_same=FALSE, start_f=0)
    ParList <- Inputs$Parameters
    obj <- MakeADFun(data=Inputs[["Data"]], parameters=ParList, random=Inputs[["Random"]], map=Inputs[["Map"]],inner.control=list(maxit=1e3), hessian=FALSE, DLL="LIME")  
    Derived <- Calc_derived_quants(obj)

    DataList_out$MSY <- Derived$MSY
    DataList_out$Fmsy <- Derived$Fmsy

      if(write==TRUE) saveRDS(DataList_out, file.path(iterpath, "True.rds"))
      if(write==FALSE) return(DataList_out)
      rm(DataList)
      rm(DataList_out)
      rm(iterpath)

}

  if(write==TRUE) return(paste0(length(itervec), " iterates of data generated in ", modpath))

}


#' Set up model paths
#'
#' \code{model_paths} sets up directory paths for model combinations
#'
#' @param res_dir directory to store results
#' @param modcombos data frame of all model combinations

#' @return vector of directory paths
#' @export
model_paths <- function(res_dir, modcombos){

    devo_path <- function(combo, res_dir){
        old <- res_dir
        for(i in 1:length(combo)){
            new <- file.path(old, combo[i])
            old <- new
            dir.create(old, showWarnings=FALSE)
        }
        return(old)
    }

    alldirs <- sapply(1:nrow(modcombos), function(x) devo_path(combo=modcombos[x,], res_dir=res_dir))

    return(alldirs)
}


#' Plot relative error
#'
#' \code{plotRE} boxplots comparing relative error across several models
#'
#' @param modpath_vec directory to find model results
#' @param itervec number of iterations
#' @param modnames model names for x axis
#' @param xaxt default "n" to remove x axis, follows same rules as plot function
#' @param yaxt default "n" to remove y axis, follows same rules as plot function
#' @param value relative error for specified reference point, "SPR" or "FFref" (associated with F/F40)
#' @param ylim set y axis limits
#' @param col.plot color for boxplots
#' @param col.line color for line designating 0 relative error
#' @param yr year to plot relative error


#' @return vector of directory paths
#' @export
plotRE <- function(modpath_vec, itervec, modnames=NULL, xaxt="n", yaxt="n", value, ylim=c(-1,1.5), col.plot="steelblue", col.line="goldenrod", yr){

    RE <- t(sapply(itervec, function(x) calcRE(modpath_vec=modpath_vec, iter=x, value=value, yr=yr)))

    colnames(RE) <- modnames

    boxplot(RE, ylim=ylim, xaxt=xaxt, yaxt=yaxt, col=col.plot)
    abline(h=0, col=col.line, lwd=2)
    
    return(RE)
}


#' model fits or kobe plots
#'
#' \code{plotResults} plot results for time series of estimated/derived parameters or kobe plots
#'
#' @param Data read RDS from model directory
#' @param Report read RDS from model directory
#' @param Sdreport read RDS from model directory
#' @param Derived_quants read RDS from model directory
#' @param flag_convergence TRUE- flag in directory that model didn't converge, FALSE = no flag in directory
#' @param parameter which parameter to plot results for
#' @param xaxt plot xaxis, TRUE or FALSE
#' @param ylab plot ylab, TRUE or FALSE
#' @param simulation get different items if these are application or simulation results

#' @return displays plot
#' 
#' @details possible values for parameter argument for model fits: "B" biomass, "N" abundance, "ML" mean length, "R" recruitment, "F" fishing mortality, "D" relative biomass, "C" catch, "I" abundance index, for kobe plot, "kobe" will show SPR compared with F/F30 and F/F40
#' @export
plotResults <- function(Data, Report, Sdreport, Derived_quants, flag_convergence, parameter, xaxt=TRUE, ylab=FALSE, simulation=TRUE){
        
        FUN = function(InputMat, log=TRUE){
          if(log==TRUE) return(c( exp(InputMat[,1]-1.96*InputMat[,2]), rev(exp(InputMat[,1]+1.96*InputMat[,2]))))
          if(log==FALSE) return(c( InputMat[,1]-1.96*InputMat[,2], rev(InputMat[,1]+1.96*InputMat[,2])))
        } 

      if(simulation==TRUE) years <- 1:length(Data$SB_t)
      if(simulation==FALSE){
        years <- Data$years_i
        years_real <- Data$years
      }
          if(parameter=="B"){
            ## Biomass
            if(simulation==TRUE) Mat <- cbind("Year"=years, "Data"=Data$SB_t, "Est"=Report$SB_t_hat)
            if(simulation==FALSE) Mat <- cbind("Year"=years, "Data"=rep(NA, length(years)), "Est"=Report$SB_t_hat)
            ymax <- 3
            matplot(y=Mat[,c("Data", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE) if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lSB_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            mtext("Biomass", side=3, line=-2, font=2)
            if(ylab==TRUE) mtext("Biomass", side=2, line=3, font=2)
          }
          if(parameter=="N"){
            ## Abundance
            if(simulation==TRUE) Mat <- cbind("Year"=years, "Data"=Data$N_t, "Est"=Report$N_t_hat)
            if(simulation==FALSE) Mat <- cbind("Year"=years, "Data"=rep(NA, length(years)), "Est"=Report$N_t_hat)
            ymax <- 4
            matplot(y=Mat[,c("Data", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lN_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            mtext("Abundance", side=3, line=-2, font=2)
            if(ylab==TRUE) mtext("Abundance", side=2, line=3, font=2)

          }
          if(parameter=="ML"){
            ## Average Length
            tml <- Data$ML_t
            pml <- rep(NA, length(years))
              names(pml) <- years
            pml[which(names(pml) %in% names(Data$ML_t))] <- tml
            Mat <- cbind("Year"=years, "Data"=pml, "Est"=Report$L_t_hat)
            ymax <- max(Report$L_t_hat)*1
            matplot(y=Mat[,c("Data", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="L_t_hat"),], log=FALSE), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            mtext("Mean Length", side=3, line=-2, font=2)
            if(ylab==TRUE) mtext("Mean Length", side=2, line=3, font=2)
          }
          if(parameter=="R"){       
            ## Recruitment
            if(simulation==TRUE) Mat <- cbind("Year"=years, "Data"=Data$R_t, "Est"=Report$R_t_hat)
            if(simulation==FALSE) Mat <- cbind("Year"=years, "Data"=rep(NA, length(years)), "Est"=Report$R_t_hat)
            ymax <- 3
            matplot(y=Mat[,c("Data", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE) if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            mtext("Recruitment", side=3, line=-2, font=2)
            if(ylab==TRUE) mtext("Recruitment", side=2, line=3, font=2)

          }       
          if(parameter=="F"){
            ## Fishing Mortality
            NAs <- rep(NA, length(years))
            if(simulation==TRUE) Mat <- cbind("Year"=years, "Data"=Data$F_t, "Est"=Report$F_t_hat)
            if(simulation==FALSE) Mat <- cbind("Year"=years, "Data"=NAs, "Est"=Report$F_t)
            ymax <- 2
            matplot(y=Mat[,c("Data", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            mtext("Fishing Mortality", side=3, line=-2, font=2)
            if(ylab==TRUE) mtext("Fishing Mortality", side=2, line=3, font=2)

          }
          if(parameter=="D"){     
            ## Relative abundance
            if(simulation==TRUE) Mat <- cbind("Year"=years, "Data"=Data$D_t, "Est"=Report$Depl)
            if(simulation==FALSE) Mat <- cbind("Year"=years, "Data"=rep(NA, length(years)), "Est"=Report$Depl)
            ymax <- 1
            matplot(y=Mat[,c("Data", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lD_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            mtext("Relative abundance", side=3, line=-2, font=2)  
            if(ylab==TRUE) mtext("Relative abundance", side=2, line=3, font=2)
          }
          if(parameter=="C"){      
            ## Catch
            tcatch <- Data$C_t
            pcatch <- rep(NA, length(years))
              names(pcatch) <- years
            pcatch[which(names(pcatch) %in% names(Data$C_t))] <- tcatch
            Mat <- cbind("Year"=years, "Data"=pcatch, "Est"=Report$C_t_hat)
            ymax <- 5
            matplot(y=Mat[,c("Data", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lC_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            mtext("Catch", side=3, line=-2, font=2)   
            if(ylab==TRUE) mtext("Catch", side=2, line=3, font=2)
          }
          if(parameter=="I"){     
            ## Index
            tindex <- Data$I_t
            pindex <- rep(NA, length(years))
              names(pindex) <- years
            pindex[which(names(pindex) %in% names(Data$I_t))] <- tindex
            Mat <- cbind("Year"=years, "Data"=pindex, "Est"=Report$I_t_hat)
            ymax <- 3
            matplot(y=Mat[,c("Data", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lI_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            mtext("Index", side=3, line=-2, font=2)
            if(ylab==TRUE) mtext("Index", side=2, line=3, font=2)
        }
        if(xaxt==TRUE){
          axis(1, at=seq(1,20, by=5), labels=years_real[seq(1,20,by=5)])
          mtext(side=1, "Year", line=4)  
        }

        if(parameter=="kobe"){
          plot(x=1, y=1, type="n", xaxs="i", yaxs="i", xlim=c(0, 1), ylim=c(0, 3), xaxt="n", yaxt="n")
          abline(h=1, lty=2, lwd=2)
          abline(v=0.3, lty=2, col="forestgreen", lwd=2)
          abline(v=0.4, lty=2, col="steelblue", lwd=2)
          points(x=Derived_quants$SPR, y=Derived_quants$FF30, col="forestgreen", pch=19, cex=2)
          points(x=Derived_quants$SPR, y=Derived_quants$FF40, col="steelblue", pch=19, cex=2)
          axis(1, at=seq(0, 1, by=0.2))
          axis(2, at=seq(0, 3, by=0.2))
          mtext(side=1, "SPR", line=2)
          mtext(side=2, "F/Fref", line=2)
          legend("topright", col=c("forestgreen", "steelblue"), legend=c("30%", "40%"), title="Target SPR", pch=19)
        }
        if(flag==TRUE) mtext(side=3, "model not converged", col="red", font=2, outer=TRUE, line=2, cex=2)
  

}


#' run LIME model - simulation mode
#'
#' \code{runModel} run length-based integrated mixed-effects model with generated data
#'
#' @param modpath model directory
#' @param itervec number of iterations of data to generate
#' @param data_avail other settings for data availability
#' @param est_sigma which variance parameters to estimate, match parameter names
#' @param biascorrect bias correction for recruitment deviatiosn on (TRUE) or off (FALSE)
#' @param sensitivity_inputs (in development) named list (parameters) with matrix of 2 rows (low, high) and # of columns for 'life histories' (artifact of testing multiple life histories in the simulation, would only have 1 column for an assessment)
#' @param sensitivity_ESS (in development) will do sensitivity analysis for effective sample size of length composition
#' @param REML default off (FALSE)
#' @param estimate_same TRUE=estimate least-common-denominator parameters for all data availability scenarios, FALSE=estimate parameters specific to data availability
#' @param lh_list list of life history information
#' @param rewrite if results already exist in the directory, should we rewrite them? TRUE or FALSE
#' @param start_f year (in numbers, not actual year) to start estimating fishing mortality (e.g. year 11 out of 20 to get estimates for last 10 years); the value of F in this year will be used as the estimate and SE for all previous years. 0=estimate all years.
#' @param simulation is this a simulation? default TRUE, FALSE means you are using real data (no need for iterations or multiple life history inputs)
#' @param input_data use this to input data for a real-world application (not simulation)
#' @useDynLib LIME

#' @return prints how many iterations were run in model directory
#' 
#' @details need to adjust to run with real data
#' @export
runModel <- function(modpath, itervec, estimate_same=FALSE, REML=FALSE, est_sigma, biascorrect=TRUE, data_avail, lh_list, sensitivity_inputs=NULL, sensitivity_ESS=NULL, rewrite, start_f, simulation=TRUE, input_data=NULL){

  if(simulation==TRUE){
    lh_num <- ifelse(grepl("LH1", modpath), 1, ifelse(grepl("LH2", modpath), 2, ifelse(grepl("LH3", modpath), 3, ifelse(grepl("LH4", modpath), 4, ifelse(grepl("LH5", modpath), 5, ifelse(grepl("CRSNAP", modpath), "CRSNAP", ifelse(grepl("SIGSUT", modpath), "SIGSUT", ifelse(grepl("HAKE", modpath), "HAKE", stop("No match to life history number")))))))))
    lh_choose <- lh_list[[lh_num]]
  }
  if(simulation==FALSE){
    lh_choose <- lh_list
  }

  if(simulation==FALSE) itervec <- 1

for(iter in itervec){

    if(simulation==TRUE) iterpath <- file.path(modpath, iter)
    if(simulation==FALSE) iterpath <- modpath

    if(rewrite==FALSE){
      if(file.exists(file.path(iterpath, "Derived_quants.rds"))) next
      # if(any(grepl("LBSPR_results", list.files(path=iterpath)))) next
      if(file.exists(file.path(iterpath, "NAs_final_gradient.txt"))) next
      if(file.exists(file.path(iterpath, "high_final_gradient.txt"))) next
    }

    if(simulation==TRUE){
        DataList <- readRDS(file.path(iterpath, "True.rds"))
        modname <- DataList$DataScenario 
    }

    if(simulation==TRUE) if(grepl("MixedEffects", modname)) modname <- strsplit(modname, "_")[[1]][2]
    if(simulation==FALSE) modname <- data_avail

    ## copies life history information with any adjustments for sensitivity analyses
    if(is.null(sensitivity_inputs)){
      param <- FALSE
      val <- FALSE
    }
    if(is.null(sensitivity_inputs)==FALSE){
      param_set <- c("M", "linf", "vbk", "CVlen", "SigmaR")
      param_set_input <- paste0("sens_", param_set)
      param <- param_set[which(sapply(1:length(param_set), function(x) grepl(param_set_input[x], modpath)))]
      val_index <- ifelse(grepl("Low", modpath), 1, ifelse(grepl("High", modpath), 2, stop("Not set up for specified level of sensitivity")))
      if(simulation==TRUE) val <- as.numeric(sensitivity_inputs[[param]][val_index, lh_num])
      if(simulation==FALSE) val <- as.numeric(sensitivity_inputs[[param]][val_index])
    }
    if(simulation==TRUE) inits <- create_inputs(lh_list=lh_choose, data_avail_list=data_avail[[modname]], param=param, val=val)
    if(simulation==FALSE) inits <- create_inputs(lh_list=lh_choose, data_avail_list=input_data, param=param, val=val) 
    Nyears <- inits$Nyears 

    if(simulation==TRUE) obs_per_yr <- inits$obs_per_yr
    if(simulation==FALSE) obs_per_yr <- input_data$obs_per_year

    if(simulation==FALSE) DataList <- input_data

  if(grepl("LBSPR", modpath)==FALSE){
    
    if(biascorrect==FALSE) vec <- 1
    if(biascorrect==TRUE) vec <- 1:2
    Sdreport <- NA
    ParList <- NA  
    df <- NULL

    for(bb in vec){
      if(all(is.na(Sdreport))) RecDev_biasadj <- rep(0, Nyears)
      if(all(is.na(Sdreport))==FALSE){
          SD <- summary(Sdreport)
          RecDev_biasadj <- 1 - SD[which(rownames(SD)=="Nu_input"), "Std. Error"]^2 / Report$sigma_R^2    
      }
      if(all(is.na(RecDev_biasadj))) RecDev_biasadj <- rep(0, Nyears)
      TmbList <- FormatInput_LB(Nyears=Nyears, DataList=DataList, linf=inits$linf, vbk=inits$vbk, t0=inits$t0, M=inits$M, AgeMax=inits$AgeMax, lbhighs=inits$highs, lbmids=inits$mids, Mat_a=inits$Mat_a, lwa=inits$lwa, lwb=inits$lwb, log_sigma_C=inits$log_sigma_C, log_sigma_I=inits$log_sigma_I, log_CV_L=inits$log_CV_L, F1=inits$F1, SigmaR=inits$SigmaR, qcoef=inits$qcoef, R0=inits$R0, S50=inits$S50, model=as.character(modname), RecDev_biasadj=RecDev_biasadj,SigmaF=inits$SigmaF, Fpen=1, Dpen=0, Dprior=c(0.8, 0.2), SigRpen=1, SigRprior=c(inits$SigmaR, 0.2), obs_per_yr=obs_per_yr, RecType=0, FType=0, LType=1, h=inits$h, SelexTypeDesc="asymptotic", est_sigma=est_sigma, REML=REML, site=1, estimate_same=FALSE, start_f=start_f)
      if(bb==1) saveRDS(TmbList, file.path(iterpath, "Inputs1.rds")) 
      if(bb==2) saveRDS(TmbList, file.path(iterpath, "Inputs2.rds"))  

      # dyn.load(paste0(run_exe, "\\", dynlib(version)))     

      if(all(is.na(ParList))) ParList <- TmbList[["Parameters"]]  

      ## create objects to save best results
      if(bb==1){
        obj_save <- NULL
        jnll <- NULL
        opt_save <- NULL
        opt_save[["final_gradient"]] <- NA
      }

      ## first run
      obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList, random=TmbList[["Random"]], map=TmbList[["Map"]],inner.control=list(maxit=1e3), hessian=FALSE, DLL="LIME")  
        # if(bb==1) check_id <- Check_Identifiable2(obj)[[4]]
        # fix_f <- grep("Bad", check_id[which(check_id[,"Param"]=="log_F_t_input"),3])      
        # good_f <- c(1:Nyears)[which(1:Nyears %in% fix_f == FALSE)] 
        # TmbList$Map[["log_F_t_input"]] = 1:length(TmbList$Parameters[["log_F_t_input"]])
        # TmbList$Map[["log_F_t_input"]][fix_f] <- NA
        # TmbList$Map[["log_F_t_input"]] <- factor(TmbList$Map[["log_F_t_input"]])
        # if(length(fix_f)>0){
        #   TmbList$Data$fix_f <- fix_f
        #   TmbList$Data$fill_f <- good_f[length(good_f)]
        # }


      ## Settings
      obj$env$inner.control <- c(obj$env$inner.control, "step.tol"=c(1e-8,1e-12)[1], "tol10"=c(1e-6,1e-8)[1], "grad.tol"=c(1e-8,1e-12)[1]) 
        Upr = rep(Inf, length(obj$par))
        Upr[match("log_sigma_R",names(obj$par))] = log(2)
        # Upr[match("logS95", names(obj$par))] = log(inits$AgeMax)
        Upr[match("log50", names(obj$par))] = log(inits$AgeMax)
        Upr[which(names(obj$par)=="log_F_t_input")] = log(5)
        Upr[match("log_F_sd", names(obj$par))] <- log(2)
        Lwr <- rep(-Inf, length(obj$par))

        ## Run optimizer
        opt <- tryCatch( nlminb( start=obj$par, objective=obj$fn, gradient=obj$gr, upper=Upr, lower=Lwr, control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10) ), error=function(e) NA)    
        jnll <- obj$report()$jnll   
        if(all(is.na(opt))==FALSE){
          opt[["final_gradient"]] = obj$gr( opt$par ) 
          opt_save <- opt
          obj_save <- obj
          jnll_save <- obj_save$report()$jnll
        }      


        ## loop to try to get opt to run
          for(i in 1:10){
            if(all(is.na(opt))){
              obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList,
                            random=TmbList[["Random"]], map=TmbList[["Map"]], 
                            inner.control=list(maxit=1e3), hessian=FALSE, DLL="LIME")
                    obj$env$inner.control <- c(obj$env$inner.control, "step.tol"=c(1e-8,1e-12)[1], "tol10"=c(1e-6,1e-8)[1], "grad.tol"=c(1e-8,1e-12)[1]) 
                opt <- tryCatch( nlminb( start= obj$env$last.par.best[-obj$env$random] + rnorm(length(obj$par),0,0.2), 
                  objective=obj$fn, gradient=obj$gr, upper=Upr, lower=Lwr, 
                  control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10)), error=function(e) NA)
                jnll <- obj$report()$jnll
            }
            if(all(is.na(opt))==FALSE){
              opt[["final_gradient"]] = obj$gr( opt$par )       
              opt_save <- opt
              obj_save <- obj
              jnll_save <- jnll
              break
            }
          }
          

        ## if opt ran: 
        if(all(is.na(opt_save))==FALSE){  

          ## check convergence -- don't let it become NA after it has had a high final gradient
          for(i in 1:5){
            if(all(is.na(opt_save[["final_gradient"]])==FALSE)){
               if(abs(min(opt_save[["final_gradient"]]))>0.01){
                obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList,
                            random=TmbList[["Random"]], map=TmbList[["Map"]], 
                            inner.control=list(maxit=1e3), hessian=FALSE, DLL="LIME")
                      obj$env$inner.control <- c(obj$env$inner.control, "step.tol"=c(1e-8,1e-12)[1], "tol10"=c(1e-6,1e-8)[1], "grad.tol"=c(1e-8,1e-12)[1]) 
                opt <- tryCatch( nlminb( start= obj$env$last.par.best[-obj$env$random] + rnorm(length(obj$par),0,0.2), 
                  objective=obj$fn, gradient=obj$gr, upper=Upr, lower=Lwr, 
                  control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10)), error=function(e) NA)
                jnll <- obj$report()$jnll
              }
            }
              if(all(is.na(opt))==FALSE & jnll<=jnll_save){
                opt[["final_gradient"]] = obj$gr( opt$par )       
                opt_save <- opt
                obj_save <- obj
                jnll_save <- jnll
              }
            if(abs(min(opt_save[["final_gradient"]]))<=0.01) break
          }
        }
        if(all(is.na(opt_save))==FALSE)  df <- data.frame(opt_save$final_gradient, names(obj_save$par), opt_save$par)


        ## write error message in directory if opt wouldn't run
        if(bb==length(vec)){
          if(all(is.null(opt_save))) write("NAs final gradient", file.path(iterpath, "NAs_final_gradient.txt"))
          if(all(is.null(opt_save)==FALSE)) if(abs(min(opt_save[["final_gradient"]]))>0.01) write(opt_save[["final_gradient"]], file.path(iterpath, "high_final_gradient.txt"))
        }

        ParList <- obj_save$env$parList( x=obj_save$par, par=obj_save$env$last.par.best )
        
        ## Standard errors
        Report = tryCatch( obj_save$report(), error=function(x) NA)
        if(bb==length(vec)) saveRDS(Report, file.path(iterpath, "Report.rds"))  

        Sdreport = tryCatch( sdreport(obj_save), error=function(x) NA )
        if(bb==length(vec)) saveRDS(Sdreport, file.path(iterpath, "Sdreport.rds"))

        if(bb==length(vec)){
          Derived = Calc_derived_quants( Obj=obj_save )
          if(bb==length(vec)) saveRDS(Derived, file.path(iterpath, "Derived_quants.rds"))
        }
  

        # dyn.unload( paste0(run_exe,"\\", dynlib(version)) )  
    } 
        if(iter==1) write.csv(df, file.path(modpath, "df.csv"))  

        rm(Report)
        rm(Sdreport)
        rm(TmbList)
        rm(opt)
        rm(obj)
        rm(df)  
        rm(opt_save)
        rm(obj_save)
  }

}

return(paste0(max(itervec), " iterates run in ", modpath))

}

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
#' @param Amat Age at 50 percent maturity
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
    W_a, L_a, linf, vbk, Mat_a, Amat, comp_sample, Nyears_comp, alt_yrs=FALSE, sample=FALSE,
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
    SB_t <- F_t <- R_t <- rep(NA, tyears)                               
    Cn_at <- N_at <- matrix(NA, nrow=AgeMax+1, ncol=tyears)

    #####################################
    ## Fishing and recruitment dynamics
    #####################################   

    if(Fdynamics=="Ramp") Framp_t <- c(rep(0.01, nburn), "rampup"=seq(F1, Fmax, length=floor(Nyears/2)), 
        "peak"=rep(Fmax, floor((Nyears-floor(Nyears/2))/2)), 
        "managed"=rep(Fmax/3, Nyears-floor(Nyears/2)-floor((Nyears-floor(Nyears/2))/2)))
    if(Fdynamics=="Constant") Fconstant_t <- c(rep(0.01, nburn), rep(Fequil, Nyears))
    if(Fdynamics=="Increasing") Finc_t <- c(rep(0.01, nburn), seq(0.01, Fmax, length=Nyears))

    if(Rdynamics=="Pulsed") Rpulse_t <- c(rep(R0, nburn), "initial"=rep(R0, floor(Nyears/3)),
        "pulse_down"=rep(R0/3, floor(Nyears/3)), "pulse_up"=rep(R0, Nyears-floor(Nyears/3)))
    if(Rdynamics=="Pulsed_up") Rpulse_t <- c(rep(R0, nburn), "initial"=rep(R0, floor(Nyears/3)), "pulse_up"=rep(R0*3, floor(Nyears/3)), "pulse_down"=rep(R0, Nyears-floor(Nyears/3)))
    if(Rdynamics=="Constant") Rconstant_t <- rep(R0, tyears)

    ##########################
    ## Initialization
    ##########################
    if(Fdynamics=="Endogenous") F_t[1] <- 0.01
    if(Fdynamics=="Ramp") F_t[1] <- Framp_t[1]
    if(Fdynamics=="Constant") F_t[1] <- Fconstant_t[1]
    if(Fdynamics=="Increasing") F_t[1] <- Finc_t[1]
    if(Fdynamics=="None") F_t[1] <- 0

    R_t[1] <- R0

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
    SB_t[1] <- sum(N_at[,1] * W_a * S_a)
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
    SB0 <- sum(Na0[-1]*Mat_a[-1]*W_a[-1])

    for(y in 2:tyears){
        ## fishing effort and recruitment, not dependent on age structure
        if(Fdynamics=="Endogenous"){
            if(y <= nburn) F_t[y] <- 0.01
            if(y > nburn) F_t[y] <- F_t[y-1]*(SB_t[y-1]/(Fequil*SB0))^Frate * exp(FishDev[y])
        }
        if(Fdynamics=="Ramp"){
            F_t[y] <- Framp_t[y] * exp(FishDev[y])
        }
        if(Fdynamics=="Constant"){
            F_t[y] <- Fconstant_t[y] * exp(FishDev[y])
        }
        if(Fdynamics=="Increasing"){
            F_t[y] <- Finc_t[y] * exp(FishDev[y])
        }
        if(Fdynamics=="None"){
            F_t[y] <- 0
        }
        if(Rdynamics=="Constant"){
            R_t[y] <- Rconstant_t[y] * exp(RecDev[y])
        }
        if(Rdynamics=="Pulsed"){
            R_t[y] <- Rpulse_t[y] * exp(RecDev[y])
        }
        if(Rdynamics=="Pulsed_up"){
            R_t[y] <- Rpulse_t[y] * exp(RecDev[y])
        }
        if(Rdynamics=="BH"){
            h <- 0.7
            R_t[y] <- (4 * h * R0 * SB_t[y-1] / ( SB0*(1-h) + SB_t[y-1]*(5*h-1))) * exp(RecDev[y])
        }
        
        ## age-structured dynamics
        for(a in 1:length(L_a)){
            if(a==1){
                N_at[a,y] <- R_t[y]
            }
            if(a>1 & a<length(L_a)){
                N_at[a,y] <- N_at[a-1,y-1]*exp(-M-F_t[y-1]*S_a[a-1])
            }
            if(a==length(L_a)){
                N_at[a,y] <- (N_at[a-1,y-1] + N_at[a,y-1])*exp(-M-F_t[y-1]*S_a[a-1])
            }
        }
        ## spawning biomass
        SB_t[y] <- sum((N_at[,y] * W_a * Mat_a)[-1])
        ## catch
        Cn_at[,y] <- N_at[,y] * (1-exp(-M-F_t[y]*S_a)) * (F_t[y]*S_a)/(M+F_t[y]*S_a)
    }
    Cn_t <- colSums(Cn_at)
    N_t <- colSums(N_at[-1,])
    D_t <- SB_t/SB0

    if(sample!=FALSE){
        C_t <- Cn_t*sample
    }
    if(sample==FALSE){
        C_t <- Cn_t
    }
    CPUE_t <- qcoef * SB_t * exp(EffDev)

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
            names(C_tout) <- names(CPUE_tout) <- 1:Nyears

    LFout <- LF[-c(1:nburn),]
        rownames(LFout) <- 1:Nyears

    R_tout <- R_t[-c(1:nburn)]
    N_tout <- N_t[-c(1:nburn)]
    SB_tout <- SB_t[-c(1:nburn)]
    D_tout <- D_t[-c(1:nburn)]
    F_tout <- F_t[-c(1:nburn)]
    L_tout <- L_t[-c(1:nburn)]
    N_atout <- N_at[,-c(1:nburn)]

    if(alt_yrs==TRUE){
        yrs <- Nyears:1
        index <- rep(c(1,1,0,0), Nyears/4)
        yr_vec <- rev(yrs[which(index==1)])

        C_tout <- C_tout[which(names(C_tout) %in% yr_vec)]
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

    DataList <- list("I_t"=CPUE_tout, "C_t"=C_tout, "DataScenario"=modname,
        "LF"=LFout, "SigmaR"=SigmaR, "R_t"=R_tout, "N_t"=N_tout, "SB_t"=SB_tout, "D_t"=D_tout, "F_t"=F_tout, 
        "L_t"=L_tout, "N_at"=N_atout, "Amat"=Amat, "Mat_a"=Mat_a, "SB0"=SB0, "Nyears"=Nyears, "L_a"=L_a,
        "W_a"=W_a, "AgeMax"=AgeMax, "M"=M, "S_a"=S_a, "plb"=plb, "plba"=plba, "page"=page, "R0"=R0, 
        "SPR"=SPR, "SPR_t"=SPR_t,
        "ML_t"=ML_t, "nlbins"=length(mids))

    return(DataList)
}



#' Simulate spatial variation in growth
#'
#' \code{spatialgrowth_sim} simulate 1D Spatial variation in asymptotic length across sites
#'
#' @param n_i number of sites
#' @param Scale Gaussian scale parameter, default 2
#' @param Sigma2 variance in asymptotic length across sites, default 1
#' @param SD_spatial Gaussian variation (standard deviation), default 0.1
#' @param linf average asymptotic length across sites
#' @param beta_y trend in asymptotic length over sites, default 0.02

#' @return data frame of asymptotic length at each site
#' @export
spatialgrowth_sim <- function(n_i, Scale=2, Sigma2=1, SD_spatial=0.1, linf, beta_y=0.02){
    # require(RandomFields)

    ## sample locations
    lat_min <- -4
    lat_max <- -1
    y_i <- runif(n=n_i, min=lat_min, max=lat_max)

    ## simulate spatial process
    RMmodel <- RMgauss(var=SD_spatial^2, scale=Scale)
    linf_i1 <- linf * exp(RFsimulate(model=RMmodel, x=rep(0, n_i), y=y_i)@data[,1] - Sigma2/2) * exp( beta_y*(y_i - mean(c(lat_min, lat_max))))
    linf_i <- linf_i1 + (linf-mean(linf_i1))

    df <- data.frame(linf_i=linf_i, y_i=y_i)
    return(df)
}

#' Calculate AIC
#'
#' \code{calc_AIC} calculate AIC for model runs
#'
#' @param modpath_vec path to model run

#' @return matrix with AIC, AICc, deltaAIC, and deltaAICc (in that order by column) for each model (down the rows)
#' @export
calc_AIC <- function(modpath_vec){
    aic_mat <- matrix(NA, nrow=length(modpath_vec), ncol=6)
    colnames(aic_mat) <- c("AIC", "AICc", "deltaAIC", "deltaAICc", "relLikeAIC", "relLikeAICc")

    for(i in 1:length(modpath_vec)){
        input <- readRDS(file.path(modpath_vec[i], "Inputs2.rds"))
        report <- readRDS(file.path(modpath_vec[i], "Report.rds"))

        nll <- report$jnll
        params <- length(as.vector(unlist(input$Parameters))) - length(which(grepl("Nu_input", names(unlist(input$Parameters)))))
        sampsize <- length(which(input$Data$I_t>0)) + sum(input$Data$obs_per_yr) + length(which(input$Data$C_t>0))

        AIC <- 2*params + 2*nll
        AICc <- AIC + (2*params*(params + 1))/(sampsize - params - 1)

        aic_mat[i,1] <- AIC
        aic_mat[i,2] <- AICc
    }

    aic_mat[,"deltaAIC"] <- aic_mat[,"AIC"] - min(aic_mat[,"AIC"])
    aic_mat[,"deltaAICc"] <- aic_mat[,"AICc"] - min(aic_mat[,"AICc"])
    aic_mat[,"relLikeAIC"] <- sapply(1:nrow(aic_mat), function(x) exp((min(aic_mat[,"AIC"]) -aic_mat[x,"AIC"])/2))
    aic_mat[,"relLikeAICc"] <- sapply(1:nrow(aic_mat), function(x) exp((min(aic_mat[,"AICc"]) - aic_mat[x,"AICc"])/2))

    return(aic_mat)
}