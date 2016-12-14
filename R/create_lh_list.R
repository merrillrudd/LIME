#' Create new life history list
#'
#' \code{create_lh_list} Creates list of life history information
#'
#' @param vbk von Bertalanffy k Brody growth coefficient
#' @param linf von Bertalanffy Linf asymptotic length
#' @param lwa length-weight scaling parameter
#' @param lwb length-weight allometric parameter
#' @param S50 starting value for age or length at 50 percent selectivity (will be estimated in LIME method)
#' @param M50 age or length at 50 percent maturity
#' @param selex_input specify whether argument S50 is an age or a length (default length)
#' @param maturity_input specify whether argument M50 is an age or a length (default length)
#' @param binwidth width of length bins (default = 1)
#' @param t0 theoretical age at length=0 (default = -0.01); avoid fixing to zero due to some issues with the first age/length bin
#' @param CVlen CV of the growth curve (default = 0.1)
#' @param SigmaC standard deviation - observation error of catch data (default = 0.2)
#' @param SigmaI standard deviation - observation error of index data (default = 0.2)
#' @param SigmaR standard deviation - process error for recruitment time series (default = 0.6 -- starting value, will be estimated)
#' @param SigmaF standard deviation - process error for fishing mortality time series (default = 0.3)
#' @param R0 equilibrium recruitment (default = 1); when no information on scale is available, will estimate relative deviations around equilibrium 1
#' @param h steepness parameter (default = 1)
#' @param qcoef starting value for catchability coefficient (when index data is available, default = 1e-5)
#' @param M value for natural mortality if there has been a study (default = NULL, calculated internally from vbk)
#' @param F1 starting value for initial fishing mortality. Default = 0.2, do not start at zero because this is used to set the initial values for estimating annual fishing mortality in log space, thus setting to zero would cause an error. 
#' @param Fequil equilibrium fishing mortality rate (used for simulation; default=0.2)
#' @param Frate parameter used to simulate fishing moratality time series (default=NULL)
#' @param Fmax maximum F used in simulation (default=NULL)
#' @param start_ages age to start (either 0 or 1; default = 0)
#' @param rho first-order autocorrelation in recruitment residuals parameter, default=0 (recruitment not autocorrelated)
#' @return List, a tagged list of life history traits
#' @export
create_lh_list <- function(vbk, linf, lwa, lwb, S50, M50, selex_input="length", maturity_input="length", binwidth=1, t0=-0.01, CVlen=0.1, SigmaC=0.2, SigmaI=0.2, SigmaR=0.6, SigmaF=0.3, R0=1,  h=1, qcoef=1e-5, M=NULL, F1=0.2, Fequil=0.2, Frate=0.2, Fmax=0.7, start_ages=0, rho=0){
            
    ## mortality
    if(is.null(M)) M <- 1.5*vbk  ## based on vbk if not specified 
    AgeMax <- ceiling(-log(0.01)/M)
    ages <- start_ages:AgeMax

    if(selex_input=="length"){
        SL50 <- S50
        S50 <- ceiling(t0-log(1-(SL50/linf))/vbk)

    }
    if(selex_input=="age"){
        SL50 <- ceiling(linf*(1-exp(-vbk*(S50-t0))))
    }

    if(maturity_input=="length"){
        ML50 <- M50
        M50 <- ceiling(t0-log(1-(ML50/linf))/vbk)
    }
    if(maturity_input=="age"){
        ML50 <- ceiling(linf*(1-exp(-vbk*(M50-t0))))
    }
    
    ## length bins
    mids <- seq((binwidth/2), linf*1.5, by=binwidth) 
    highs <- mids + (binwidth/2)
    lows <- mids - (binwidth)/2

    ## growth at age
    L_a <- linf*(1-exp(-vbk*(ages - t0)))
    W_a <- lwa*L_a^lwb  

    ## maturity
    if(maturity_input=="length"){
        Mat_l <- 1 / (1 + exp(ML50 - mids))
        Mat_ages <- ceiling(t0-log(1-(mids/linf))/vbk)
        names(Mat_l) <- Mat_ages
        Mat_a <- rep(NA, length(ages))
        for(a in 1:length(ages)){
            if(start_ages==0){
                if(a==1) Mat_a[a] <- 1e-20
                if(a>1){
                    fill <- Mat_l[which(names(Mat_l)==(a-1))][length(Mat_l[which(names(Mat_l)==(a-1))])]
                    if(length(fill)==1) Mat_a[a] <- fill
                    if(length(fill)==0) Mat_a[a] <- Mat_a[a-1]
                }
            }
            if(start_ages!=0){
                fill <- Mat_l[which(names(Mat_l)==a)][length(Mat_l[which(names(Mat_l)==a)])]
                if(length(fill)==1) Mat_a[a] <- fill
                if(length(fill)==0) Mat_a[a] <- Mat_a[a-1]
            }
        }       
    }
    if(maturity_input=="age"){
        Mat_a <- rep(NA, length(ages))
        if(start_ages==0) Mat_a <- c(1e-20, 1/(1+exp(M50 - ages[-1])))
        if(start_ages!=0) Mat_a <- 1/(1+exp(M50 - ages))
    }

    ## selectivity 
    if(start_ages==0) S_a <- c(1e-20, 1 / (1 + exp(S50 - ages[-1]))) # Selectivity at age
    if(start_ages!=0) S_a <- 1/(1+exp(S50-ages))
        
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
    Outs$SigmaF <- SigmaF
    Outs$R0 <- R0
    Outs$lwa <- lwa
    Outs$lwb <- lwb
    Outs$S50 <- S50
    Outs$SL50 <- SL50
    Outs$h <- h
    Outs$qcoef <- qcoef
    Outs$M <- M
    Outs$F1 <- F1
    Outs$AgeMax <- AgeMax
    Outs$mids <- mids
    Outs$highs <- highs
    Outs$lows <- lows
    Outs$S_a <- S_a
    Outs$L_a <- L_a
    Outs$W_a <- W_a
    Outs$M50 <- M50
    Outs$ML50 <- ML50
    Outs$Mat_a <- Mat_a
    Outs$Fequil <- Fequil
    Outs$Frate <- Frate
    Outs$Fmax <- Fmax
    Outs$rho <- rho
    return(Outs)
}
