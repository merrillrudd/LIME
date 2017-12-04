#' Create new life history list
#'
#' \code{create_lh_list} Creates list of life history information
#'
#' @author M.B. Rudd
#' @param vbk von Bertalanffy k Brody growth coefficient
#' @param linf von Bertalanffy Linf asymptotic length
#' @param lwa length-weight scaling parameter
#' @param lwb length-weight allometric parameter
#' @param M value for natural mortality if there has been a study (default = NULL, calculated internally from vbk)
#' @param AgeMax option to specify maximum age; default=NULL will calculate as the age at which 1 percent of individuals are left in the unfished condition
#' @param F1 starting value for initial fishing mortality. Default = 0.2, do not start at zero because this is used to set the initial values for estimating annual fishing mortality in log space, thus setting to zero would cause an error. 
#' @param S50 starting value for age or length at 50 percent selectivity (will be estimated in LIME method)
#' @param M50 age or length at 50 percent maturity
#' @param S95 default=NULL for one-parameter logistic model; starting value for age or length at 95 percent selectivity
#' @param M95 default=NULL for one-parameter logistic model; age or length at 50 percent maturity
#' @param Sslope default=NULL, option to specify slope of logistic curve for length-at-selectivity
#' @param Mslope default=NULL option to specify slope of logistic curve for length-at-maturity
#' @param selex_input specify whether argument S50 is an age or a length (default length)
#' @param maturity_input specify whether argument M50 is an age or a length (default length)
#' @param selex_type default="logistic" for 1-parameter logistic selex, alternate="dome" for dome-shaped selectivity and must specify dome-params LV and RV.
#' @param dome_sd standard deviation of normal distribution to the right side of the fully selected age/length 
#' @param binwidth width of length bins (default = 1)
#' @param t0 theoretical age at length=0 (default = -0.01); avoid fixing to zero due to some issues with the first age/length bin
#' @param R0 equilibrium recruitment (default = 1); when no information on scale is available, will estimate relative deviations around equilibrium 1
#' @param h steepness parameter (default = 1)
#' @param CVlen CV of the growth curve (default = 0.1)
#' @param SigmaC standard deviation - observation error of catch data (default = 0.2)
#' @param SigmaI standard deviation - observation error of index data (default = 0.2)
#' @param SigmaR standard deviation - process error for recruitment time series (default = 0.6 -- starting value, will be estimated)
#' @param SigmaF standard deviation - process error for fishing mortality time series (default = 0.3)
#' @param qcoef starting value for catchability coefficient (when index data is available, default = 1e-5)
#' @param Fequil equilibrium fishing mortality rate (used for simulation; default=0.2)
#' @param Frate parameter used to simulate fishing moratality time series (default=NULL)
#' @param Fmax maximum F used in simulation (default=NULL)
#' @param start_ages age to start (either 0 or 1; default = 0)
#' @param rho first-order autocorrelation in recruitment residuals parameter, default=0 (recruitment not autocorrelated)
#' @param theta dirichlet-multinomial parameter related to effective sample size. default to 10, will not be used if length frequency distribution LFdist is set to multinomial (0). Only used if distribution is dirichlet-multinomial (LFdist=1)
#' @param nseasons specify number of sub-time periods per year; default=1 (instantaneous sampling)
#' @importFrom stats pnorm
#' 
#' 
#' @return List, a tagged list of life history traits
#' @export
create_lh_list <- function(vbk, linf, lwa, lwb, S50, M50, S95=NULL, M95=NULL, Sslope=NULL, Mslope=NULL, selex_input="length", maturity_input="length", selex_type="logistic", dome_sd=NULL, binwidth=1, t0=-0.01, CVlen=0.1, SigmaC=0.001, SigmaI=0.001, SigmaR=0.737, SigmaF=0.2, R0=1,  h=1, qcoef=1e-5, M=NULL, AgeMax=NULL, F1=0.2, Fequil=0.5, Frate=0.2, Fmax=0.7, start_ages=0, rho=0, theta=10, nseasons=1){
            
    ## mortality
    if(is.null(M)) M <- 1.5*vbk  ## based on vbk if not specified 
    if(is.null(AgeMax)) AgeMax <- ceiling(-log(0.01)/M)
    ages <- seq(start_ages, to=(AgeMax+1 - (1/nseasons)), by=(1/nseasons))

    ## monthly mortality
    M <- M/nseasons

    ## length bins
    mids <- seq((binwidth/2), linf*1.5, by=binwidth) 
    highs <- mids + (binwidth/2)
    lows <- mids - (binwidth)/2


    ## growth at age
    L_a <- linf*(1-exp(-vbk*(ages - t0)))
    W_a <- lwa*L_a^lwb  
    W_l <- lwa*mids^lwb

    ##probability of being in a length bin given age
    lbprobs <- function(mnl,sdl) return(pnorm(highs,mnl,sdl)-pnorm(lows,mnl,sdl))
    vlprobs <- Vectorize(lbprobs,vectorize.args=c("mnl","sdl"))
    plba_a <- t(vlprobs(L_a,L_a*CVlen))
    plba_a <- plba_a/rowSums(plba_a)

    ## maturity and selectivity
    if(is.null(M95)) mat_param <- 1
    if(is.null(M95)==FALSE) mat_param <- 2
    if(is.null(Mslope)==FALSE) mat_param <- 3
    if(is.null(S95)) sel_param <- 1
    if(is.null(S95)==FALSE) sel_param <- 2
    if(is.null(Sslope)==FALSE) sel_param <- 3


    if(selex_input=="length"){
        SL50 <- S50
        S50 <- ceiling(t0-log(1-(SL50/linf))/vbk)
        if(is.null(S95)==FALSE){
            SL95 <- S95
            S95 <- ceiling(t0-log(1-(SL95/linf))/vbk)
        }
    }
    if(selex_input=="age"){
        SL50 <- ceiling(linf*(1-exp(-vbk*(S50-t0))))
        if(is.null(S95)==FALSE) SL95 <- ceiling(linf*(1-exp(-vbk*(S95-t0))))
    }

    if(maturity_input=="length"){
        ML50 <- M50
        M50 <- ceiling(t0-log(1-(ML50/linf))/vbk)
        if(is.null(M95)==FALSE){
            ML95 <- M95
            M95 <- ceiling(t0-log(1-(ML95/linf))/vbk)
        }
    }
    if(maturity_input=="age"){
        ML50 <- ceiling(linf*(1-exp(-vbk*(M50-t0))))
        if(is.null(M95)==FALSE) ML95 <- ceiling(linf*(1-exp(-vbk*(M95-t0))))
    }
    
    ## maturity
    if(mat_param==1){
        ## specified length bins
        Mat_l <- (1 / (1 + exp(ML50 - mids)))
        ## 1 cm length bins default
        # Mat_l2 <- 1 / (1 + exp(ML50 - mids2))
    }
    if(mat_param==2){
        Mat_l <- (1 /(1 + exp(-log(19)*(mids-ML50)/(ML95-ML50)))) # Maturity at length
        # Mat_l2 <- 1 /(1 + exp(-log(19)*(mids2-ML50)/(ML95-ML50)))
    }
    if(mat_param==3){
        Mat_l <- (1 /(1 + exp(-((mids-ML50)/Mslope))))
    }

    if(maturity_input=="length"){
        Mat_a <- apply(t(plba_a)*Mat_l, 2, sum)
    }
    if(maturity_input=="age"){
        Mat_a <- rep(NA, length(ages))
        if(mat_param==1){
           Mat_a <- 1/(1+exp(M50 - ages))
        }
        if(mat_param==2){
           Mat_a <- 1/(1+exp(-log(19)*(ages-M50)/(M95-M50)))
        }
    }
    if(is.null(M95)){
        id_L95 <- which(round(Mat_a, 2) %in% seq(from=0.92,to=1.00,by=0.01))[1]
        ML95 <- L_a[id_L95]
        M95 <- ceiling(t0-log(1-(ML95/linf))/vbk)
    }


    ## selectivity
    if(sel_param==1){
        S_l <- (1 / (1 + exp(SL50 - mids)))
    }
    if(sel_param==2){
        S_l <- (1 /(1 + exp(-log(19)*(mids-SL50)/(SL95-SL50)))) # Selectivity-at-Length
    }
    if(sel_param==3){
        S_l <- (1 /(1 + exp(-((mids-SL50)/Sslope))))
    }

    if(selex_input=="length"){
        if(selex_type=="dome"){
            Sfull <- which(round(S_l,2)==1.00)[1]
            if(is.na(Sfull)) Sfull <- which(round(S_l,1)==1.00)[1]
            index <- (Sfull+1):length(S_l)
            S_l[index] <- exp((-(index-Sfull)^2)/(2*dome_sd^2))
        }
        S_a <- apply(t(plba_a)*S_l, 2, sum)
    }
    if(selex_input=="age"){
        S_a <- rep(NA, length(ages))
        if(sel_param==1){
           S_a <- 1/(1+exp(S50 - ages))
        }
        if(sel_param==2){
            S_a <- 1/(1+exp(-log(19)*(ages-S50)/(S95-S50)))
        }
        if(selex_type=="dome"){
            Sfull <- which(round(S_a,1)==1.00)[1]
            index <- (Sfull+1):length(S_a)
            S_a[index] <- exp((-(index-Sfull)^2)/(2*dome_sd^2))
        }
    }
    if(is.null(S95)){
        id_L95 <- which(round(S_a, 2) %in% seq(from=0.92,to=1.00,by=0.01))[1]
        SL95 <- L_a[id_L95]
        S95 <- ceiling(t0-log(1-(SL95/linf))/vbk)
    }
    if(selex_type!="dome") Sfull <- NULL

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
    Outs$S95 <- S95
    Outs$SL50 <- SL50
    Outs$SL95 <- SL95
    Outs$Sfull <- Sfull
    Outs$dome_sd <- dome_sd
    Outs$selex_type <- selex_type
    Outs$selex_input <- selex_input
    Outs$maturity_input <- maturity_input
    Outs$h <- h
    Outs$qcoef <- qcoef
    Outs$M <- M
    Outs$F1 <- F1
    Outs$AgeMax <- AgeMax
    Outs$ages <- ages
    Outs$mids <- mids
    Outs$highs <- highs
    Outs$lows <- lows
    Outs$S_a <- S_a
    Outs$S_l <- S_l
    Outs$L_a <- L_a
    Outs$W_a <- W_a
    Outs$W_l <- W_l
    Outs$M50 <- M50
    Outs$M95 <- M95
    Outs$ML50 <- ML50
    Outs$ML95 <- ML95
    Outs$Mat_a <- Mat_a
    Outs$Mat_l <- Mat_l
    Outs$Fequil <- Fequil
    Outs$Frate <- Frate
    Outs$Fmax <- Fmax
    Outs$rho <- rho
    Outs$theta <- theta
    Outs$nseasons <- nseasons
    return(Outs)
}
