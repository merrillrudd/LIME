#' Saved life history values
#'
#' \code{choose_lh_list} Fixed lists of life history/starting value information for pre-loaded species
#'
#' @param species species code name
#' @param selex "asymptotic"= assume asymptotic selectivity, "dome"= assume dome-shaped selectivity
#' 
#' @return List, a tagged list of life history/starting value information
#' @export
choose_lh_list <- function(species, selex, param_adjust=FALSE, val=FALSE, start_ages=0){

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
        Fequil <- 0.25 #0.34 rate from previous studies
        Frate <- 0.2
        Fmax <- 0.7

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
        F1 <- 0.05  ## change to 0.34 in real assessment

        ## sensitivities
        if("linf" %in% param_adjust) linf <- val[which(param_adjust=="linf")]
        if("vbk" %in% param_adjust) vbk <- val[which(param_adjust=="vbk")]
        if("M" %in% param_adjust) M <- val[which(param_adjust=="M")]
        if("CVlen" %in% param_adjust) CVlen <- val[which(param_adjust=="CVlen")]
        if("SigmaR" %in% param_adjust) SigmaR <- val[which(param_adjust=="SigmaR")]
        if("SigmaF" %in% param_adjust) SigmaF <- val[which(param_adjust=="SigmaF")]
        if("ML50" %in% param_adjust) ML50 <- val[which(param_adjust=="ML50")]
        if("R0" %in% param_adjust) R0 <- val[which(param_adjust=="R0")]
        if("binwidth" %in% param_adjust) binwidth <- val[which(param_adjust=="binwidth")]


        ## derived
        mids <- seq((binwidth/2), linf*1.5, by=binwidth) # from 120 cm
        highs <- mids + (binwidth/2)
        lows <- mids - (binwidth)/2
        
        ages <- start_ages:AgeMax
        M50 <- round(t0-log(1-(ML50/linf))/vbk)
        A95 <- M50+1
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

        ## selectivity 
        if(selex=="asymptotic"){
            if(start_ages==0) S_a <- c(1e-20, 1 / (1 + exp(S50 - ages[-1]))) # Selectivity at age
            if(start_ages!=0) S_a <- 1/(1+exp(S50-ages))
            Syoung <- NA
            Sold <- NA
        }
        if(selex=="dome"){
            S_a_calc <- rep(NA, length(ages))
            Syoung <- S50
            Sold <- ceiling(quantile(c(AgeMax,Syoung), probs=0.75))
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
        F1 <- 0.05
        Fmax <- 0.7

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

        ## derived
        M <- 1.49

        ## sensitivities
        if("linf" %in% param_adjust) linf <- val[which(param_adjust=="linf")]
        if("vbk" %in% param_adjust) vbk <- val[which(param_adjust=="vbk")]
        if("M" %in% param_adjust) M <- val[which(param_adjust=="M")]
        if("CVlen" %in% param_adjust) CVlen <- val[which(param_adjust=="CVlen")]
        if("SigmaR" %in% param_adjust) SigmaR <- val[which(param_adjust=="SigmaR")]
        if("SigmaF" %in% param_adjust) SigmaF <- val[which(param_adjust=="SigmaF")]
        if("ML50" %in% param_adjust) ML50 <- val[which(param_adjust=="ML50")]
        if("R0" %in% param_adjust) R0 <- val[which(param_adjust=="R0")]
        if("binwidth" %in% param_adjust) binwidth <- val[which(param_adjust=="binwidth")]

        AgeMax <- round(-log(0.01)/M)
        ages <- start_ages:AgeMax
        M50 <- round(t0-log(1-(ML50/linf))/vbk)
        A95 <- M50+1
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

        ## selectivity 
        if(selex=="asymptotic"){
            S_l <- 1/(1+exp(-log(19)*(mids-SL50)/(SL95-SL50)))
            S_ages <- round(t0-log(1-(mids/linf))/vbk)
            names(S_l) <- S_ages
            S_a <- rep(NA, length(ages))
            for(a in 1:length(ages)){
                if(start_ages==0){
                    if(a==1) S_a[a] <- 1e-20
                    if(a>1){
                        fill <- S_l[which(names(S_l)==(a-1))][length(S_l[which(names(S_l)==(a-1))])]
                        if(length(fill)==1) S_a[a] <- fill
                        if(length(fill)==0) S_a[a] <- S_a[a-1]
                    }
                }
                if(start_ages!=0){
                    fill <- S_l[which(names(S_l)==a)][length(S_l[which(names(S_l)==a)])]
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
            Sold <- ceiling(quantile(c(AgeMax,Syoung), probs=0.75))
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
        vbk <- 0.14
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
        F1 <- 0.05
        Fmax <- 0.7

        ## variation terms
        SigmaF <- 0.3
        SigmaC <- 0.2
        SigmaI <- 0.2
        SigmaR <- 0.6

        ## bins
        binwidth <- 1

        ## selectivity
        ML50 <- 70
        S50 <- 2
        SL50 <- round(t0-log(1-(2/linf))/vbk)

        ## derived
        M <- 0.15

        ## sensitivities
        if("linf" %in% param_adjust) linf <- val[which(param_adjust=="linf")]
        if("vbk" %in% param_adjust) vbk <- val[which(param_adjust=="vbk")]
        if("M" %in% param_adjust) M <- val[which(param_adjust=="M")]
        if("CVlen" %in% param_adjust) CVlen <- val[which(param_adjust=="CVlen")]
        if("SigmaR" %in% param_adjust) SigmaR <- val[which(param_adjust=="SigmaR")]
        if("SigmaF" %in% param_adjust) SigmaF <- val[which(param_adjust=="SigmaF")]
        if("ML50" %in% param_adjust) ML50 <- val[which(param_adjust=="ML50")]
        if("R0" %in% param_adjust) R0 <- val[which(param_adjust=="R0")]
        if("binwidth" %in% param_adjust) binwidth <- val[which(param_adjust=="binwidth")]

        AgeMax <- round(-log(0.01)/M)
        ages <- start_ages:AgeMax
        M50 <- round(t0-log(1-(ML50/linf))/vbk)
        A95 <- M50+1
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
        if(start_ages==0) Mat_a <- c(1e-20, 1 / (1 + exp(M50 - ages[-1])))
        if(start_ages!=0) Mat_a <- 1/(1+exp(M50 - ages))

        ## selectivity 
        if(selex=="asymptotic"){
            if(start_ages==0) S_a <- c(1e-20, 1 / (1 + exp(S50 - ages[-1]))) # Selectivity at age
            if(start_ages!=0) S_a <- 1/(1+exp(S50-ages))
            Syoung <- NA
            Sold <- NA
        }
        if(selex=="dome"){
            S_a_calc <- rep(NA, length(ages))
            Syoung <- S50
            Sold <- ceiling(quantile(c(AgeMax,Syoung), probs=0.75))
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
    Outs$SigmaF <- SigmaF
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
    Outs$ML95 <- ML95
    Outs$Mat_a <- Mat_a
    Outs$Fequil <- Fequil
    Outs$Frate <- Frate
    Outs$Fmax <- Fmax

    return(Outs)

}
