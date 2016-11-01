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

        if(length(Mat_a)==AgeMax) Mat_a <- c(1e-20, Mat_a)
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