#' TMB input formatting
#'
#' \code{format_input} Formats data, parameters, random effects, and mapped parameters for TMB input
#'
#' @param input tagged list of LIME inputs. Output from create_inputs.
#' @param data_avail types of data included, must at least include LCX where X is the number of years of length composition data. May also include "Catch" or "Index" separated by underscore. For example, "LC10", "Catch_LC1", "Index_Catch_LC20".
#' @param Fpen penalty on fishing mortality 0= off, 1=on
#' @param SigRpen penalty on sigmaR, 0=off, 1=on
#' @param SigRprior vector with prior info for sigmaR penalty, first term is the mean and second term is the standard deviation
#' @param est_sigma list of variance parameters to estimate, must match parameter names: log_sigma_R, log_sigma_C, log_sigma_I, log_CV_L, log_sigma_F
#' @param REML FALSE==off, TRUE==on
#' @param fix_f year (by index, not actual year) to start estimating fishing mortality (e.g. year 11 out of 20 to get estimates for last 10 years); the value of F in this year will be used as the estimate and SE for all previous years. 0=estimate all years.
#' @param f_startval default=NULL and F starting values are at 0 for all years. Can also specify vector of F starting values for all years to be modeled (can start at truth for debugging).
#' @param fix_param default=FALSE - parameters are fixed depending on the data available. Can also list vector of parameter names to fix at their starting values (use param_adjust and val_adjust to set these adjustments)
#' @param C_opt  default=0, NO catch data available. Copt=1 means the catch is in numbers, Copt2 means the catch is in weight. 

#' @return List, a tagged list of Data, Parameters, Random, Map
#' @export
format_input <- function(input, data_avail, Fpen, SigRpen, SigRprior, est_sigma, REML, fix_f, f_startval, fix_param=FALSE, C_opt=0){

    with(input, {
        ## data-rich model
        if(grepl("Index",data_avail) & grepl("Catch",data_avail) & grepl("LC",data_avail)){
            Data <- list(n_t=Nyears, n_lb=ncol(LF), 
                n_c=length(C_t),
                n_i=length(I_t), 
                n_lc=nrow(LF),
                n_ml=0, fix_f=as.vector(fix_f),
                T_yrs=1:Nyears, C_yrs=as.numeric(names(C_t)),
                I_yrs=as.numeric(names(I_t)),
                LC_yrs=as.numeric(rownames(LF)),
                ML_yrs=as.vector(0),
                obs_per_yr=obs_per_year,
                I_t=I_t, C_t=C_t, C_opt=C_opt,
                ML_t=as.vector(0), LF=LF,
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=highs, lbmids=mids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen, SigRpen=SigRpen, SigRprior=SigRprior)       
        }

        ## same as above but fit to mean length data instead of length composition data
        if(grepl("Index",data_avail) & grepl("Catch",data_avail) & grepl("ML",data_avail)){
            Data <- list(n_t=Nyears, n_lb=ncol(LF), 
                n_c=length(C_t),
                n_i=length(I_t), 
                n_lc=0,
                n_ml=nrow(LF), fix_f=as.vector(fix_f),
                T_yrs=1:Nyears, C_yrs=as.numeric(names(C_t)),
                I_yrs=as.numeric(names(I_t)),
                LC_yrs=as.vector(0),
                ML_yrs=as.numeric(rownames(LF)),
                obs_per_yr=obs_per_year,
                I_t=I_t, C_t=C_t, C_opt=C_opt,
                ML_t=rowMeans(LF), LF=as.matrix(0),
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=highs, lbmids=mids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen, SigRpen=SigRpen, SigRprior=SigRprior)       
        }

        ## index and length composition data
        if(grepl("Index",data_avail) & grepl("LC",data_avail) & grepl("Catch",data_avail)==FALSE){
            if(is.matrix(LF)){
                n_lc <- nrow(LF)
                LC_yrs <- as.numeric(rownames(LF))
                LF <- as.matrix(LF)
            }
            if(is.vector(LF)){
                n_lc <- 1
                LC_yrs <- Nyears
                LF <- t(as.matrix(LF))
            }
            Data <- list(n_t=Nyears, n_lb=ncol(LF), 
                n_c=0,
                n_i=length(I_t), 
                n_lc=n_lc,
                n_ml=0, fix_f=as.vector(fix_f),
                T_yrs=1:Nyears, C_yrs=as.vector(0),
                I_yrs=as.numeric(names(I_t)),
                LC_yrs=LC_yrs,
                ML_yrs=as.vector(0), 
                obs_per_yr=obs_per_year,
                I_t=I_t, C_t=as.vector(0), C_opt=C_opt,
                ML_t=as.vector(0), LF=LF,
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=highs, lbmids=mids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen, SigRpen=SigRpen, SigRprior=SigRprior)       
        }

        ## index and mean length data
        if(grepl("Index",data_avail) & grepl("ML",data_avail) & grepl("Catch",data_avail)==FALSE){
            if(is.matrix(LF)){
                n_ml <- nrow(LF)
                ML_yrs <- as.numeric(rownames(LF))
                LF <- as.matrix(LF)
            }
            if(is.vector(LF)){
                n_ml <- 1
                ML_yrs <- Nyears
                LF <- t(as.matrix(LF))
            }
            Data <- list(n_t=Nyears, n_lb=ncol(LF), 
                n_c=0,
                n_i=length(I_t), 
                n_lc=0,
                n_ml=n_ml, fix_f=as.vector(fix_f),
                T_yrs=1:Nyears, C_yrs=as.vector(0),
                I_yrs=as.numeric(names(I_t)),
                LC_yrs=as.vector(0),
                ML_yrs=ML_yrs,
                obs_per_yr=obs_per_year,
                I_t=I_t, C_t=as.vector(0), C_opt=C_opt,
                ML_t=rowMeans(LF), LF=as.matrix(0),
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=highs, lbmids=mids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen,  SigRpen=SigRpen, SigRprior=SigRprior)       
        }

        ## catch and length composition data
        if(grepl("Catch",data_avail) & grepl("LC",data_avail) & grepl("Index",data_avail)==FALSE){
            if(is.matrix(LF)){
                n_lc <- nrow(LF)
                LC_yrs <- as.numeric(rownames(LF))
                LF <- as.matrix(LF)
            }
            if(is.vector(LF)){
                n_lc <- 1
                LC_yrs <- Nyears
                LF <- t(as.matrix(LF))
            }
            Data <- list(n_t=Nyears, n_lb=ncol(LF), 
                n_c=length(C_t),
                n_i=0, 
                n_lc=n_lc,
                n_ml=0, fix_f=as.vector(fix_f),
                T_yrs=1:Nyears, C_yrs=as.numeric(names(C_t)),
                I_yrs=as.vector(0),
                LC_yrs=LC_yrs,
                ML_yrs=as.vector(0), 
                obs_per_yr=obs_per_year,
                I_t=as.vector(0), C_t=C_t, C_opt=C_opt,
                ML_t=as.vector(0), LF=LF,
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=highs, lbmids=mids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen, SigRpen=SigRpen, SigRprior=SigRprior)       
        }


        ## catch and mean length data
        if(grepl("Catch",data_avail) & grepl("ML",data_avail) & grepl("Index",data_avail)==FALSE){
            if(is.matrix(LF)){
                n_ml <- nrow(LF)
                ML_yrs <- as.numeric(rownames(LF))
                LF <- as.matrix(LF)
            }
            if(is.vector(LF)){
                n_ml <- 1
                ML_yrs <- Nyears
                LF <- t(as.matrix(LF))
            }
            Data <- list(n_t=Nyears, n_lb=ncol(LF), 
                n_c=length(C_t),
                n_i=0, 
                n_lc=0,
                n_ml=n_ml, fix_f=as.vector(fix_f),
                T_yrs=1:Nyears, C_yrs=as.numeric(names(C_t)),
                I_yrs=s.vector(0),
                LC_yrs=as.vector(0),
                ML_yrs=ML_yrs,
                obs_per_yr=obs_per_year,
                I_t=as.vector(0), C_t=C_t, C_opt=C_opt,
                ML_t=rowMeans(LF), LF=as.matrix(0),
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=highs, lbmids=mids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen,  SigRpen=SigRpen, SigRprior=SigRprior)       
        }

        ## length composition data only 
        if(grepl("LC",data_avail) & grepl("Index",data_avail)==FALSE & grepl("Catch",data_avail)==FALSE){
            if(is.matrix(LF)){
                n_lc <- nrow(LF)
                LC_yrs <- as.numeric(rownames(LF))
                LF <- as.matrix(LF)
            }
            if(is.vector(LF)){
                n_lc <- 1
                LC_yrs <- Nyears
                LF <- t(as.matrix(LF))
            }
            Data <- list(n_t=Nyears, n_lb=ncol(LF), 
                n_c=0,
                n_i=0, 
                n_lc=n_lc,
                n_ml=0, fix_f=as.vector(fix_f),
                T_yrs=1:Nyears, C_yrs=as.vector(0),
                I_yrs=as.vector(0),
                LC_yrs=LC_yrs,
                ML_yrs=as.vector(0),
                obs_per_yr=obs_per_year,
                I_t=as.vector(0), C_t=as.vector(0), C_opt=C_opt,
                ML_t=as.vector(0), LF=LF,
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=highs, lbmids=mids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen,  SigRpen=SigRpen, SigRprior=SigRprior)
        }       

        ## set input parameters - regardless of data availability 
        if(all(is.null(f_startval))) f_startval <- rep(1, Nyears)
        Parameters <- list(log_sigma_F=log(SigmaF), log_F_t_input=log(f_startval),log_q_I=log(qcoef), beta=log(R0), log_sigma_R=log(SigmaR), logS50=log(S50), log_sigma_C=log_sigma_C, log_sigma_I=log_sigma_I, log_CV_L=log_CV_L,Nu_input=rep(0,Nyears))

        ## turn off parameter estimation - depends on data availability
            Map = list()

            if("log_sigma_F" %in% est_sigma==FALSE){
                Map[["log_sigma_F"]] <- NA
                Map[["log_sigma_F"]] <- factor(Map[["log_sigma_F"]])
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

            if(all(fix_f!=0)){
                Map[["log_F_t_input"]] = 1:length(Parameters[["log_F_t_input"]])
                Map[["log_F_t_input"]][fix_f] <- NA
                Map[["log_F_t_input"]] <- factor(Map[["log_F_t_input"]])   
            }

                if(grepl("Catch",data_avail)==FALSE){        
                    Map[["beta"]] <- NA
                    Map[["beta"]] <- factor(Map[["beta"]])
                }

                if(grepl("Index",data_avail)==FALSE){
                    Map[["log_q_I"]] <- NA
                    Map[["log_q_I"]] <- factor(Map[["log_q_I"]])
                }

                if(grepl("LC",data_avail) & grepl("Index",data_avail)==FALSE & grepl("Catch",data_avail)==FALSE & grepl("Rich",data_avail)==FALSE & grepl("Moderate",data_avail)==FALSE & grepl("Sample",data_avail)==FALSE){                        
                   Map[["log_q_I"]] <- NA
                   Map[["log_q_I"]] <- factor(Map[["log_q_I"]])       
                   Map[["beta"]] <- NA
                   Map[["beta"]] <- factor(Map[["beta"]])
                }

            if(all(fix_param!=FALSE)){
                for(i in 1:length(fix_param)){
                    Map[[fix_param[i]]] <- NA
                    Map[[fix_param[i]]] <- factor(Map[[fix_param[i]]])
                }
            }

        if(length(Map)==0) Map <- NULL


        if(REML==FALSE) Random <- c("Nu_input")
        if(REML==TRUE){
            Random_vec <- c("Nu_input", "log_F_t_input", "log_q_I", "beta", "logS50") # 
            Random <- Random_vec[which(Random_vec %in% names(Map) == FALSE)]
        }
        if("log_sigma_F" %in% est_sigma) Random <- c(Random, "log_F_t_input")




    Return <- list("Parameters"=Parameters, "Data"=Data, "Random"=Random, "Map"=Map)
    return(Return)
    })

}