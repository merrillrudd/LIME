#' TMB input formatting
#'
#' \code{format_input} Formats data, parameters, random effects, and mapped parameters for TMB input
#'
#' @author M.B. Rudd
#' @param input tagged list of LIME inputs. Output from create_inputs.
#' @param data_avail types of data included, must at least include LCX where X is the number of years of length composition data. May also include "Catch" or "Index" separated by underscore. For example, "LC10", "Catch_LC1", "Index_Catch_LC20".
#' @param Fpen penalty on fishing mortality 0= off, 1=on
#' @param SigRpen penalty on sigmaR, 0=off, 1=on
#' @param SigRprior vector with prior info for sigmaR penalty, first term is the mean and second term is the standard deviation
#' @param est_sigma list of variance parameters to estimate, must match parameter names: log_sigma_R, log_sigma_C, log_sigma_I, log_CV_L, log_sigma_F
#' @param f_startval default=NULL and F starting values are at 0 for all years. Can also specify vector of F starting values for all years to be modeled (can start at truth for debugging).
#' @param fix_param default=FALSE - parameters are fixed depending on the data available. Can also list vector of parameter names to fix at their starting values (use param_adjust and val_adjust to set these adjustments)
#' @param fix_param_t default=FALSE - fix certain parameters in time series (e.g. fishing mortality, recruitment deviations) list with first item the name of the parameter and second item the numbers in the time series to be fixed. 
#' @param C_opt  default=0, NO catch data available. Copt=1 means the catch is in numbers, Copt2 means the catch is in weight. 
#' @param LFdist likelihood distribution for length composition data, default=0 for multinomial, alternate=1 for dirichlet-multinomial
#' @param S_l_input input a vector specifying selectivity-at-length, or set less than 0 to use 1-parameter logistic function for selectivity
#' @param theta_type if 0, estimate annual theta; if 1, estimate single theta for all years of length comp
#' @param randomR default = TRUE, estimate recruitment as a random effect; if FALSE, turn off random effect on recruitment (do not derive deviations)
#' 
#' @return List, a tagged list of Data, Parameters, Random, Map
#' @export
format_input <- function(input, data_avail, Fpen, SigRpen, SigRprior, est_sigma, f_startval, fix_param=FALSE, fix_param_t=FALSE, C_opt=0, LFdist, S_l_input, theta_type, randomR){

    with(input, {

        if(nseasons==1){
            S_yrs_inp <- 1:Nyears
            Nyears2 <- Nyears
        }
        if(nseasons>1){
            Nyears2 <- ceiling(Nyears/nseasons)
            S_yrs_inp <- years_i
        }

        ## data-rich model
        if(grepl("Index",data_avail) & grepl("Catch",data_avail) & grepl("LC",data_avail)){
            if(is.matrix(LF)){
                n_lc <- nrow(LF)
                LC_yrs <- as.numeric(rownames(LF))
                LF <- as.matrix(LF)
                if(is.null(input$obs_per_year)) n_inp <- rowSums(LF)
                if(is.null(input$obs_per_year)==FALSE) n_inp <- input$obs_per_year
            }
            if(is.vector(LF)){
                n_lc <- 1
                LC_yrs <- Nyears
                LF <- t(as.matrix(LF))
                if(is.null(input$obs_per_year)) n_inp <- sum(LF)
                if(is.null(input$obs_per_year)==FALSE) n_inp <- input$obs_per_year
            }   
            Data <- list(n_t=Nyears, n_lb=ncol(LF), 
                n_c=length(C_t),
                n_i=length(I_t), 
                n_lc=nrow(LF),
                n_ml=0,
                T_yrs=1:Nyears, C_yrs=as.numeric(names(C_t)),
                I_yrs=as.numeric(names(I_t)),
                LC_yrs=as.numeric(rownames(LF)),
                ML_yrs=as.vector(0),
                obs_per_yr=n_inp,
                I_t=I_t, C_t=C_t, C_opt=C_opt,
                ML_t=as.vector(0), LF=LF, LFdist=LFdist, theta_type=theta_type, lbhighs=highs, lbmids=mids,
                n_a=length(ages), ages=ages, L_a=L_a, W_a=W_a, M=M, h=h, Mat_a=Mat_a,
                Fpen=Fpen, SigRpen=SigRpen, SigRprior=SigRprior, S_l_input=S_l_input, S_yrs=S_yrs_inp, n_s=nseasons, n_y=max(S_yrs_inp))       
        }

        ## same as above but fit to mean length data instead of length composition data
        if(grepl("Index",data_avail) & grepl("Catch",data_avail) & grepl("ML",data_avail)){
            Data <- list(n_t=Nyears, n_lb=ncol(LF), 
                n_c=length(C_t),
                n_i=length(I_t), 
                n_lc=0,
                n_ml=nrow(LF),
                T_yrs=1:Nyears, C_yrs=as.numeric(names(C_t)),
                I_yrs=as.numeric(names(I_t)),
                LC_yrs=as.vector(0),
                ML_yrs=as.numeric(rownames(LF)),
                obs_per_yr=n_inp,
                I_t=I_t, C_t=C_t, C_opt=C_opt,
                ML_t=rowMeans(LF), LF=as.matrix(0), LFdist=LFdist, theta_type=theta_type, lbhighs=highs, lbmids=mids,
                n_a=length(ages), ages=ages, L_a=L_a, W_a=W_a, M=M, h=h, Mat_a=Mat_a,
                Fpen=Fpen, SigRpen=SigRpen, SigRprior=SigRprior, S_l_input=S_l_input, S_yrs=S_yrs_inp, n_s=nseasons, n_y=max(S_yrs_inp))       
        }

        ## index and length composition data
        if(grepl("Index",data_avail) & grepl("LC",data_avail) & grepl("Catch",data_avail)==FALSE){
            if(is.matrix(LF)){
                n_lc <- nrow(LF)
                LC_yrs <- as.numeric(rownames(LF))
                LF <- as.matrix(LF)
                if(is.null(input$obs_per_year)) n_inp <- rowSums(LF)
                if(is.null(input$obs_per_year)==FALSE) n_inp <- input$obs_per_year
            }
            if(is.vector(LF)){
                n_lc <- 1
                LC_yrs <- Nyears
                LF <- t(as.matrix(LF))
                if(is.null(input$obs_per_year)) n_inp <- sum(LF)
                if(is.null(input$obs_per_year)==FALSE) n_inp <- input$obs_per_year
            }  
            Data <- list(n_t=Nyears, n_lb=ncol(LF), 
                n_c=0,
                n_i=length(I_t), 
                n_lc=n_lc,
                n_ml=0,
                T_yrs=1:Nyears, C_yrs=as.vector(0),
                I_yrs=as.numeric(names(I_t)),
                LC_yrs=LC_yrs,
                ML_yrs=as.vector(0), 
                obs_per_yr=n_inp,
                I_t=I_t, C_t=as.vector(0), C_opt=C_opt,
                ML_t=as.vector(0), LF=LF, LFdist=LFdist, theta_type=theta_type, lbhighs=highs, lbmids=mids,
                n_a=length(ages), ages=ages, L_a=L_a, W_a=W_a, M=M, h=h, Mat_a=Mat_a,
                Fpen=Fpen, SigRpen=SigRpen, SigRprior=SigRprior, S_l_input=S_l_input, S_yrs=S_yrs_inp, n_s=nseasons, n_y=max(S_yrs_inp))       
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
                n_ml=n_ml,
                T_yrs=1:Nyears, C_yrs=as.vector(0),
                I_yrs=as.numeric(names(I_t)),
                LC_yrs=as.vector(0),
                ML_yrs=ML_yrs,
                obs_per_yr=n_inp,
                I_t=I_t, C_t=as.vector(0), C_opt=C_opt,
                ML_t=rowMeans(LF), LF=as.matrix(0), LFdist=LFdist, theta_type=theta_type, lbhighs=highs, lbmids=mids,
                n_a=length(ages), ages=ages, L_a=L_a, W_a=W_a, M=M, h=h, Mat_a=Mat_a,
                Fpen=Fpen,  SigRpen=SigRpen, SigRprior=SigRprior, S_l_input=S_l_input, S_yrs=S_yrs_inp, n_s=nseasons, n_y=max(S_yrs_inp))       
        }

        ## catch and length composition data
        if(grepl("Catch",data_avail) & grepl("LC",data_avail) & grepl("Index",data_avail)==FALSE){
            if(is.matrix(LF)){
                n_lc <- nrow(LF)
                LC_yrs <- as.numeric(rownames(LF))
                LF <- as.matrix(LF)
                if(is.null(input$obs_per_year)) n_inp <- rowSums(LF)
                if(is.null(input$obs_per_year)==FALSE) n_inp <- input$obs_per_year
            }
            if(is.vector(LF)){
                n_lc <- 1
                LC_yrs <- Nyears
                LF <- t(as.matrix(LF))
                if(is.null(input$obs_per_year)) n_inp <- sum(LF)
                if(is.null(input$obs_per_year)==FALSE) n_inp <- input$obs_per_year
            }  
            Data <- list(n_t=Nyears, n_lb=ncol(LF), 
                n_c=length(C_t),
                n_i=0, 
                n_lc=n_lc,
                n_ml=0,
                T_yrs=1:Nyears, C_yrs=as.numeric(names(C_t)),
                I_yrs=as.vector(0),
                LC_yrs=LC_yrs,
                ML_yrs=as.vector(0), 
                obs_per_yr=n_inp,
                I_t=as.vector(0), C_t=C_t, C_opt=C_opt,
                ML_t=as.vector(0), LF=LF, LFdist=LFdist, theta_type=theta_type, lbhighs=highs, lbmids=mids,
                n_a=length(ages), ages=ages, L_a=L_a, W_a=W_a, M=M, h=h, Mat_a=Mat_a,
                Fpen=Fpen, SigRpen=SigRpen, SigRprior=SigRprior, S_l_input=S_l_input, S_yrs=S_yrs_inp, n_s=nseasons, n_y=max(S_yrs_inp))       
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
                n_ml=n_ml,
                T_yrs=1:Nyears, C_yrs=as.numeric(names(C_t)),
                I_yrs=s.vector(0),
                LC_yrs=as.vector(0),
                ML_yrs=ML_yrs,
                obs_per_yr=n_inp,
                I_t=as.vector(0), C_t=C_t, C_opt=C_opt,
                ML_t=rowMeans(LF), LF=as.matrix(0), LFdist=LFdist, theta_type=theta_type, lbhighs=highs, lbmids=mids,
                n_a=length(ages), ages=ages, L_a=L_a, W_a=W_a, M=M, h=h, Mat_a=Mat_a,
                Fpen=Fpen,  SigRpen=SigRpen, SigRprior=SigRprior, S_l_input=S_l_input, S_yrs=S_yrs_inp, n_s=nseasons, n_y=max(S_yrs_inp))       
        }

        ## length composition data only 
        if(grepl("LC",data_avail) & grepl("Index",data_avail)==FALSE & grepl("Catch",data_avail)==FALSE){
            if(is.matrix(LF)){
                n_lc <- nrow(LF)
                LC_yrs <- as.numeric(rownames(LF))
                LF <- as.matrix(LF)
                if(is.null(input$obs_per_year)) n_inp <- rowSums(LF)
                if(is.null(input$obs_per_year)==FALSE) n_inp <- input$obs_per_year
            }
            if(is.vector(LF)){
                n_lc <- 1
                LC_yrs <- Nyears
                LF <- t(as.matrix(LF))
                if(is.null(input$obs_per_year)) n_inp <- sum(LF)
                if(is.null(input$obs_per_year)==FALSE) n_inp <- input$obs_per_year
            }  
            Data <- list(n_t=Nyears, n_lb=ncol(LF), 
                n_c=0,
                n_i=0, 
                n_lc=n_lc,
                n_ml=0,
                T_yrs=1:Nyears, C_yrs=as.vector(0),
                I_yrs=as.vector(0),
                LC_yrs=LC_yrs,
                ML_yrs=as.vector(0),
                obs_per_yr=n_inp,
                I_t=as.vector(0), C_t=as.vector(0), C_opt=C_opt,
                ML_t=as.vector(0), LF=LF, LFdist=LFdist, theta_type=theta_type, lbhighs=highs, lbmids=mids,
                n_a=length(ages), ages=ages, L_a=L_a, W_a=W_a, M=M, h=h, Mat_a=Mat_a,
                Fpen=Fpen,  SigRpen=SigRpen, SigRprior=SigRprior, S_l_input=S_l_input, S_yrs=S_yrs_inp, n_s=nseasons, n_y=max(S_yrs_inp))
        }       

        ## set input parameters - regardless of data availability 
        if(all(is.null(f_startval))) f_startval <- rep(1, Nyears2)
        if(theta_type==0) input_theta <- rep(log(theta), Data$n_lc)
        if(theta_type==1) input_theta <- log(theta)
        Parameters <- list(log_sigma_F=log(SigmaF), log_F_t_input=log(f_startval),log_q_I=log(qcoef), beta=log(R0), log_sigma_R=log(SigmaR), logS50=log(SL50), logSdelta=log(SL95-SL50), log_sigma_C=log_sigma_C, log_sigma_I=log_sigma_I, log_CV_L=log_CV_L, log_theta=input_theta, Nu_input=rep(0,Nyears2))

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

            if(LFdist==0){
                Map[["theta"]] <- NA
                Map[["theta"]] <- factor(Map[["theta"]])
            }

            if(S_l_input[1]>=0){
                Map[["logS50"]] <- NA
                Map[["logS50"]] <- factor(Map[["logS50"]])

                Map[["logSdelta"]] <- NA
                Map[["logSdelta"]] <- factor(Map[["logSdelta"]])
            }

            if(randomR==FALSE){
                Map[["log_sigma_R"]] <- NA
                Map[["log_sigma_R"]] <- factor(Map[["log_sigma_R"]])

                Map[["Nu_input"]] <- rep(NA, length(Parameters$Nu_input))
                Map[["Nu_input"]] <- factor(Map[["Nu_input"]])
            }

            if(all(fix_param_t[[1]]!=FALSE)){
                for(i in 1:length(fix_param_t)){
                    index <- 1:length(Parameters[[fix_param_t[[i]][1]]])
                    index[as.numeric(fix_param_t[[i]][2:length(fix_param_t[[1]])])] <- NA
                    Map[[fix_param_t[[i]][1]]] <- index
                    Map[[fix_param_t[[i]][1]]] <- factor(Map[[fix_param_t[[i]][1]]])
                }
            }

        if(length(Map)==0) Map <- NULL


        if(randomR==TRUE) Random <- c("Nu_input")
        if(randomR==FALSE) Random <- NULL
        # if(REML==TRUE){
        #     Random_vec <- c("Nu_input", "log_F_t_input", "log_q_I", "beta", "logS50") # 
        #     Random <- Random_vec[which(Random_vec %in% names(Map) == FALSE)]
        # }
        # if("log_sigma_F" %in% est_sigma) Random <- c(Random, "log_F_t_input")




    Return <- list("Parameters"=Parameters, "Data"=Data, "Random"=Random, "Map"=Map)
    return(Return)
    })

}