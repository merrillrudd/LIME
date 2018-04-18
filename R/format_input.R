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
#' @param LFdist likelihood distribution for length composition data, default=0 for multinomial, alternate=1 for dirichlet-multinomial
#' @param C_type  default=0, NO catch data available. Copt=1 means the catch is in numbers, Copt2 means the catch is in weight. 
#' @param est_more list of variance parameters to estimate, must match parameter names: log_sigma_R, log_sigma_C, log_sigma_I, log_CV_L, log_sigma_F
#' @param fix_more default=FALSE - parameters are fixed depending on the data available. Can also list vector of parameter names to fix at their starting values (use param_adjust and val_adjust to set these adjustments)
#' @param est_F_ft default=TRUE, otherwise 0 for off and 1 for on
#' @param f_startval_ft default=NULL and F starting values are at 0 for all years. Can also specify vector of F starting values for all years to be modeled (can start at truth for debugging).
#' @param rdev_startval_t default=NULL and Recruitment deviation starting values are at 0 for all years. Can also specify vector of recruitment deviation starting values for all years to be modeled (can start at truth for debugging)
#' @param est_selex_f default=TRUE to estimate selectivity parameters, can set to FALSE for all or multiple fleets
#' @param vals_selex_ft input selectivity-at-length (columns) by fleet (rows) - negative values in the first column indicate to estimate selectivity
#' @param Rdet default=FALSE to estimate recruitment deviations, TRUE=deterministic recruitment
#' @param mirror vector of parameter names to mirror between fleets
#' @param est_totalF TRUE estimate total F instead of by fleet
#' @param prop_f proportion of catch from each fleet
#' @return List, a tagged list of Data, Parameters, Random, Map
#' @export
format_input <- function(input, 
                        data_avail, 
                        Fpen, 
                        SigRpen, 
                        SigRprior, 
                        LFdist, 
                        C_type,
                        est_more, 
                        fix_more,
                        est_F_ft,
                        f_startval_ft, 
                        rdev_startval_t,
                        est_selex_f,
                        vals_selex_ft,
                        Rdet,
                        mirror,
                        est_totalF,
                        prop_f){

    with(input, {

        if(nseasons==1){
            S_yrs_inp <- 1:Nyears
            Nyears2 <- Nyears
        }
        if(nseasons>1){
            Nyears2 <- ceiling(Nyears/nseasons)
            S_yrs_inp <- unlist(lapply(1:Nyears2, function(x) rep(x, nseasons)))
        }

        selex_type_f <- rep(1,nfleets)
        for(i in 1:nfleets){
            if(any(vals_selex_ft[i,] > 0)) selex_type_f[i] <- 0
        }

            mirror_theta_inp <- ifelse("log_theta" %in% mirror, 1, 0)
            mirror_q_inp <- ifelse("log_q_f" %in% mirror, 1, 0)

            if(all(is.null(input$neff_ft))){
                if(is.vector(LF[,,1])==FALSE) n_inp <- t(sapply(1:nfleets, function(x) rowSums(LF[,,x])))
                if(is.vector(LF[,,1])) n_inp <- t(sapply(1:nfleets, function(x) sum(LF[,,x])))
            }
            if(all(is.null(input$neff_ft)==FALSE)) n_inp <- input$neff_ft

        if(Rdet==FALSE) Rdet <- 0
        if(Rdet==TRUE) Rdet <- 1

        # if(all(est_F_ft == TRUE)){
        #     indexF_ft <- matrix(1:Nyears2, nrow=nfleets, ncol=Nyears2)
        # }
        # if(all(est_F_ft == TRUE)==FALSE){
        #     indexF_ft <- matrix(1:Nyears2, nrow=nfleets, ncol=Nyears2)
        #     for(i in 1:nfleets){
        #         sub <- est_F_ft[i,]
        #         off <- which(sub == 0)
        #         new <- sapply(1:length(off), function(x){
        #             good <- which(sub == 1)
        #             if(off[x]>max(good)){
        #                 return(max(good))
        #             }
        #             if(off[x]<max(good)){
        #                 return(good[which(good > off[x])][1])
        #             }
        #         })
        #         indexF_ft[i,off] <- indexF_ft[i,new]
        #     }
        # }


        ## data-rich model
        if(grepl("Index",data_avail) & grepl("Catch",data_avail) & grepl("LC",data_avail)){

            # if(is.matrix(LF)){
            #     n_lc <- nrow(LF)
            #     LC_yrs <- as.numeric(rownames(LF))
            #     LF <- as.matrix(LF)
            #     if(is.null(input$neff_ft)) n_inp <- t(sapply(1:nfleets, function(x) rowSums(LF[,,x])))
            #     if(is.null(input$neff_ft)==FALSE) n_inp <- input$neff_ft
            # }
            # if(is.vector(LF)){
            #     n_lc <- 1
            #     LC_yrs <- Nyears
            #     LF <- t(as.matrix(LF))
            #     if(is.null(input$neff_ft)) n_inp <- sum(LF)
            #     if(is.null(input$neff_ft)==FALSE) n_inp <- input$neff_ft
            # }   
            Data <- list("n_t"=dim(LF)[1],
                         "n_lb"=dim(LF)[2],
                         "n_fl"=dim(LF)[3],
                         "n_a"=length(ages),
                         "LF_tlf"=LF,
                         "n_LF_ft"=n_inp,
                         "I_ft"=I_ft,
                         "C_ft"=C_ft,
                         "C_type"=C_type,
                         "ML_ft"=matrix(-999, nrow=dim(LF)[3], ncol=dim(LF)[1]),
                         "ages"=ages,
                         "match_ages"=seq(min(ages), max(ages), by=1),
                         "L_a"=L_a,
                         "W_a"=W_a,
                         "M"=M,
                         "h"=h,
                         "Mat_a"=Mat_a,
                         "lbhighs"=highs,
                         "lbmids"=mids,
                         "Fpen"=Fpen,
                         "SigRpen"=SigRpen,
                         "SigRprior"=SigRprior,
                         "selex_type_f"=selex_type_f,
                         "vals_selex_ft"=vals_selex_ft,
                         "Rdet"=Rdet,
                         "LFdist"=LFdist,
                         "S_yrs"=S_yrs_inp,
                         "n_s"=nseasons,
                         "n_y"=Nyears2,
                         "mirror_theta"=mirror_theta_inp,
                         "mirror_q"=mirror_q_inp,
                         # "indexF_ft"=indexF_ft,
                         "est_totalF"=ifelse(est_totalF==TRUE,1,0),
                         "prop_f"=prop_f)   
        }

        ## index and length composition data
        if(grepl("Index",data_avail) & grepl("LC",data_avail) & grepl("Catch",data_avail)==FALSE){
            # if(is.matrix(LF)){
            #     n_lc <- nrow(LF)
            #     LC_yrs <- as.numeric(rownames(LF))
            #     LF <- as.matrix(LF)
            #     if(is.null(input$obs_per_year)) n_inp <- rowSums(LF)
            #     if(is.null(input$obs_per_year)==FALSE) n_inp <- input$obs_per_year
            # }
            # if(is.vector(LF)){
            #     n_lc <- 1
            #     LC_yrs <- Nyears
            #     LF <- t(as.matrix(LF))
            #     if(is.null(input$obs_per_year)) n_inp <- sum(LF)
            #     if(is.null(input$obs_per_year)==FALSE) n_inp <- input$obs_per_year
            # }  

            Data <- list("n_t"=dim(LF)[1],
                         "n_lb"=dim(LF)[2],
                         "n_fl"=dim(LF)[3],
                         "n_a"=length(ages),
                         "LF_tlf"=LF,
                         "n_LF_ft"=n_inp,
                         "I_ft"=I_ft,
                         "C_ft"=matrix(-999, nrow=dim(LF)[3], ncol=dim(LF)[1]),
                         "C_type"=0,
                         "ML_ft"=matrix(-999, nrow=dim(LF)[3], ncol=dim(LF)[1]),
                         "ages"=ages,
                         "match_ages"=seq(min(ages), max(ages), by=1),
                         "L_a"=L_a,
                         "W_a"=W_a,
                         "M"=M,
                         "h"=h,
                         "Mat_a"=Mat_a,
                         "lbhighs"=highs,
                         "lbmids"=mids,
                         "Fpen"=Fpen,
                         "SigRpen"=SigRpen,
                         "SigRprior"=SigRprior,
                         "selex_type_f"=selex_type_f,
                         "vals_selex_ft"=vals_selex_ft,
                         "Rdet"=Rdet,
                         "LFdist"=LFdist,
                         "S_yrs"=S_yrs_inp,
                         "n_s"=nseasons,
                         "n_y"=Nyears2,
                         "mirror_theta"=mirror_theta_inp,
                         "mirror_q"=mirror_q_inp,
                         # "indexF_ft"=indexF_ft,
                         "est_totalF"=ifelse(est_totalF==TRUE,1,0),
                         "prop_f"=prop_f)   
        }


        ## catch and length composition data
        if(grepl("Catch",data_avail) & grepl("LC",data_avail) & grepl("Index",data_avail)==FALSE){
            # if(is.matrix(LF)){
            #     n_lc <- nrow(LF)
            #     LC_yrs <- as.numeric(rownames(LF))
            #     LF <- as.matrix(LF)
            #     if(is.null(input$obs_per_year)) n_inp <- rowSums(LF)
            #     if(is.null(input$obs_per_year)==FALSE) n_inp <- input$obs_per_year
            # }
            # if(is.vector(LF)){
            #     n_lc <- 1
            #     LC_yrs <- Nyears
            #     LF <- t(as.matrix(LF))
            #     if(is.null(input$obs_per_year)) n_inp <- sum(LF)
            #     if(is.null(input$obs_per_year)==FALSE) n_inp <- input$obs_per_year
            # }  

            Data <- list("n_t"=dim(LF)[1],
                         "n_lb"=dim(LF)[2],
                         "n_fl"=dim(LF)[3],
                         "n_a"=length(ages),
                         "LF_tlf"=LF,
                         "n_LF_ft"=n_inp,
                         "I_ft"=matrix(-999, nrow=dim(LF)[3], ncol=dim(LF)[1]),
                         "C_ft"=C_ft,
                         "C_type"=C_type,
                         "ML_ft"=matrix(-999, nrow=dim(LF)[3], ncol=dim(LF)[1]),
                         "ages"=ages,
                         "match_ages"=seq(min(ages), max(ages), by=1),
                         "L_a"=L_a,
                         "W_a"=W_a,
                         "M"=M,
                         "h"=h,
                         "Mat_a"=Mat_a,
                         "lbhighs"=highs,
                         "lbmids"=mids,
                         "Fpen"=Fpen,
                         "SigRpen"=SigRpen,
                         "SigRprior"=SigRprior,
                         "selex_type_f"=selex_type_f,
                         "vals_selex_ft"=vals_selex_ft,
                         "Rdet"=Rdet,
                         "LFdist"=LFdist,
                         "S_yrs"=S_yrs_inp,
                         "n_s"=nseasons,
                         "n_y"=Nyears2,
                         "mirror_theta"=mirror_theta_inp,
                         "mirror_q"=mirror_q_inp,
                         # "indexF_ft"=indexF_ft,
                         "est_totalF"=ifelse(est_totalF==TRUE,1,0),
                         "prop_f"=prop_f)      
        }

        ## length composition data only 
        if(grepl("LC",data_avail) & grepl("Index",data_avail)==FALSE & grepl("Catch",data_avail)==FALSE){
            # if(is.matrix(LF)){
            #     n_lc <- nrow(LF)
            #     LC_yrs <- as.numeric(rownames(LF))
            #     LF <- as.matrix(LF)
            #     if(is.null(input$obs_per_year)) n_inp <- rowSums(LF)
            #     if(is.null(input$obs_per_year)==FALSE) n_inp <- input$obs_per_year
            # }
            # if(is.vector(LF)){
            #     n_lc <- 1
            #     LC_yrs <- Nyears
            #     LF <- t(as.matrix(LF))
            #     if(is.null(input$obs_per_year)) n_inp <- sum(LF)
            #     if(is.null(input$obs_per_year)==FALSE) n_inp <- input$obs_per_year
            # }  

            Data <- list("n_t"=dim(LF)[1],
                         "n_lb"=dim(LF)[2],
                         "n_fl"=dim(LF)[3],
                         "n_a"=length(ages),
                         "LF_tlf"=LF,
                         "n_LF_ft"=n_inp,
                         "I_ft"=matrix(-999, nrow=dim(LF)[3], ncol=dim(LF)[1]),
                         "C_ft"=matrix(-999, nrow=dim(LF)[3], ncol=dim(LF)[1]),
                         "C_type"=0,
                         "ML_ft"=matrix(-999, nrow=dim(LF)[3], ncol=dim(LF)[1]),
                         "ages"=ages,
                         "match_ages"=seq(min(ages), max(ages), by=1),
                         "L_a"=L_a,
                         "W_a"=W_a,
                         "M"=M,
                         "h"=h,
                         "Mat_a"=Mat_a,
                         "lbhighs"=highs,
                         "lbmids"=mids,
                         "Fpen"=Fpen,
                         "SigRpen"=SigRpen,
                         "SigRprior"=SigRprior,
                         "selex_type_f"=selex_type_f,
                         "vals_selex_ft"=vals_selex_ft,
                         "Rdet"=Rdet,
                         "LFdist"=LFdist,
                         "S_yrs"=S_yrs_inp,
                         "n_s"=nseasons,
                         "n_y"=Nyears2,
                         "mirror_theta"=mirror_theta_inp,
                         "mirror_q"=mirror_q_inp,
                         # "indexF_ft"=indexF_ft,
                         "est_totalF"=ifelse(est_totalF==TRUE,1,0),
                         "prop_f"=prop_f)   
        }       

        ## set input parameters - regardless of data availability 
        if(all(is.null(f_startval_ft)==FALSE)){
            checkdim <- dim(f_startval_ft)
            if(est_totalF==TRUE & checkdim[1]!=1 & checkdim[2]!=Nyears2) stop("if estimating total F, provide start values for only one fleet")
        }
        if(all(is.null(f_startval_ft))){
            if(est_totalF==TRUE) f_startval_ft <- matrix(0.2, nrow=1, ncol=Nyears2)
            if(est_totalF==FALSE) f_startval_ft <- t(sapply(1:nfleets, function(x) rep(0.2, Nyears2)))
        }
        if(all(is.null(rdev_startval_t))) rdev_startval_t <- rep(0, Nyears2)
        Parameters <- list("log_F_ft"=log(f_startval_ft),
                        "log_q_f"=rep(log(qcoef), Data$n_fl),
                        "beta"=log(R0),
                        "log_sigma_R"=log(SigmaR),
                        "log_S50_f"=log(SL50),
                        "log_Sdelta_f"=log(SL95 - SL50), 
                        "log_sigma_F"=log(SigmaF), 
                        "log_sigma_C"=log(SigmaC),
                        "log_sigma_I"=log(SigmaI),
                        "log_CV_L"=log(CVlen),
                        "log_theta"=log(rep(theta, Data$n_fl)),
                        "Nu_input"=rdev_startval_t)

        ## turn off parameter estimation 
        Map = list()

            ## based on data availability
            if(grepl("Catch",data_avail)==FALSE){        
                Map[["beta"]] <- NA
                Map[["beta"]] <- factor(Map[["beta"]])
            }
            if(grepl("Index",data_avail)==FALSE){
                Map[["log_q_f"]] <- rep(NA, length(Parameters[["log_q_f"]]))
                Map[["log_q_f"]] <- factor(Map[["log_q_f"]])
            }

            ## based on function inputs
            if(all(fix_more!=FALSE)){
                for(i in 1:length(fix_more)){
                    Map[[fix_more[i]]] <- rep(NA, length(Parameters[[fix_more[i]]]))
                    Map[[fix_more[i]]] <- factor(Map[[fix_more[i]]])
                }
            }
            if("log_sigma_F" %in% est_more==FALSE){
                Map[["log_sigma_F"]] <- NA
                Map[["log_sigma_F"]] <- factor(Map[["log_sigma_F"]])
            }
            if("log_sigma_C" %in% est_more==FALSE){
                Map[["log_sigma_C"]] <- NA
                Map[["log_sigma_C"]] <- factor(Map[["log_sigma_C"]])
            }
            if("log_sigma_I" %in% est_more==FALSE){
                Map[["log_sigma_I"]] <- NA
                Map[["log_sigma_I"]] <- factor(Map[["log_sigma_I"]])
            }
            if("log_CV_L" %in% est_more==FALSE){
                Map[["log_CV_L"]] <- NA
                Map[["log_CV_L"]] <- factor(Map[["log_CV_L"]])
            }
            if(all(est_more==FALSE)){
                Map[["log_sigma_F"]] <- NA
                Map[["log_sigma_F"]] <- factor(Map[["log_sigma_F"]])
                Map[["log_CV_L"]] <- NA
                Map[["log_CV_L"]] <- factor(Map[["log_CV_L"]]) 
                Map[["log_sigma_C"]] <- NA
                Map[["log_sigma_C"]] <- factor(Map[["log_sigma_C"]])  
                Map[["log_sigma_I"]] <- NA
                Map[["log_sigma_I"]] <- factor(Map[["log_sigma_I"]])
            }

            if(LFdist==0){
                Map[["log_theta"]] <- rep(NA, length(Parameters[["log_theta"]]))
                Map[["log_theta"]] <- factor(Map[["log_theta"]])
            }

            if(any(est_selex_f == FALSE)){
                Map[["log_S50_f"]] <- Parameters$log_S50_f
                if(all(est_selex_f==FALSE) & nfleets > 1) Map[["log_S50_f"]][which(est_selex_f==FALSE)] <- NA
                if(all(est_selex_f==FALSE)) Map[["log_S50_f"]] <- rep(NA, length(Parameters[["log_S50_f"]]))
                Map[["log_S50_f"]] <- factor(Map[["log_S50_f"]])

                Map[["log_Sdelta_f"]] <- Parameters$log_Sdelta_f
                if(all(est_selex_f==FALSE) & nfleets > 1) Map[["log_Sdelta_f"]][which(est_selex_f==FALSE)] <- NA
                if(all(est_selex_f==FALSE)) Map[["log_Sdelta_f"]] <- rep(NA, length(Parameters[["log_Sdelta_f"]]))
                Map[["log_Sdelta_f"]] <- factor(Map[["log_Sdelta_f"]])
            }

            if(any(vals_selex_ft >= 0)){
                Map[["log_S50_f"]] <- Parameters$log_S50_f
                Map[["log_S50_f"]][which(vals_selex_ft[,1] >= 0)] <- NA
                Map[["log_S50_f"]] <- factor(Map[["log_S50_f"]])

                Map[["log_Sdelta_f"]] <- Parameters$log_Sdelta_f
                Map[["log_Sdelta_f"]][which(vals_selex_ft[,1] >= 0)] <- NA
                Map[["log_Sdelta_f"]] <- factor(Map[["log_Sdelta_f"]])                
            }

            if(Rdet==1){
                Map[["log_sigma_R"]] <- NA
                Map[["log_sigma_R"]] <- factor(Map[["log_sigma_R"]])

                Map[["Nu_input"]] <- rep(NA, length(Parameters$Nu_input))
                Map[["Nu_input"]] <- factor(Map[["Nu_input"]])
            }

            if("log_theta" %in% mirror){
                Map[["log_theta"]] <- c(Parameters$log_theta[1], rep(NA, (length(Parameters$log_theta)-1)))
                Map[["log_theta"]] <- factor(Map[["log_theta"]])
            }

            if(any(est_F_ft == 0)){
                fstart <- Parameters$log_F_ft + matrix(rnorm(length(Parameters$log_F_ft), 0, 2), nrow=nrow(Parameters$log_F_ft), ncol=ncol(Parameters$log_F_ft))
                for(i in 1:nfleets){
                    fstart[i,which(est_F_ft[i,]==0)] <- NA
                }
                Map[["log_F_ft"]] <- fstart
                Map[["log_F_ft"]] <- factor(Map[["log_F_ft"]])
            }

        if(length(Map)==0) Map <- NULL


        if(Rdet==0 | Rdet==FALSE) Random <- c("Nu_input")
        if(Rdet==1 | Rdet==TRUE) Random <- NULL
        # if(REML==TRUE){
        #     Random_vec <- c("Nu_input", "log_F_t_input", "log_q_I", "beta", "logS50") # 
        #     Random <- Random_vec[which(Random_vec %in% names(Map) == FALSE)]
        # }
        # if("log_sigma_F" %in% est_sigma) Random <- c(Random, "log_F_t_input")




    Return <- list("Parameters"=Parameters, "Data"=Data, "Random"=Random, "Map"=Map)
    return(Return)
    })

}