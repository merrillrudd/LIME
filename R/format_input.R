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
#' @param f_startval_ft default=NULL and F starting values are at 0 for all years. Can also specify vector of F starting values for all years to be modeled (can start at truth for debugging).
#' @param rdev_startval_t default=NULL and Recruitment deviation starting values are at 0 for all years. Can also specify vector of recruitment deviation starting values for all years to be modeled (can start at truth for debugging)
#' @param est_selex_f default=TRUE to estimate selectivity parameters, can set to FALSE for all or multiple fleets
#' @param randomR default = TRUE, estimate recruitment as a random effect; if FALSE, turn off random effect on recruitment (do not derive deviations)
#' 
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
                        f_startval_ft, 
                        rdev_startval_t,
                        est_selex_f,
                        randomR){

    with(input, {

        if(nseasons==1){
            S_yrs_inp <- 1:Nyears
            Nyears2 <- Nyears
        }
        if(nseasons>1){
            Nyears2 <- ceiling(Nyears/nseasons)
            S_yrs_inp <- years_i
        }
        selex_type_f <- sapply(1:nfleets, function(x) ifelse(selex_type[x]=="logistic",1, ifelse(selex_type[x]=="dome",2,0)))
        if(any(selex_type_f==0)) stop("specify selex_type in create_lh_list")


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
            if(is.null(input$neff_ft)) n_inp <- t(sapply(1:nfleets, function(x) rowSums(LF[,,x])))
            if(is.null(input$neff_ft)==FALSE) n_inp <- input$neff_ft

            Data <- list("n_t"=dim(LF)[1],
                         "n_lb"=dim(LF)[2],
                         "n_fl"=dim(LF)[3],
                         "n_a"=length(ages),
                         "LF_tlf"=LF,
                         "n_lc_ft"=n_inp,
                         "I_ft"=I_ft,
                         "C_ft"=C_ft,
                         "C_type"=C_type,
                         "ML_ft"=as.matrix(0),
                         "ages"=ages,
                         "L_a"=L_a,
                         "W_a"=W_a,
                         "M"=M,
                         "h"=h,
                         "Mat_a"=Mat_a,
                         "lbhighs"=seq(binwidth, by=binwidth, length=dim(LF)[2]),
                         "lbmids"=seq(binwidth/2, by=binwidth, length=dim(LF)[2]),
                         "Fpen"=Fpen,
                         "SigRpen"=SigRpen,
                         "SigRprior"=SigRprior,
                         "selex_type_f"=selex_type_f,
                         "LFdist"=LFdist,
                         "S_yrs"=S_yrs_inp,
                         "n_s"=nseasons,
                         "n_y"=length(Nyears2))   
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
            if(is.null(input$neff_ft)) n_inp <- t(sapply(1:nfleets, function(x) rowSums(LF[,,x])))
            if(is.null(input$neff_ft)==FALSE) n_inp <- input$neff_ft

            Data <- list("n_t"=dim(LF)[1],
                         "n_lb"=dim(LF)[2],
                         "n_fl"=dim(LF)[3],
                         "n_a"=length(ages),
                         "LF_tlf"=LF,
                         "n_lc_ft"=n_inp,
                         "I_ft"=I_ft,
                         "C_ft"=as.matrix(0),
                         "C_type"=0,
                         "ML_ft"=as.matrix(0),
                         "ages"=ages,
                         "L_a"=L_a,
                         "W_a"=W_a,
                         "M"=M,
                         "h"=h,
                         "Mat_a"=Mat_a,
                         "lbhighs"=seq(binwidth, by=binwidth, length=dim(LF)[2]),
                         "lbmids"=seq(binwidth/2, by=binwidth, length=dim(LF)[2]),
                         "Fpen"=Fpen,
                         "SigRpen"=SigRpen,
                         "SigRprior"=SigRprior,
                         "selex_type_f"=selex_type_f,
                         "LFdist"=LFdist,
                         "S_yrs"=S_yrs_inp,
                         "n_s"=nseasons,
                         "n_y"=length(Nyears2))   
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
            if(is.null(input$neff_ft)) n_inp <- t(sapply(1:nfleets, function(x) rowSums(LF[,,x])))
            if(is.null(input$neff_ft)==FALSE) n_inp <- input$neff_ft

            Data <- list("n_t"=dim(LF)[1],
                         "n_lb"=dim(LF)[2],
                         "n_fl"=dim(LF)[3],
                         "n_a"=length(ages),
                         "LF_tlf"=LF,
                         "n_lc_ft"=n_inp,
                         "I_ft"=as.matrix(0),
                         "C_ft"=C_ft,
                         "C_type"=C_type,
                         "ML_ft"=as.matrix(0),
                         "ages"=ages,
                         "L_a"=L_a,
                         "W_a"=W_a,
                         "M"=M,
                         "h"=h,
                         "Mat_a"=Mat_a,
                         "lbhighs"=seq(binwidth, by=binwidth, length=dim(LF)[2]),
                         "lbmids"=seq(binwidth/2, by=binwidth, length=dim(LF)[2]),
                         "Fpen"=Fpen,
                         "SigRpen"=SigRpen,
                         "SigRprior"=SigRprior,
                         "selex_type_f"=selex_type_f,
                         "LFdist"=LFdist,
                         "S_yrs"=S_yrs_inp,
                         "n_s"=nseasons,
                         "n_y"=length(Nyears2))      
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
            if(is.null(input$neff_ft)) n_inp <- t(sapply(1:nfleets, function(x) rowSums(LF[,,x])))
            if(is.null(input$neff_ft)==FALSE) n_inp <- input$neff_ft

            Data <- list("n_t"=dim(LF)[1],
                         "n_lb"=dim(LF)[2],
                         "n_fl"=dim(LF)[3],
                         "n_a"=length(ages),
                         "LF_tlf"=LF,
                         "n_lc_ft"=n_inp,
                         "I_ft"=as.matrix(0),
                         "C_ft"=as.matrix(0),
                         "C_type"=0,
                         "ML_ft"=as.matrix(0),
                         "ages"=ages,
                         "L_a"=L_a,
                         "W_a"=W_a,
                         "M"=M,
                         "h"=h,
                         "Mat_a"=Mat_a,
                         "lbhighs"=seq(binwidth, by=binwidth, length=dim(LF)[2]),
                         "lbmids"=seq(binwidth/2, by=binwidth, length=dim(LF)[2]),
                         "Fpen"=Fpen,
                         "SigRpen"=SigRpen,
                         "SigRprior"=SigRprior,
                         "selex_type_f"=selex_type_f,
                         "LFdist"=LFdist,
                         "S_yrs"=S_yrs_inp,
                         "n_s"=nseasons,
                         "n_y"=length(Nyears2))   
        }       

        ## set input parameters - regardless of data availability 
        if(all(is.null(f_startval_ft))) f_startval_ft <- t(sapply(1:nfleets, function(x) rep(0.2, Nyears2)))
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
                Map[["theta"]] <- NA
                Map[["theta"]] <- factor(Map[["theta"]])
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

            if(randomR==FALSE){
                Map[["log_sigma_R"]] <- NA
                Map[["log_sigma_R"]] <- factor(Map[["log_sigma_R"]])

                Map[["Nu_input"]] <- rep(NA, length(Parameters$Nu_input))
                Map[["Nu_input"]] <- factor(Map[["Nu_input"]])
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