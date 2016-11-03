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
                        # next
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
