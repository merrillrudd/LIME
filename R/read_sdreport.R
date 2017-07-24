#' read TMB Sdreport 
#'
#' \code{read_sdreport} read summary(Sdreport) to calculate confidence intervals for a parameter 
#' 
#' @param InputMat summary(Sdreport) which rownames have specific parameter name
#' @param log is the parameter in logspace?

#' @return vector with lower and upper 95 confidence intervals. Upper intervals are reversed for ease plotting using polygon function. 
#' 
#' @export
read_sdreport <- function(InputMat, log=TRUE){
          index <- which(is.na(InputMat[,2])==FALSE)
          if(log==TRUE) return(c( exp(InputMat[index,1]-1.96*InputMat[index,2]), rev(exp(InputMat[index,1]+1.96*InputMat[index,2]))))
          if(log==FALSE) return(c( InputMat[index,1]-1.96*InputMat[index,2], rev(InputMat[index,1]+1.96*InputMat[index,2])))
} 