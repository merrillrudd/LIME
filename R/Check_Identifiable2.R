#' Check for identifiability of fixed effects -- frm TMBhelpers (slightly adjusted by MBR)
#'
#' \code{Check_Identifiable} calculates the matrix of second-derivatives of the marginal likelihood w.r.t. fixed effects, to see if any linear combinations are unidentifiable
#'
#' @param obj, The compiled TMB object
#'
#' @return A tagged list of the hessian and the message
#' @details Slight adjustment made in the RowMax calculation - when the eigenvector was a vector, the apply function was not working 

#' @export
Check_Identifiable2 = function( obj ){
  # Finite-different hessian
  ParHat = TMBhelper:::extract_fixed( obj )
  List = NULL
  List[["Hess"]] = optimHess( par=ParHat, fn=obj$fn, gr=obj$gr )

  # Check eigendecomposition
  List[["Eigen"]] = eigen( List[["Hess"]] )
  List[["WhichBad"]] = which( List[["Eigen"]]$values < sqrt(.Machine$double.eps) )

  # Check for parameters
  RowMax = tryCatch(apply( List[["Eigen"]]$vectors[,List[["WhichBad"]]], MARGIN=1, FUN=function(vec){max(abs(vec))} ), error=function(e) NA)
  if(all(is.na(RowMax))) RowMax <- sapply(1:length(List[["Eigen"]]$vectors[,List[["WhichBad"]]]), function(x) max(abs(vec)))
  List[["BadParams"]] = data.frame("Param"=names(obj$par), "MLE"=ParHat, ifelse(RowMax>0.1, "Bad","OK"))

  # Message
  if( length(List[["WhichBad"]])==0 ){
    message( "All parameters are identifiable" )
  }else{
    print( List[["BadParams"]] )
  }

  # Return
  return( invisible(List) )
}
