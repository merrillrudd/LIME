#' Plot relative error
#'
#' \code{plotRE} boxplots comparing relative error across several models
#'
#' @param modpath_vec directory to find model results
#' @param itervec number of iterations
#' @param modnames model names for x axis
#' @param xaxt default "n" to remove x axis, follows same rules as plot function
#' @param yaxt default "n" to remove y axis, follows same rules as plot function
#' @param value relative error for specified reference point, "SPR" or "FFref" (associated with F/F40)
#' @param ylim set y axis limits
#' @param col.plot color for boxplots
#' @param col.line color for line designating 0 relative error
#' @param yr year to plot relative error


#' @return vector of directory paths
#' @export
plotRE <- function(modpath_vec, itervec, modnames=NULL, xaxt="n", yaxt="n", value, ylim=c(-1,1.5), col.plot="steelblue", col.line="goldenrod", yr){

    RE <- t(sapply(itervec, function(x) calcRE(modpath_vec=modpath_vec, iter=x, value=value, yr=yr)))

    colnames(RE) <- modnames

    boxplot(RE, ylim=ylim, xaxt=xaxt, yaxt=yaxt, col=col.plot)
    abline(h=0, col=col.line, lwd=2)
    
    return(RE)
}
