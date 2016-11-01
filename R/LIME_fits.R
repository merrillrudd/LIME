#' Plot LIME model fits
#'
#' \code{LIME_fits} create figure of LIME model fits for several population parameters
#'
#' @param Inputs list of inputs for model run
#' @param Report report file from model run
#' @param Sdreport sdreport file from model run
#' @param obsData observed data
#' @param save TRUE if you want to save the figure to directory, FALSE if you want to view figure

#' @return figure of LIME model fits for several population parameters
#' @export
LIME_fits <- function(Inputs, Report, Sdreport, obsData, save){

        FUN = function(InputMat, log=TRUE){
          if(log==TRUE) return(c( exp(InputMat[,1]-1.96*InputMat[,2]), rev(exp(InputMat[,1]+1.96*InputMat[,2]))))
          if(log==FALSE) return(c( InputMat[,1]-1.96*InputMat[,2], rev(InputMat[,1]+1.96*InputMat[,2])))
        } 

      Nyears <- Inputs$Data$n_t
  if(any(is.na(Sdreport))) next
  
  par(mfrow=c(3,2), mar=c(1,5,0,0), omi=c(0.7,0.5,0.5,0.5))
    Mat <- cbind("Year"=1:Nyears, "Est"=Report$F_t)
    ymax <- 5
    matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax), lwd=c(2,5), xaxt="n", yaxt="n", cex.axis=1.5, xaxs="i", yaxs="i", ylab="")
    axis(2, las=2, at=pretty(c(0,ymax-1)), cex.axis=1.5)
    if(all(is.na(Sdreport))==FALSE){
        if( !("condition" %in% names(attributes(Sdreport)))){
            polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        }
    }
    mtext(side=2, "Fishing mortality", line=3, cex=1.5)
  
    Mat <- cbind("Year"=1:Nyears, "Est"=Report$R_t_hat)
    ymax <- 3
    matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax), lwd=c(2,5), xaxt="n", yaxt="n", cex.axis=1.5, xaxs="i", yaxs="i", ylab="")
    axis(2, las=2, at=pretty(c(0,ymax-1)), cex.axis=1.5)
    if(all(is.na(Sdreport))==FALSE){
        if( !("condition" %in% names(attributes(Sdreport)))){
            polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        }
    }
    mtext(side=2, "Relative recruitment", line=3, cex=1.5)

    Mat <- cbind("Year"=1:Nyears, "Est"=Report$Depl)
    ymax <- 3
    matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax), lwd=c(2,5), xaxt="n", yaxt="n", cex.axis=1.5, xaxs="i", yaxs="i", ylab="")
    axis(2, las=2, at=pretty(c(0,ymax-1)), cex.axis=1.5)
    if(all(is.na(Sdreport))==FALSE){
        if( !("condition" %in% names(attributes(Sdreport)))){
            polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lD_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)    
        }
    }
    mtext(side=2, "Relative biomass", line=3, cex=1.5)

  Mat <- cbind("Year"=1:Nyears, "Est"=Report$SPR_t)
  ymax <- 1
    matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax), lwd=c(2,5), xaxt="n", yaxt="n", cex.axis=1.5, xaxs="i", yaxs="i", ylab="")
    axis(2, las=2, at=pretty(c(0,ymax)), cex.axis=1.5)
    if(all(is.na(Sdreport))==FALSE){
        if( !("condition" %in% names(attributes(Sdreport)))){
            polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),], log=FALSE), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)  
        }
    }
    mtext(side=2, "SPR", line=3, cex=1.5)

  Mat <- cbind("Year"=1:Nyears, "Est"=Report$I_t)
  ymax <- 0.1
    matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax), lwd=c(2,5), yaxt="n", cex.axis=1.5, xaxs="i", yaxs="i", ylab="")
    axis(2, las=2, at=pretty(c(0,ymax)), cex.axis=1.5)
    if(all(is.na(Sdreport))==FALSE){
        if( !("condition" %in% names(attributes(Sdreport)))){
            polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lI_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)  
        }
    }
    points(x=Inputs$Data$I_yrs, y=Inputs$Data$I_t, pch=19, cex=1.5)
    mtext(side=2, "CPUE index", line=3.5, cex=1.5)

  Mat <- cbind("Year"=1:Nyears, "Est"=Report$L_t_hat)
  ymax <- 50
    matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax), lwd=c(2,5), yaxt="n", cex.axis=1.5, xaxs="i", yaxs="i", ylab="")
    axis(2, las=2, at=pretty(c(0,ymax-1)), cex.axis=1.5)
    if(all(is.na(Sdreport))==FALSE){
        if( !("condition" %in% names(attributes(Sdreport)))){
            polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="L_t_hat"),], log=FALSE), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA) 
        }
    } 
    points(x=Inputs$Data$LC_yrs, y=obsData$ML_t, pch=19, cex=1.5)
    mtext(side=2, "Mean length", line=3, cex=1.5)

    mtext("Year", outer=TRUE, side=1, line=2, cex=1.5)

}