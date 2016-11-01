#' model fits or kobe plots
#'
#' \code{plotResults} plot results for time series of estimated/derived parameters or kobe plots
#'
#' @param Data read RDS from model directory
#' @param Report read RDS from model directory
#' @param Sdreport read RDS from model directory
#' @param Derived_quants read RDS from model directory
#' @param flag_convergence TRUE- flag in directory that model didn't converge, FALSE = no flag in directory
#' @param parameter which parameter to plot results for
#' @param xaxt plot xaxis, TRUE or FALSE
#' @param ylab plot ylab, TRUE or FALSE
#' @param simulation get different items if these are application or simulation results

#' @return displays plot
#' 
#' @details possible values for parameter argument for model fits: "B" biomass, "N" abundance, "ML" mean length, "R" recruitment, "F" fishing mortality, "D" relative biomass, "C" catch, "I" abundance index, for kobe plot, "kobe" will show SPR compared with F/F30 and F/F40
#' @export
plotResults <- function(Data, Report, Sdreport, Derived_quants, flag_convergence, parameter, xaxt=TRUE, ylab=FALSE, simulation=TRUE){
        
        FUN = function(InputMat, log=TRUE){
          if(log==TRUE) return(c( exp(InputMat[,1]-1.96*InputMat[,2]), rev(exp(InputMat[,1]+1.96*InputMat[,2]))))
          if(log==FALSE) return(c( InputMat[,1]-1.96*InputMat[,2], rev(InputMat[,1]+1.96*InputMat[,2])))
        } 

      if(simulation==TRUE) years <- 1:length(Data$SB_t)
      if(simulation==FALSE){
        years <- Data$years_i
        years_real <- Data$years
      }
          if(parameter=="B"){
            ## Biomass
            if(simulation==TRUE) Mat <- cbind("Year"=years, "Data"=Data$SB_t, "Est"=Report$SB_t_hat)
            if(simulation==FALSE) Mat <- cbind("Year"=years, "Data"=rep(NA, length(years)), "Est"=Report$SB_t_hat)
            ymax <- 3
            matplot(y=Mat[,c("Data", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE) if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lSB_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            mtext("Biomass", side=3, line=-2, font=2)
            if(ylab==TRUE) mtext("Biomass", side=2, line=3, font=2)
          }
          if(parameter=="N"){
            ## Abundance
            if(simulation==TRUE) Mat <- cbind("Year"=years, "Data"=Data$N_t, "Est"=Report$N_t_hat)
            if(simulation==FALSE) Mat <- cbind("Year"=years, "Data"=rep(NA, length(years)), "Est"=Report$N_t_hat)
            ymax <- 4
            matplot(y=Mat[,c("Data", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lN_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            mtext("Abundance", side=3, line=-2, font=2)
            if(ylab==TRUE) mtext("Abundance", side=2, line=3, font=2)

          }
          if(parameter=="ML"){
            ## Average Length
            tml <- Data$ML_t
            pml <- rep(NA, length(years))
              names(pml) <- years
            pml[which(names(pml) %in% names(Data$ML_t))] <- tml
            Mat <- cbind("Year"=years, "Data"=pml, "Est"=Report$L_t_hat)
            ymax <- max(Report$L_t_hat)*1
            matplot(y=Mat[,c("Data", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="L_t_hat"),], log=FALSE), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            mtext("Mean Length", side=3, line=-2, font=2)
            if(ylab==TRUE) mtext("Mean Length", side=2, line=3, font=2)
          }
          if(parameter=="R"){       
            ## Recruitment
            if(simulation==TRUE) Mat <- cbind("Year"=years, "Data"=Data$R_t, "Est"=Report$R_t_hat)
            if(simulation==FALSE) Mat <- cbind("Year"=years, "Data"=rep(NA, length(years)), "Est"=Report$R_t_hat)
            ymax <- 3
            matplot(y=Mat[,c("Data", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE) if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            mtext("Recruitment", side=3, line=-2, font=2)
            if(ylab==TRUE) mtext("Recruitment", side=2, line=3, font=2)

          }       
          if(parameter=="F"){
            ## Fishing Mortality
            NAs <- rep(NA, length(years))
            if(simulation==TRUE) Mat <- cbind("Year"=years, "Data"=Data$F_t, "Est"=Report$F_t_hat)
            if(simulation==FALSE) Mat <- cbind("Year"=years, "Data"=NAs, "Est"=Report$F_t)
            ymax <- 2
            matplot(y=Mat[,c("Data", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            mtext("Fishing Mortality", side=3, line=-2, font=2)
            if(ylab==TRUE) mtext("Fishing Mortality", side=2, line=3, font=2)

          }
          if(parameter=="D"){     
            ## Relative abundance
            if(simulation==TRUE) Mat <- cbind("Year"=years, "Data"=Data$D_t, "Est"=Report$Depl)
            if(simulation==FALSE) Mat <- cbind("Year"=years, "Data"=rep(NA, length(years)), "Est"=Report$Depl)
            ymax <- 1
            matplot(y=Mat[,c("Data", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lD_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            mtext("Relative abundance", side=3, line=-2, font=2)  
            if(ylab==TRUE) mtext("Relative abundance", side=2, line=3, font=2)
          }
          if(parameter=="C"){      
            ## Catch
            tcatch <- Data$C_t
            pcatch <- rep(NA, length(years))
              names(pcatch) <- years
            pcatch[which(names(pcatch) %in% names(Data$C_t))] <- tcatch
            Mat <- cbind("Year"=years, "Data"=pcatch, "Est"=Report$C_t_hat)
            ymax <- 5
            matplot(y=Mat[,c("Data", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lC_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            mtext("Catch", side=3, line=-2, font=2)   
            if(ylab==TRUE) mtext("Catch", side=2, line=3, font=2)
          }
          if(parameter=="I"){     
            ## Index
            tindex <- Data$I_t
            pindex <- rep(NA, length(years))
              names(pindex) <- years
            pindex[which(names(pindex) %in% names(Data$I_t))] <- tindex
            Mat <- cbind("Year"=years, "Data"=pindex, "Est"=Report$I_t_hat)
            ymax <- 3
            matplot(y=Mat[,c("Data", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lI_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            mtext("Index", side=3, line=-2, font=2)
            if(ylab==TRUE) mtext("Index", side=2, line=3, font=2)
        }
        if(xaxt==TRUE){
          axis(1, at=seq(1,20, by=5), labels=years_real[seq(1,20,by=5)])
          mtext(side=1, "Year", line=4)  
        }

        if(parameter=="kobe"){
          plot(x=1, y=1, type="n", xaxs="i", yaxs="i", xlim=c(0, 1), ylim=c(0, 3), xaxt="n", yaxt="n")
          abline(h=1, lty=2, lwd=2)
          abline(v=0.3, lty=2, col="forestgreen", lwd=2)
          abline(v=0.4, lty=2, col="steelblue", lwd=2)
          points(x=Derived_quants$SPR, y=Derived_quants$FF30, col="forestgreen", pch=19, cex=2)
          points(x=Derived_quants$SPR, y=Derived_quants$FF40, col="steelblue", pch=19, cex=2)
          axis(1, at=seq(0, 1, by=0.2))
          axis(2, at=seq(0, 3, by=0.2))
          mtext(side=1, "SPR", line=2)
          mtext(side=2, "F/Fref", line=2)
          legend("topright", col=c("forestgreen", "steelblue"), legend=c("30%", "40%"), title="Target SPR", pch=19)
        }
        if(flag==TRUE) mtext(side=3, "model not converged", col="red", font=2, outer=TRUE, line=2, cex=2)
  

}
