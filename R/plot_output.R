#' plot LIME output
#'
#' \code{plot_output} plot output from LIME or LB-SPR
#'
#' @author M.B. Rudd
#' @param all_years vector of all years modeled
#' @param lc_years vector of years with length composition data
#' @param Inputs LIME input file
#' @param Report LIME report file
#' @param Sdreport LIME standard error file
#' @param LBSPR LBSPR results - either straight from LBSPR output or LIME-edited
#' @param lh life history list
#' @param true_years vector of true years (in case all_years and lc_years are 1:20 instead of 1998:2017)
#' @param True default=NULL, possible to specify true list from generated data if simulation
#' @param plot options for plotting include "Fish"=fishing mortality, "Rec"=recruitment (LIME only), "SPR"=spawning potential ratio, "ML"=mean length (including observed values; LIME only), "SB"=spawning biomass (LIME only), "Selex"=selectivity-at-length
#' @importFrom graphics points polygon segments
#' @importFrom grDevices rgb
#' 
#' @return figure with length composition data and model fits if Report or LBSPR are specified
#' 
#' @export
plot_output <- function(all_years, lc_years, Inputs=NULL, Report=NULL, Sdreport=NULL, LBSPR=NULL, lh, true_years, True=NULL, plot=c("Fish","Rec","SPR","ML","SB","Selex")){

    if(is.null(LBSPR)==FALSE){
        if(isS4(LBSPR)){
          LBSPR_outs <- list()
          LBSPR_outs$pLF <- LBSPR@pLCatch
          LBSPR_outs$SL50 <- LBSPR@Ests[,"SL50"]
          LBSPR_outs$SL95 <- LBSPR@Ests[,"SL95"]
          LBSPR_outs$FM <- LBSPR@Ests[,"FM"]
          LBSPR_outs$SPR <- LBSPR@Ests[,"SPR"]
          LBSPR_outs$SPR_Var <- LBSPR@Vars[,"SPR"]
          LBSPR_outs$SL50_Var <- LBSPR@Vars[,"SL50"]
          LBSPR_outs$SL95_Var <- LBSPR@Vars[,"SL95"]
          LBSPR_outs$FM_Var <- LBSPR@Vars[,"FM"]

          LBSPR <- LBSPR_outs
        }
    }

if(length(plot)==1) dim <- c(1,1)
if(length(plot)==2) dim <- c(2,1)
if(length(plot)==3) dim <- c(3,1)
if(length(plot)==4) dim <- c(2,2)
if(length(plot)==5 | length(plot)==6) dim <- c(2,3)

      by <- round(length(true_years)/5)
    lab <- rev(seq(from=true_years[length(true_years)], to=min(true_years), by=-by))
    ilab <- which(true_years %in% lab)

if(all(is.null(Inputs))==FALSE){
  if(Inputs$Data$n_s==1){
    xY <- seq_along(all_years)
    xLC <- which(all_years %in% lc_years)
  }
  if(Inputs$Data$n_s>1){
    xY <- 1:Inputs$Data$n_y
    xLC <- unique(Inputs$Data$S_yrs[which(all_years %in% lc_years)])
  }
  ilab2 <- sapply(1:length(ilab), function(x){
    sub <- which(Inputs$Data$S_yrs %in% ilab[x])
    return(sub[length(sub)])
  })
}
if(all(is.null(Inputs))){
  if(all(is.null(LBSPR))) stop("Must specify either LIME or LBSPR output")
  xY <- seq_along(all_years)
  xLC <- which(all_years %in% lc_years)
  ilab2 <- sapply(1:length(ilab), function(x){
    sub <- which(seq_along(all_years) %in% ilab[x])
    return(sub[length(sub)])
  })
}

par(mfrow=dim, mar=c(4,5,2,2))

if(all(is.null(Inputs))==FALSE){
  F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=lh$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.4)$root, error=function(e) NA)
  F30 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=lh$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, ref=0.3)$root, error=function(e) NA)
}

if("Fish" %in% plot){
  if(all(is.null(Sdreport))==FALSE){
    ylim <- c(0, max(Report$F_t)*1.1)
    if(all(is.na(Sdreport))==FALSE){
      sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_y"),]
      sd[,2][which(is.na(sd[,2]))] <- 0
      # ylim <- c(0, max(max(read_sdreport(sd, log=TRUE))*1.2))#, ymax))
    }
  }

  if(all(is.null(Report))==FALSE){
    plot(x=xY, y=Report$F_y, lwd=2, col="blue", ylim=ylim, type="l", xaxt="n", xaxs="i", yaxs="i", cex.axis=2, ylab="Fishing mortality", xlab="Year", cex.lab=2, xlim=c(min(xY),max(xY)+1))
    points(x=xLC, y=Report$F_y[xLC], col="blue", pch=19, xpd=NA, cex=2)
    if(all(is.na(Sdreport))==FALSE){
      polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)  
    }
  }

  if(all(is.null(Inputs))==FALSE){
    abline(h=F40*Inputs$Data$n_s, lwd=2, lty=2)
    abline(h=F30*Inputs$Data$n_s, lwd=2, lty=3)
  }
  if(all(is.null(True))==FALSE) lines(True$F_t, lwd=2)
  if(all(is.null(LBSPR))==FALSE & all(is.null(Report))==FALSE){
    par(new=TRUE)
    xxLC <- which(true_years %in% LBSPR$years)
    plot(x=xxLC, LBSPR$FM*(lh$M*lh$nseasons), xaxs="i", yaxs="i", xlab="", ylab="", ylim=ylim, xaxt="n", yaxt="n", lwd=2, col="red", type="p", pch=19, xpd=NA, cex=2, xlim=c(min(xY),max(xY)+1))
    index <- which(is.na(LBSPR$FM_Var)==FALSE)
    ignore <- sapply(1:length(xxLC), function(x) segments(x0=xxLC[x],x1=xxLC[x],y0=LBSPR$FM[index[x]]*(lh$M*lh$nseasons)-1.96*sqrt(LBSPR$FM_Var[index[x]]), y1=LBSPR$FM[index[x]]*(lh$M*lh$nseasons)+1.96*sqrt(LBSPR$FM_Var[index[x]]), lwd=4, col=rgb(1,0,0,alpha=0.4)))
    # lines(x=xxLC, y=LBSPR$Smooth[,"FM"]*(lh$M*lh$nseasons), lwd=2, col="red")
  }
  if(all(is.null(LBSPR))==FALSE & all(is.null(Report))){
    xxLC <- which(true_years %in% LBSPR$years)
    ylim <- c(0, max(LBSPR$FM*(lh$M*lh$nseasons))*1.5)
    plot(x=xxLC, LBSPR$FM*(lh$M*lh$nseasons), xaxs="i", yaxs="i", xlab="Year", ylab="Fishing mortality", cex.lab=2, ylim=ylim, xaxt="n", yaxt="n", lwd=2, col="red", type="p", pch=19, xpd=NA, cex=2, xlim=c(min(xY),max(xY)+1))
    index <- which(is.na(LBSPR$FM_Var)==FALSE)
    ignore <- sapply(1:length(xxLC), function(x) segments(x0=xxLC[x],x1=xxLC[x],y0=LBSPR$FM[index[x]]*(lh$M*lh$nseasons)-1.96*sqrt(LBSPR$FM_Var[index[x]]), y1=LBSPR$FM[index[x]]*(lh$M*lh$nseasons)+1.96*sqrt(LBSPR$FM_Var[index[x]]), lwd=4, col=rgb(1,0,0,alpha=0.4)))
    axis(2, cex.axis=2)
    # lines(x=xxLC, y=LBSPR$Smooth[,"FM"]*(lh$M*lh$nseasons), lwd=2, col="red")
  }
    axis(1, cex.axis=2, at=ilab, labels=lab)
}

if("Rec" %in% plot){
  if(all(is.na(Sdreport))==FALSE){
      sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]
      sd[,2][which(is.na(sd[,2]))] <- 0
  }
  ylim <- c(0, max(read_sdreport(sd, log=TRUE)))
  plot(x=seq_along(all_years), y=Report$R_t, lwd=2, col="blue", ylim=c(0, max(Report$R_t)*1.5), type="l", xaxt="n", ylab="Recruitment", xlab="Year", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2, xlim=c(min(seq_along(all_years)),max(seq_along(all_years))+lh$nseasons))
  points(x=which(all_years %in% lc_years), y=Report$R_t[which(all_years %in% lc_years)], col="blue", pch=19, xpd=NA, cex=2)
  if(all(is.na(Sdreport))==FALSE){
    polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
  }
  axis(1, cex.axis=2, at=ilab2, labels=lab)
  if(all(is.null(True))==FALSE){
    par(new=TRUE)
    plot(True$R_t, xaxs="i", yaxs="i", xlab="", ylab="", ylim=c(0, max(Report$R_t)*1.5), xaxt="n", yaxt="n", lwd=2, type="l", xlim=c(min(seq_along(all_years)),max(seq_along(all_years))+lh$nseasons))
  }
}

if("SPR" %in% plot){
  if(all(is.null(Report))==FALSE){
    plot(x=seq_along(all_years), y=Report$SPR_t, lwd=2, col="blue", ylim=c(0, 1), type="l", xaxt="n", ylab="SPR", xlab="Year", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2, xlim=c(min(seq_along(all_years)),max(seq_along(all_years))+lh$nseasons))
    points(x=which(all_years %in% lc_years), y=Report$SPR_t[which(all_years%in%lc_years)], col="blue", pch=19, xpd=NA, cex=2)
  }
  if(all(is.null(Sdreport))==FALSE){
    if(all(is.na(Sdreport))==FALSE){
    sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
    sd[,2][which(is.na(sd[,2]))] <- 0
    polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
    } 
  }
  if(all(is.null(True))==FALSE){
      par(new=TRUE)
    plot(True$SPR_t, xaxs="i", yaxs="i", xlab="", ylab="", ylim=c(0, 1), xaxt="n", yaxt="n", lwd=2, type="l", xlim=c(min(seq_along(all_years)),max(seq_along(all_years))+lh$nseasons))
  }
  if(all(is.null(LBSPR))==FALSE & all(is.null(Report))==FALSE){
    par(new=TRUE)
    xxLC <- which(true_years %in% LBSPR$years)
    plot(x=xxLC, LBSPR$SPR, xaxs="i", yaxs="i", xlab="", ylab="", ylim=c(0,1), xaxt="n", yaxt="n", lwd=2, col="red", type="p", pch=19, xpd=NA, cex=2, xlim=c(min(xY),max(xY)+1))
    index <- which(is.na(LBSPR$SPR_Var)==FALSE)
    ignore <- sapply(1:length(xxLC), function(x) segments(x0=xxLC[x],x1=xxLC[x],y0=LBSPR$SPR[index[x]]-1.96*sqrt(LBSPR$SPR_Var[index[x]]), y1=LBSPR$SPR[index[x]]+1.96*sqrt(LBSPR$SPR_Var[index[x]]), lwd=4, col=rgb(1,0,0,alpha=0.4)))
    lines(x=xxLC, y=LBSPR$Smooth[,"SPR"], lwd=2, col="red")
  }
  if(all(is.null(Report))){
    xxLC <- which(true_years %in% LBSPR$years)
    plot(x=xxLC, LBSPR$SPR, xaxs="i", yaxs="i", xlab="Year", ylab="SPR", ylim=c(0,1), xaxt="n", cex.axis=2, lwd=2, col="red", type="p", pch=19, xpd=NA, cex=2, xlim=c(min(xY),max(xY)+1), cex.lab=2)
    index <- which(is.na(LBSPR$SPR_Var)==FALSE)
    ignore <- sapply(1:length(xxLC), function(x) segments(x0=xxLC[x],x1=xxLC[x],y0=LBSPR$SPR[index[x]]-1.96*sqrt(LBSPR$SPR_Var[index[x]]), y1=LBSPR$SPR[index[x]]+1.96*sqrt(LBSPR$SPR_Var[index[x]]), lwd=4, col=rgb(1,0,0,alpha=0.4)))
    lines(x=xxLC, y=LBSPR$Smooth[,"SPR"], lwd=2, col="red")    
  }
      abline(h=0.4, lwd=2, lty=2)
    abline(h=0.3, lwd=2, lty=3)

    axis(1, cex.axis=2, at=ilab2, labels=lab)
}

if("ML" %in% plot){
  plot(x=seq_along(all_years), y=Report$L_t_hat, lwd=2, col="blue", ylim=c(0, max(Report$L_t_hat)*1.5), type="l", xaxt="n", ylab="Mean length", xlab="Year", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2, xlim=c(min(seq_along(all_years)), max(seq_along(all_years))+lh$nseasons))
  ML_obs <- sapply(1:nrow(Inputs$Data$LF), function(x) sum(Inputs$Data$LF[x,]*Inputs$Data$lbmids)/sum(Inputs$Data$LF[x,]))
  points(x=which(all_years %in% lc_years), y=ML_obs, pch=17, cex=2)
  points(x=which(all_years %in% lc_years), y=Report$L_t_hat[which(all_years %in% lc_years)], col="blue", pch=19, xpd=NA, cex=2)
  if(all(is.na(Sdreport))==FALSE){
  sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="L_t_hat"),]
  sd[,2][which(is.na(sd[,2]))] <- 0
    polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
  }
  axis(1, cex.axis=2, at=ilab2, labels=lab)
  if(all(is.null(True))==FALSE){
    par(new=TRUE)
    plot(True$ML_t, , xaxs="i", yaxs="i", xlab="", ylab="", ylim=c(0, max(Report$L_t_hat)*1.5), xaxt="n", yaxt="n", lwd=2, type="l")
  }
}

if("SB" %in% plot){
  if(all(is.na(Sdreport))==FALSE){
      sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lD_t"),]
      sd[,2][which(is.na(sd[,2]))] <- 0
  }
  ylim <- c(0, max(read_sdreport(sd, log=TRUE)))
  plot(x=seq_along(all_years), y=Report$D_t, lwd=2, col="blue", ylim=c(0, max(Report$D_t)*1.5), type="l", xaxt="n", ylab="Relative spawning biomass", xlab="Year", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2, xlim=c(min(seq_along(all_years)), max(seq_along(all_years))+lh$nseasons))
  points(x=which(all_years %in% lc_years), y=Report$D_t[which(all_years %in% lc_years)], col="blue", pch=19, xpd=NA, cex=2)
  if(all(is.na(Sdreport))==FALSE){
    polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
  }
  axis(1, cex.axis=2, at=ilab2, labels=lab)
  if(all(is.null(True))==FALSE){
      par(new=TRUE)
    plot(True$D_t, , xaxs="i", yaxs="i", xlab="", ylab="", ylim=c(0, max(Report$D_t)*1.5), xaxt="n", yaxt="n", lwd=2, type="l")
  }
}
    
if("Selex" %in% plot){


    xlabs <- pretty(seq_along(lh$highs))
    plot_labs <- rep(NA, length(xlabs))
    if(xlabs[1]!=0) warning("Should start length bin labels at 0")
    plot_labs[1] <- 0
    elabs <- as.numeric(lh$highs[xlabs][which(is.na(lh$highs[xlabs])==FALSE)])
    plot_labs[2:(length(elabs)+1)] <- elabs
    if(is.na(plot_labs[length(plot_labs)])) plot_labs[length(plot_labs)] <- max(elabs) + elabs[1]



  mids <- lh$mids
  if(all(is.null(Report))==FALSE){
    plot(x=1:length(mids), y=Report$S_l, lwd=2, col="blue", ylim=c(0, 1.1), type="l", xaxt="n", ylab="Selectivity at length", xlab="Length (cm)", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2)
    if(all(is.na(Sdreport))==FALSE)  polygon( y=read_sdreport(summary(Sdreport)[which(rownames(summary(Sdreport))=="S_l"),], log=FALSE), x=c(which(is.na(summary(Sdreport)[which(rownames(summary(Sdreport))=="S_l"),2])==FALSE), rev(which(is.na(summary(Sdreport)[which(rownames(summary(Sdreport))=="S_l"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
  }
  if(all(is.null(LBSPR))==FALSE & all(is.null(Report))==FALSE){
    for(i in 1:length(lc_years)){
      SL50 <- LBSPR$SL50[i]
      SL95 <- LBSPR$SL95[i]
      S_l2 <- 1.0/(1+exp(-log(19)*(mids-SL50)/(SL95-SL50))) # Selectivity-at-Length
      lines(x=1:length(mids), y=S_l2, col="#AA000050", lwd=2)
    }
  legend("bottomright", col=c("blue", "red", "black", "black","black"), lwd=2, legend=c("LIME", "LB-SPR", "SPR 40%", "SPR 30%", "Observed"), cex=1.7, lty=c(1,1,2,3,0), pch=c(19,19,NA,NA,17))
  }
  if(all(is.null(LBSPR))==FALSE & all(is.null(Report))){
    for(i in 1:length(lc_years)){
      SL50 <- LBSPR$SL50[i]
      SL95 <- LBSPR$SL95[i]
      S_l2 <- 1.0/(1+exp(-log(19)*(mids-SL50)/(SL95-SL50))) # Selectivity-at-Length
      if(i==1) plot(x=1:length(mids), y=S_l2, col="#AA000050", lwd=2, ylim=c(0,1.1), type="l", xaxt="n", ylab="Selectivity at length", xlab="Length (cm)", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2)
      if(i>1) lines(x=1:length(mids), y=S_l2, col="#AA000050", lwd=2)
      legend("bottomright", col=c("red", "black", "black","black"), lwd=2, legend=c("LB-SPR", "SPR 40%", "SPR 30%"), cex=1.7, lty=c(1,2,3), pch=c(19,NA,NA))
    }
  }
  if(all(is.null(True))==FALSE) lines(True$S_l, lwd=2)
    axis(1, cex.axis=2, at=xlabs, labels=plot_labs)
}

}