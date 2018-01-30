#' plot LIME output
#'
#' \code{plot_output} plot output from LIME or LB-SPR
#'
#' @author M.B. Rudd
#' @param Inputs LIME input file
#' @param Report LIME report file
#' @param Sdreport LIME standard error file
#' @param LBSPR LBSPR results - either straight from LBSPR output or LIME-edited
#' @param lh life history list
#' @param true_years vector of true years (in case all_years and lc_years are 1:20 instead of 1998:2017)
#' @param True default=NULL, possible to specify true list from generated data if simulation
#' @param plot options for plotting include "Fish"=fishing mortality, "Rec"=recruitment (LIME only), "SPR"=spawning potential ratio, "ML"=mean length (including observed values; LIME only), "SB"=spawning biomass (LIME only), "Selex"=selectivity-at-length
#' @param set_ylim list to set ylim for each of the 'plot' arguments
#' @importFrom graphics points polygon segments
#' @importFrom grDevices rgb
#' 
#' @return figure with length composition data and model fits if Report or LBSPR are specified
#' 
#' @export
plot_output <- function(Inputs=NULL, Report=NULL, Sdreport=NULL, LBSPR=NULL, lh, true_years=NULL, True=NULL, plot=c("Fish","Rec","SPR","ML","SB","Selex"), set_ylim=list("Fish" = c(0,1), "SPR" = c(0,1))){

  all_years <- Inputs$Data$T_yrs
  lc_years <- Inputs$Data$LC_yrs
  if(all(is.null(true_years))) true_years <- all_years

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
          LBSPR_outs$years <- LBSPR@Years

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

    # nf <- Inputs$Data$n_f
    nf <- 1
    ns <- Inputs$Data$n_s

if(all(is.null(Inputs))==FALSE){
  if(ns==1){
    xY <- seq_along(all_years)
    xLC <- which(all_years %in% lc_years)
  }
  if(ns>1){
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
  F40 <- tryCatch(sapply(1:nf, function(x) uniroot(calc_ref, lower=0, upper=200, ages=lh$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_fa[x,], ref=0.4)$root), error=function(e) NA)
  F30 <- tryCatch(sapply(1:nf, function(x) uniroot(calc_ref, lower=0, upper=200, ages=lh$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_fa[x,], ref=0.3)$root), error=function(e) NA)
}

if("Fish" %in% plot){
    if("Fish" %in% names(set_ylim) ==FALSE) ylim <- c(0, max(Report$F_y)*2)
    if("Fish" %in% names(set_ylim)) ylim <- set_ylim[["Fish"]]

  if(all(is.null(Sdreport))==FALSE){
    if(all(is.na(Sdreport))==FALSE){
      sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_y"),]
      sd[,2][which(is.na(sd[,2]))] <- 0
      # ylim <- c(0, max(max(read_sdreport(sd, log=TRUE))*1.2))#, ymax))
    }
  }

  if(all(is.null(Report))==FALSE){
    if(nf>1){
      colfn <- colorRampPalette(c("red","blue"))
      cols <- colfn(nf)
    }
    if(nf==1) cols <- "#228B22"
    for(f in 1:nf){
      if(f==1) plot(x=xY, y=Report$F_y, lwd=2, col=cols[f], ylim=ylim, type="l", xaxt="n", xaxs="i", yaxs="i", cex.axis=2, ylab="Fishing mortality", xlab="Year", cex.lab=2, xlim=c(min(xY),max(xY)))
      if(f>1) lines(x=xY, y=Report$F_y, lwd=2, col=cols[f])
      points(x=xLC, y=Report$F_y[xLC], col=cols[f], pch=19, cex=2, xpd=NA)
      index <- seq(f, nrow(sd), by=nf)
      if(all(is.na(Sdreport))==FALSE){
        polygon( y=read_sdreport(sd[index,], log=TRUE), x=c(which(is.na(sd[index,2])==FALSE), rev(which(is.na(sd[index,2])==FALSE))), col=paste0(cols[f],"20"), border=NA)  
      }
    }
  }

  if(all(is.null(Inputs))==FALSE){
    if(nf>1){
      colfn <- colorRampPalette(c("red","blue"))
      cols <- colfn(nf)
    }
    if(nf==1) cols <- "#228B22"    
    makelines <- sapply(1:nf, function(x){
      abline(h=F40[x]*ns, lwd=2, lty=2, col=cols[x])
      abline(h=F30[x]*ns, lwd=2, lty=3, col=cols[x])
    })
  }

  if(all(is.null(True))==FALSE) lines(True$F_t, lwd=2)
  if(all(is.null(LBSPR))==FALSE & all(is.null(Report))==FALSE){
    par(new=TRUE)
    xxLC <- which(true_years %in% LBSPR$years)
    plot(x=xxLC, LBSPR$FM*(lh$M*lh$nseasons), xaxs="i", yaxs="i", xlab="", ylab="", ylim=ylim, xaxt="n", yaxt="n", lwd=2, col=gray(0.3), type="p", pch=19, cex=2, xlim=c(min(xY),max(xY)))
    index <- which(is.na(LBSPR$FM_Var)==FALSE)
    ignore <- sapply(1:length(xxLC), function(x) segments(x0=xxLC[x],x1=xxLC[x],y0=LBSPR$FM[index[x]]*(lh$M*lh$nseasons)-1.96*sqrt(LBSPR$FM_Var[index[x]]), y1=LBSPR$FM[index[x]]*(lh$M*lh$nseasons)+1.96*sqrt(LBSPR$FM_Var[index[x]]), lwd=4, col=paste0(gray(0.3),"40")))
    # lines(x=xxLC, y=LBSPR$Smooth[,"FM"]*(lh$M*lh$nseasons), lwd=2, col=gray(0.3))
  }
  if(all(is.null(LBSPR))==FALSE & all(is.null(Report))){
    xxLC <- which(true_years %in% LBSPR$years)
    ylim <- c(0, max(LBSPR$FM*(lh$M*lh$nseasons))*1.5)
    plot(x=xxLC, LBSPR$FM*(lh$M*lh$nseasons), xaxs="i", yaxs="i", xlab="Year", ylab="Fishing mortality", cex.lab=2, ylim=ylim, xaxt="n", yaxt="n", lwd=2, col=gray(0.3), type="p", pch=19, cex=2, xlim=c(min(xY),max(xY)))
    index <- which(is.na(LBSPR$FM_Var)==FALSE)
    ignore <- sapply(1:length(xxLC), function(x) segments(x0=xxLC[x],x1=xxLC[x],y0=LBSPR$FM[index[x]]*(lh$M*lh$nseasons)-1.96*sqrt(LBSPR$FM_Var[index[x]]), y1=LBSPR$FM[index[x]]*(lh$M*lh$nseasons)+1.96*sqrt(LBSPR$FM_Var[index[x]]), lwd=4, col=paste0(gray(0.3),"40")))
    axis(2, cex.axis=2)
    # lines(x=xxLC, y=LBSPR$Smooth[,"FM"]*(lh$M*lh$nseasons), lwd=2, col=gray(0.3))
  }
    axis(1, cex.axis=2, at=ilab, labels=lab)
}

if("Rec" %in% plot){
  if(all(is.na(Sdreport))==FALSE){
      sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]
      sd[,2][which(is.na(sd[,2]))] <- 0
  }
  if("Rec" %in% names(set_ylim) == FALSE) ylim <- c(0, max(read_sdreport(sd, log=TRUE)))
  if("Rec" %in% names(set_ylim)) ylim <- set_ylim[["Rec"]]

  plot(x=seq_along(all_years), y=Report$R_t, lwd=2, col="#228B22", ylim=c(0, max(Report$R_t)*1.5), type="l", xaxt="n", ylab="Recruitment", xlab="Year", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2, xlim=c(min(seq_along(all_years)),max(seq_along(all_years))), xpd=NA)
  points(x=which(all_years %in% lc_years), y=Report$R_t[which(all_years %in% lc_years)], col="#228B22", pch=19, cex=2, xpd=NA)
  if(all(is.na(Sdreport))==FALSE){
    polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col="#228B2220", border=NA)
  }
  axis(1, cex.axis=2, at=ilab2, labels=lab)
  if(all(is.null(True))==FALSE) lines(True$R_t, lwd=2)

}

if("SPR" %in% plot){
  if(all(is.null(Report))==FALSE){
  if("SPR" %in% names(set_ylim) == FALSE) ylim <- c(0,1)
  if("SPR" %in% names(set_ylim)) ylim <- set_ylim[["SPR"]]
    plot(x=seq_along(all_years), y=Report$SPR_t, lwd=2, col="#228B22", ylim=c(0, 1), type="l", xaxt="n", ylab="SPR", xlab="Year", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2, xlim=c(min(seq_along(all_years)),max(seq_along(all_years))), xpd=NA)
    points(x=which(all_years %in% lc_years), y=Report$SPR_t[which(all_years%in%lc_years)], col="#228B22", pch=19, cex=2, xpd=NA)
  }
  if(all(is.null(Sdreport))==FALSE){
    if(all(is.na(Sdreport))==FALSE){
    sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
    sd[,2][which(is.na(sd[,2]))] <- 0
    polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col="#228B2220", border=NA)
    } 
  }
  if(all(is.null(True))==FALSE) lines(True$SPR_t, lwd=2)
  if(all(is.null(LBSPR))==FALSE & all(is.null(Report))==FALSE){
    par(new=TRUE)
    xxLC <- which(true_years %in% LBSPR$years)
    plot(x=xxLC, LBSPR$SPR, xaxs="i", yaxs="i", xlab="", ylab="", ylim=c(0,1), xaxt="n", yaxt="n", lwd=2, col=gray(0.3), type="p", pch=19, cex=2, xlim=c(min(xY),max(xY)))
    index <- which(is.na(LBSPR$SPR_Var)==FALSE)
    ignore <- sapply(1:length(xxLC), function(x) segments(x0=xxLC[x],x1=xxLC[x],y0=LBSPR$SPR[index[x]]-1.96*sqrt(LBSPR$SPR_Var[index[x]]), y1=LBSPR$SPR[index[x]]+1.96*sqrt(LBSPR$SPR_Var[index[x]]), lwd=4, col=paste0(gray(0.3),"40")))
    lines(x=xxLC, y=LBSPR$Smooth[,"SPR"], lwd=2, col=gray(0.3))
  }
  if(all(is.null(Report))){
    xxLC <- which(true_years %in% LBSPR$years)
    plot(x=xxLC, LBSPR$SPR, xaxs="i", yaxs="i", xlab="Year", ylab="SPR", ylim=c(0,1), xaxt="n", cex.axis=2, lwd=2, col=gray(0.3), type="p", pch=19, cex=2, xlim=c(min(xY),max(xY)), cex.lab=2)
    index <- which(is.na(LBSPR$SPR_Var)==FALSE)
    ignore <- sapply(1:length(xxLC), function(x) segments(x0=xxLC[x],x1=xxLC[x],y0=LBSPR$SPR[index[x]]-1.96*sqrt(LBSPR$SPR_Var[index[x]]), y1=LBSPR$SPR[index[x]]+1.96*sqrt(LBSPR$SPR_Var[index[x]]), lwd=4, col=paste0(gray(0.3),"40")))
    lines(x=xxLC, y=LBSPR$Smooth[,"SPR"], lwd=2, col=gray(0.3))    
  }
      abline(h=0.4, lwd=2, lty=2)
    abline(h=0.3, lwd=2, lty=3)

    axis(1, cex.axis=2, at=ilab2, labels=lab)
}

if("ML" %in% plot){

  if("ML" %in% names(set_ylim) == FALSE) ylim <- c(0, max(Report$L_t_hat)*1.5)
  if("ML" %in% names(set_ylim)) ylim <- set_ylim[["ML"]]

  if(all(is.na(Sdreport))==FALSE){
    sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="L_t_hat"),]
    sd[,2][which(is.na(sd[,2]))] <- 0
  }
    if(nf>1){
      colfn <- colorRampPalette(c("red","blue"))
      cols <- colfn(nf)
    }
    if(nf==1) cols <- "#228B22"
  for(f in 1:nf){
    ML_obs <- sapply(1:nrow(Inputs$Data$LF_tbf[,,f]), function(x) sum(Inputs$Data$LF_tbf[x,,f]*Inputs$Data$lbmids)/sum(Inputs$Data$LF_tbf[x,,f]))

    if(f==1)   plot(x=seq_along(all_years), y=Report$L_t_hat, lwd=2, col=cols[f], ylim=ylim, type="l", xaxt="n", ylab="Mean length", xlab="Year", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2, xlim=c(min(seq_along(all_years)), max(seq_along(all_years))))
    if(f>1) lines(x=seq_along(all_years), y=Report$L_t_hat, lwd=2, col=cols[f])
    points(x=which(all_years %in% lc_years), y=Report$L_t_hat[which(all_years %in% lc_years)], col=cols[f], pch=19, cex=2, xpd=NA)
    points(x=which(all_years %in% lc_years), y=ML_obs, cex=2.5, xpd=NA, col=cols[f], lwd=2)
    index <- seq(f,nrow(sd), by=nf)
    polygon( y=read_sdreport(sd[index,], log=FALSE), x=c(which(is.na(sd[index,2])==FALSE), rev(which(is.na(sd[index,2])==FALSE))), col=paste0(cols[f],"20"), border=NA)
  }



  axis(1, cex.axis=2, at=ilab2, labels=lab)
  if(all(is.null(True))==FALSE) lines(True$ML_t, lwd=2)

}

if("SB" %in% plot){
  if(all(is.na(Sdreport))==FALSE){
      sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lD_t"),]
      sd[,2][which(is.na(sd[,2]))] <- 0
  }
  if("SB" %in% names(set_ylim) == FALSE) ylim <- c(0, max(Report$D_t)*1.5)
  if("SB" %in% names(set_ylim)) ylim <- set_ylim[["SB"]]

  plot(x=seq_along(all_years), y=Report$D_t, lwd=2, col="#228B22", ylim=ylim, type="l", xaxt="n", ylab="Relative spawning biomass", xlab="Year", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2, xlim=c(min(seq_along(all_years)), max(seq_along(all_years))+lh$nseasons))
  points(x=which(all_years %in% lc_years), y=Report$D_t[which(all_years %in% lc_years)], col="#228B22", pch=19, cex=2, xpd=NA)
  if(all(is.na(Sdreport))==FALSE){
    polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col="#228B2220", border=NA)
  }
  axis(1, cex.axis=2, at=ilab2, labels=lab)
  if(all(is.null(True))==FALSE) lines(True$D_t, lwd=2)

}
    
if("Selex" %in% plot){


    xlabs <- pretty(seq_along(lh$highs))
    plot_labs <- rep(NA, length(xlabs))
    if(xlabs[1]!=0) warning("Should start length bin labels at 0")
    plot_labs[1] <- 0
    elabs <- as.numeric(lh$highs[xlabs][which(is.na(lh$highs[xlabs])==FALSE)])
    plot_labs[2:(length(elabs)+1)] <- elabs
    if(is.na(plot_labs[length(plot_labs)])) plot_labs[length(plot_labs)] <- max(elabs) + elabs[1]


  if(all(is.null(Sdreport))==FALSE){
    if(all(is.na(Sdreport))==FALSE){
      sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="S_fl"),]
      sd[,2][which(is.na(sd[,2]))] <- 0
      # ylim <- c(0, max(max(read_sdreport(sd, log=TRUE))*1.2))#, ymax))
    }
  }

  mids <- lh$mids
  if(all(is.null(Report))==FALSE){
    if(nf>1){
      colfn <- colorRampPalette(c("red","blue"))
      cols <- colfn(nf)
    }
    if(nf==1) cols <- "#228B22"

    for(f in 1:nf){
      if(f==1) plot(x=1:length(mids), y=Report$S_fl[f,], lwd=2, col=cols[f], ylim=c(0, 1.1), type="l", xaxt="n", ylab="Selectivity at length", xlab="Length (cm)", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2)
      if(f>1) lines(x=1:length(mids), y=Report$S_fl[f,], lwd=2, col=cols[f])
      if(all(is.na(Sdreport))==FALSE){
        index <- seq(f,nrow(sd),by=nf)
        polygon( y=read_sdreport(sd[index,], log=FALSE), x=c(which(is.na(sd[index,2])==FALSE), rev(which(is.na(sd[index,2])==FALSE))), col=paste0(cols[f],"20"), border=NA)  
      }
    }
  }
  if(all(is.null(LBSPR))==FALSE & all(is.null(Report))==FALSE){
    for(i in 1:length(lc_years)){
      SL50 <- LBSPR$SL50[i]
      SL95 <- LBSPR$SL95[i]
      S_l2 <- 1.0/(1+exp(-log(19)*(mids-SL50)/(SL95-SL50))) # Selectivity-at-Length
      lines(x=1:length(mids), y=S_l2, col=paste0(gray(0.3),"50"), lwd=2)
    }
  # legend("bottomright", col=c("#228B22", gray(0.3), "black", "black","black"), lwd=2, legend=c("LIME", "LB-SPR", "SPR 40%", "SPR 30%", "Observed"), cex=1.7, lty=c(1,1,2,3,0), pch=c(19,19,NA,NA,17))
  }
  if(all(is.null(LBSPR))==FALSE & all(is.null(Report))){
    for(i in 1:length(lc_years)){
      SL50 <- LBSPR$SL50[i]
      SL95 <- LBSPR$SL95[i]
      S_l2 <- 1.0/(1+exp(-log(19)*(mids-SL50)/(SL95-SL50))) # Selectivity-at-Length
      if(i==1) plot(x=1:length(mids), y=S_l2, col=paste0(gray(0.3),"50"), lwd=2, ylim=c(0,1.1), type="l", xaxt="n", ylab="Selectivity at length", xlab="Length (cm)", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2)
      if(i>1) lines(x=1:length(mids), y=S_l2, col=paste0(gray(0.3),"50"), lwd=2)
      # legend("bottomright", col=c(gray(0.3), "black", "black","black"), lwd=2, legend=c("LB-SPR", "SPR 40%", "SPR 30%"), cex=1.7, lty=c(1,2,3), pch=c(19,NA,NA))
    }
  }
  if(all(is.null(True))==FALSE) lines(True$S_l[1,], lwd=2)
    axis(1, cex.axis=2, at=xlabs, labels=plot_labs)
}

if("Fish" %in% plot | "Selex" %in% plot | "ML" %in% plot & nf > 1){
  legend("bottomright", col=c("#228B22", cols), legend=c("Total", sapply(1:nf, function(x) paste0("Fleet ", x))), lty=1, lwd=2)
}

}