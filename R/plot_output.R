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
plot_output <- function(Inputs=NULL, Report=NULL, Sdreport=NULL, LBSPR=NULL, lh, true_years=NULL, True=NULL, plot=c("Fish","Rec","SPR","ML","SB","Selex"), set_ylim=list("SPR" = c(0,1)), legend=FALSE){

    nf <- Inputs$Data$n_f
    ns <- Inputs$Data$n_s

  all_years <- 1:Inputs$Data$n_t
  lc_years <- lapply(1:nf, function(x){
    sub <- Inputs$Data$LF_tlf[,,x]
    if(is.vector(sub) & sum(sub)>0) out <- all_years
    if(is.matrix(sub)) out <- which(rowSums(sub) > 0)
    return(out)
  })
  if(all(is.null(true_years))) true_years <- all_years

    if(is.null(LBSPR)==FALSE){
        if(isS4(LBSPR)){
          LBSPR_outs <- list()
          LBSPR_outs$pLF <- LBSPR@pLCatch
          LBSPR_outs$SL50 <- LBSPR@Ests[,"SL50"]
          LBSPR_outs$SL95 <- LBSPR@Ests[,"SL95"]
          LBSPR_outs$FM<-LBSPR@FM
          LBSPR_outs$SPR<-LBSPR@SPR
          LBSPR_outs$FM_smooth <- LBSPR@Ests[,"FM"]
          LBSPR_outs$SPR_smooth <- LBSPR@Ests[,"SPR"]
          LBSPR_outs$var_FM <- LBSPR@Vars[,"FM"] 
          LBSPR_outs$var_SPR <- LBSPR@Vars[,"SPR"]
          LBSPR_outs$var_S50 <- LBSPR@Vars[,"SL50"]
          LBSPR_outs$var_S95 <- LBSPR@Vars[,"SL95"]
          LBSPR_outs$years <- LBSPR@Years

          xLC_lbspr <- which(true_years %in% LBSPR@Years)
          if(length(xLC_lbspr)==0) xLC_lbspr <- which(seq_along(true_years) %in% LBSPR@Years)

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
  if(ns==1){
    xY <- seq_along(all_years)
    xLC <- lapply(1:nf, function(x) which(all_years %in% lc_years[[x]]))
  }
  if(ns>1){
    xY <- 1:Inputs$Data$n_y
    xLC <- lapply(1:nf, function(x) unique(Inputs$Data$S_yrs[which(all_years %in% lc_years[[x]])]))
  }
  ilab2 <- sapply(1:length(ilab), function(x){
    sub <- which(Inputs$Data$S_yrs %in% ilab[x])
    return(sub[length(sub)])
  })
}
if(all(is.null(Inputs))){
  if(all(is.null(LBSPR))) stop("Must specify either LIME or LBSPR output")
  xY <- seq_along(all_years)
  xLC <- lapply(1:nf, function(x) which(all_years %in% lc_years[[x]]))
  ilab2 <- sapply(1:length(ilab), function(x){
    sub <- which(seq_along(all_years) %in% ilab[x])
    return(sub[length(sub)])
  })
}

par(mfrow=dim, mar=c(4,5,2,2))

if(all(is.null(Inputs))==FALSE){
  F50 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=lh$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_fa=Report$S_fa, ref=0.5)$root, error=function(e) NA)
  F30 <- tryCatch(uniroot(calc_ref, lower=0, upper=200, ages=lh$ages, Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_fa=Report$S_fa, ref=0.3)$root, error=function(e) NA)
}
    col_total <- "#228B22"
    if(nf>1){
      colfn <- colorRampPalette(c("red","blue"))
      cols <- colfn(nf)
    }
    if(nf==1) cols <- col_total


if("Fish" %in% plot){
    if("Fish" %in% names(set_ylim) ==FALSE) ylim <- c(0, max(Report$F_y)*2)
    if("Fish" %in% names(set_ylim)) ylim <- set_ylim[["Fish"]]

  if(all(is.null(Sdreport))==FALSE){
    if(all(is.na(Sdreport))==FALSE){
      sd_total <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_y"),]
      sd_fleet <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_fy"),]
      if(is.vector(sd_fleet)) sd_fleet <- as.matrix(sd_fleet)
      # sd[,2][which(is.na(sd[,2]))] <- 0
      # ylim <- c(0, max(max(read_sdreport(sd, log=TRUE))*1.2))#, ymax))
    }
  }

  if(all(is.null(Report))==FALSE){
    plot(x=1, y=1, type="n", xaxt="n", xaxs="i", yaxs="i", ylim=ylim, cex.axis=2, ylab="Fishing mortality", xlab="Year", cex.lab=2, xlim=c(min(xY),max(xY)))

      # index <- seq(f, nrow(sd), by=nf)
      if(all(is.na(Sdreport))==FALSE){
        yvals <- read_sdreport(sd_total, log=TRUE)
        yvals[which(is.finite(yvals)==FALSE)] <- max(yvals[-which(is.finite(yvals)==FALSE)])
        polygon( y=yvals, x=c(which(is.na(sd_total[,2])==FALSE), rev(which(is.na(sd_total[,2])==FALSE))), col=paste0(col_total,"40"), border=NA)  
      }

    for(f in 1:nf){
      if(is.null(nrow(sd_fleet))==FALSE){
        if(nrow(sd_fleet)==length(xY)) index <- xY
        if(nrow(sd_fleet) > length(xY)) index <- seq(from=f, to=nrow(sd_fleet), by=nf)
        lines(x=xY, y=exp(sd_fleet[index,1]), lwd=2, col=cols[f])
        points(x=xLC[[f]], y=exp(sd_fleet[index[xLC[[f]]],1]), col=cols[f], pch=19, cex=2, xpd=NA)
      }
    }
    lines(x=xY, y=exp(sd_total[,1]), lwd=2, col=col_total)
    points(x=unique(unlist(xLC)), y=exp(sd_total[unique(unlist(xLC)),1]), lwd=2, col=col_total, pch=19, cex=2, xpd=NA)
  }

  if(all(is.null(Inputs))==FALSE){  
    makelines <- sapply(1:nf, function(x){
      abline(h=F50[x]*ns, lwd=2, lty=2)
      abline(h=F30[x]*ns, lwd=2, lty=3)
    })
  }

  if(all(is.null(True))==FALSE) lines(True$F_y, lwd=2)
  if(all(is.null(LBSPR))==FALSE & all(is.null(Report))==FALSE){
    points(x=xLC_lbspr, LBSPR$FM*(lh$M*lh$nseasons),col="#AA00AA", pch=19, cex=2, xpd=NA)
    index <- which(is.na(LBSPR$var_FM)==FALSE)
    ignore <- sapply(1:length(xLC_lbspr), function(x) segments(x0=xLC_lbspr[x],x1=xLC_lbspr[x],y0=LBSPR$FM[index[x]]*(lh$M*lh$nseasons)-1.96*sqrt(LBSPR$var_FM[index[x]]), y1=LBSPR$FM[index[x]]*(lh$M*lh$nseasons)+1.96*sqrt(LBSPR$var_FM[index[x]]), lwd=4, col=paste0("#AA00AA","40")))
    lines(x=xLC_lbspr, y=LBSPR$FM_smooth*(lh$M*lh$nseasons), col="#AA00AA", lw=2)
  }
}

if("Rec" %in% plot){

  if(all(is.na(Sdreport))==FALSE){
      sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]
      sd[,2][which(is.na(sd[,2]))] <- 0
      sd <- sd[seq(1,by=ns,length.out=Inputs$Data$n_y),]
      r_est <- exp(sd[,1])

      if("Rec" %in% names(set_ylim) == FALSE) ylim <- c(0, max(read_sdreport(sd, log=TRUE)))
      if("Rec" %in% names(set_ylim)) ylim <- set_ylim[["Rec"]]
      
      plot(x=1, y=1, type="n", xaxt="n", xaxs="i", yaxs="i", ylim=ylim, cex.axis=2, ylab="Recruitment", xlab="Year", cex.lab=2, xlim=c(min(xY),max(xY)))    
        polygon( y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_total, "40"), border=NA)
        lines(x=seq_along(xY), y=r_est, lwd=2, col=col_total, ylim=ylim, xpd=NA)
        points(x=unique(unlist(xLC)), y=r_est[unique(unlist(xLC))], col=col_total, pch=19, cex=2, xpd=NA)
        axis(1, cex.axis=2, at=ilab2, labels=lab)
        if(all(is.null(True))==FALSE) lines(True$R_t[seq(1,by=ns,length.out=Inputs$Data$n_y)], lwd=2)
  }

}

if("SPR" %in% plot){

  if("SPR" %in% names(set_ylim) == FALSE) ylim <- c(0,1)
  if("SPR" %in% names(set_ylim)) ylim <- set_ylim[["SPR"]]
  plot(x=1, y=1, type="n", xaxt="n", ylab="SPR", xlab="Year", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2, xlim=c(min(seq_along(xY)), max(seq_along(xY))), ylim=ylim)

  if(all(is.null(Sdreport))==FALSE){
    if(all(is.na(Sdreport))==FALSE){
      sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),]
      sd[,2][which(is.na(sd[,2]))] <- 0
      sd <- sd[seq(1,by=ns,length.out=Inputs$Data$n_y),]  
      polygon( y=read_sdreport(sd, log=FALSE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_total, "40"), border=NA)
      lines(x=seq_along(xY), y=sd[,1], lwd=2, col=col_total)
      points(x=unique(unlist(xLC)), y=sd[unique(unlist(xLC)),1], pch=19, cex=2, xpd=NA, col=col_total)
    } 
  }
  if(all(is.null(True))==FALSE) lines(True$SPR_t[seq(1,by=ns,length.out=Inputs$Data$n_y)], lwd=2)
  if(all(is.null(LBSPR))==FALSE){
    points(x=xLC_lbspr, LBSPR$SPR, col="#AA00AA", pch=19, cex=2, xpd=NA)
    index <- which(is.na(LBSPR$var_SPR)==FALSE)
    ignore <- sapply(1:length(xLC_lbspr), function(x) segments(x0=xLC_lbspr[x],x1=xLC_lbspr[x],y0=LBSPR$SPR[index[x]]-1.96*sqrt(LBSPR$var_SPR[index[x]]), y1=LBSPR$SPR[index[x]]+1.96*sqrt(LBSPR$var_SPR[index[x]]), lwd=4, col=paste0("#AA00AA","40")))
    lines(x=xLC_lbspr, y=LBSPR$SPR_smooth, lwd=2, col="#AA00AA")
  }
      abline(h=0.4, lwd=2, lty=2)
    abline(h=0.3, lwd=2, lty=3)

    axis(1, cex.axis=2, at=ilab2, labels=lab)
}

if("ML" %in% plot){

  if("ML" %in% names(set_ylim) == FALSE) ylim <- c(0, max(Report$ML_ft_hat)*1.5)
  if("ML" %in% names(set_ylim)) ylim <- set_ylim[["ML"]]
  
  plot(x=1,y=1, type="n", xaxt="n", ylab="Mean length", xlab="Year", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2, xlim=c(min(seq_along(xY)), max(seq_along(xY))), ylim=ylim)

  if(all(is.na(Sdreport))==FALSE){
    sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="ML_ft_hat"),]
    sd[,2][which(is.na(sd[,2]))] <- 0
    for(f in 1:nf){
      # sd <- sd[seq(f,by=nf,length.out=Inputs$Data$n_y),]  
      sdf <- sd[seq(f*ns,by=nf*ns,length.out=Inputs$Data$n_y),]  

      ML_obs <- lapply(1:nf, function(y){
        subLF <- Inputs$Data$LF_tlf[,,y]
        ml <- sapply(1:nrow(subLF), function(x){
          sum(subLF[x,]*Inputs$Data$lbmids)/sum(subLF[x,])
        })
        return(ml[seq(1,by=ns,length.out=Inputs$Data$n_y)])
      })
      polygon(y=read_sdreport(sdf, log=FALSE), x=c(which(is.na(sdf[,2])==FALSE), rev(which(is.na(sdf[,2])==FALSE))), col=paste0(cols[f],"40"), border=NA)
      lines(x=seq_along(xY), y=sdf[,1], lwd=2, col=cols[f])
      points(x=which(is.na(ML_obs[[f]])==FALSE), y=sdf[which(is.na(ML_obs[[f]])==FALSE),1], pch=19, col=cols[f], xpd=NA, cex=2)
      lines(x=seq_along(xY), y=ML_obs[[f]], lwd=2, xpd=NA)
    }
  }

  axis(1, cex.axis=2, at=ilab2, labels=lab)

}

if("SB" %in% plot){
  if("SB" %in% names(set_ylim) == FALSE) ylim <- c(0, max(Report$D_t)*1.5)
  if("SB" %in% names(set_ylim)) ylim <- set_ylim[["SB"]]

  plot(x=1, y=1, type="n", xaxt="n", ylab="Relative spawning biomass", xlab="Year", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2, xlim=c(min(seq_along(xY)), max(seq_along(xY))), ylim=ylim)

  if(all(is.na(Sdreport))==FALSE){
      sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lD_t"),]
      sd[,2][which(is.na(sd[,2]))] <- 0
      sd <- sd[seq(1,by=ns,length.out=Inputs$Data$n_y),]  

      polygon(y=read_sdreport(sd, log=TRUE), x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_total, "40"), border=NA)
      lines(x=seq_along(xY), y=exp(sd[,1]), lwd=2, col=col_total)
      points(x=unique(unlist(xLC)), y=exp(sd[unique(unlist(xLC)),1]), pch=19, col=col_total, xpd=NA, cex=2)

  }
  axis(1, cex.axis=2, at=ilab2, labels=lab)
  if(all(is.null(True))==FALSE) lines(True$D_t[seq(1,by=ns,length.out=Inputs$Data$n_y)], lwd=2)

}
    
if("Selex" %in% plot){

plot(x=1, y=1, type="n", xlim=c(min(mids),max(mids)), ylim=c(0, 1.1), ylab="Selectivity at length", xlab="Length (cm)", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2)
  if(all(is.null(LBSPR))==FALSE){
    for(i in 1:length(xLC_lbspr)){
      SL50 <- LBSPR$SL50[i]
      SL95 <- LBSPR$SL95[i]
      sd50 <- sqrt(LBSPR$var_S50[i])
      sd95 <- sqrt(LBSPR$var_S95[i])

      S_l2 <- 1.0/(1+exp(-log(19)*(mids-SL50)/(SL95-SL50))) # Selectivity-at-Length
      S_l2_low <- 1.0/(1+exp(-log(19)*(mids-(SL50-1.96*sd50))/((SL95-1.96*sd95)-(SL50-1.96*sd50)))) 
      S_l2_up <- 1.0/(1+exp(-log(19)*(mids-(SL50+1.96*sd50))/((SL95+1.96*sd95)-(SL50+1.96*sd50)))) 
      polygon(x=c(mids, rev(mids)), y=c(S_l2_low, rev(S_l2_up)), col="#AA00AA40", border=NA)
      lines(x=mids, y=S_l2, col="#AA00AA", lwd=2)
    }
  # legend("bottomright", col=c("#228B22", "#AA00AA", "black", "black","black"), lwd=2, legend=c("LIME", "LB-SPR", "SPR 40%", "SPR 30%", "Observed"), cex=1.7, lty=c(1,1,2,3,0), pch=c(19,19,NA,NA,17))
  }
  if(all(is.null(Sdreport))==FALSE){
    if(all(is.na(Sdreport))==FALSE){
      sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="S_fl"),]
      sd[,2][which(is.na(sd[,2]))] <- 0
      # ylim <- c(0, max(max(read_sdreport(sd, log=TRUE))*1.2))#, ymax))
    }
  }

  mids <- Inputs$Data$lbmids
  if(all(is.null(Report))==FALSE){
    if(nf>1){
      colfn <- colorRampPalette(c("red","blue"))
      cols <- colfn(nf)
    }
    if(nf==1) cols <- "#228B22"

    for(f in 1:nf){
      lines(x=mids, y=Report$S_fl[f,], lwd=2, col=cols[f])
      if(nf > 1) lty <- f+1
      if(nf ==1) lty <- 1
      if(all(is.null(True))==FALSE) lines(True$S_fl[f,], lwd=2, lty=lty)
      if(all(is.na(Sdreport))==FALSE){
        index <- seq(f,nrow(sd),by=nf)
        polygon( y=read_sdreport(sd[index,], log=FALSE), x=c(mids[which(is.na(sd[index,2])==FALSE)], rev(mids[which(is.na(sd[index,2])==FALSE)])), col=paste0(cols[f],"40"), border=NA)  
      }
    }
  }

  if(all(is.null(True))==FALSE) lines(True$S_fl[1,], lwd=2)
    # axis(1, cex.axis=2, at=xlabs, labels=plot_labs)
}

if(nf > 1 & legend==TRUE){
  legend("bottomright", col=c("#228B22", cols), legend=c("Total", sapply(1:nf, function(x) paste0("Fleet ", x))), lty=1, lwd=2)  
}

}