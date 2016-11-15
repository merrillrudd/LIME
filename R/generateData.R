#' Generate data for simulation testing
#'
#' \code{generateData} Generates data from the operating model for use in simulation testing
#'
#' @param modpath directory to save generated data
#' @param modname name of model (to identify differences between different model runs)
#' @param itervec number of iterations of data to generate
#' @param spatial does asymptotic length vary by spatially? TRUE or FALSE
#' @param Fdynamics Specify name of pattern of fishing mortality dynamics, "Constant", "Endogenous", "Ramped", "Increasing", or "None"
#' @param Rdynamics Specify name of pattern of recruitment dynamics, "Constant", "Pulsed", "Pulsed_up", or "BH"
#' @param LType 1 (default): pool site-specific length composition into 1 dataset; 0: keep length composition collected at different sites separate by site
#' @param plotML plot of generated mean length data, default=FALSE
#' @param plotLF_compare plot length comp pooled vs not pooled, default=FALSE
#' @param plotLF plot length composition at chosen LType, default=FALSE
#' @param selex "asymptotic" for asympotitc selectivity in length composition generation (currently only option programmed)
#' @param write write generated dataset? default=TRUE. FALSE helpful for plotting.
#' @param lh_list list of life history inputs
#' @param data_avail_list list of other model settings
#' @param param_adjust vector of names of parameters to adjust true value, default FALSE to include no parameter
#' @param val vector of values aligning with names in param_adjust to adjust value to, default FALSE to include no parameter
#' @param rewrite TRUE will re-run OM and observation model. FALSE will skip if it's already written in directory.

#' @return print how many iterations were written into the model directory
#' @export
generateData <- function(modpath, modname, itervec, spatial, Fdynamics, Rdynamics, LType=1, plotML=FALSE, plotLF_compare=FALSE, plotLF=FALSE, selex="asymptotic", write=TRUE, lh_list, data_avail_list, rewrite, param_adjust=NULL, val=NULL){

    lh_num <- ifelse(grepl("LH1", modpath), 1, ifelse(grepl("LH2", modpath), 2, ifelse(grepl("LH3", modpath), 3, ifelse(grepl("LH4", modpath), 4, ifelse(grepl("LH5", modpath), 5, ifelse(grepl("CRSNAP", modpath), "CRSNAP", ifelse(grepl("SIGSUT", modpath), "SIGSUT", ifelse(grepl("HAKE", modpath), "HAKE", stop("No match to life history number")))))))))
  lh_choose <- lh_list[[lh_num]]
  if(is.null(param_adjust)==FALSE){
      for(pp in 1:length(param_adjust)){
        lh_choose[[param_adjust[pp]]] <- val[pp]
      }    
  }
  
  Nyears_comp <- data_avail_list$Nyears_comp
  Nyears <- data_avail_list$Nyears

  for(iter in itervec){

    iterpath <- file.path(modpath, iter)
    if(write==TRUE) dir.create(iterpath, showWarnings=FALSE) 
    if(rewrite==FALSE){
      if(file.exists(file.path(iterpath, "True.rds"))) next
    }

    ## simulated data with no spatial structure in growth
    DataList <- with(c(lh_choose, data_avail_list), SimData_LB(Nyears=Nyears, AgeMax=AgeMax, 
      M=M, F1=F1, h=h, S_a=S_a, qcoef=qcoef, Frate=Frate, Fequil=Fequil, 
      SigmaF=SigmaF, Fdynamics=Fdynamics, Rdynamics=Rdynamics, 
      R0=R0, Fmax=Fmax, CVlen=CVlen, mids=mids, highs=highs, 
      lows=lows, W_a=W_a, L_a=L_a, Mat_a=Mat_a, M50=M50, 
      comp_sample=comp_sample, SigmaR=SigmaR, Nyears_comp=Nyears_comp, 
      alt_yrs=FALSE, sample=FALSE, nburn=20, seed=iter, modname=modname)) 

    # simulated data with spatial structure
    if(spatial==TRUE){

      set.seed(max(itervec)+iter)
      ## spatial process in Linf - varies with each iteration
      spatial_sim <- spatialgrowth_sim(n_i=20, linf=lh_choose$linf)
      if(write==TRUE) saveRDS(spatial_sim, file.path(iterpath, "spatial_sim.rds"))  

      ## life history - truth with spatial structure - varies with each iteration
      lh_spatial <- lapply(1:nrow(spatial_sim), function(x) choose_lh_list(species=lh_num, param_adjust=c("linf","ML50","binwidth","SigmaR","SigmaF"), val=c(spatial_sim[x,"linf_i"], lh_choose$ML50*(spatial_sim[x,"linf_i"]/lh_choose$linf), lh_choose$binwidth, lh_choose$SigmaR, lh_choose$SigmaF), selex="asymptotic", start_ages=ifelse(lh_choose$AgeMax < length(lh_choose$S_a), 0, 1))) 


      ## simulated data with spatial structure in growth
      DataList_site <- lapply(1:length(lh_spatial), function(x) with(c(lh_spatial[[x]], data_avail_list), SimData_LB(Nyears=Nyears, AgeMax=AgeMax,
          M=M, F1=F1, h=h, S_a=S_a, qcoef=qcoef, Frate=Frate, Fequil=Fequil, 
          SigmaF=SigmaF, Fdynamics=Fdynamics, Rdynamics=Rdynamics,
          R0=R0, Fmax=Fmax, CVlen=CVlen, mids=mids, highs=highs,
          lows=lows, W_a=W_a, L_a=L_a, Mat_a=Mat_a, M50=M50,
          comp_sample=comp_sample/nrow(spatial_sim), SigmaR=SigmaR, Nyears_comp=Nyears_comp,
          alt_yrs=FALSE, sample=FALSE, nburn=20, seed=iter, modname=modname)))  
      SPR_site <- sapply(1:length(DataList_site), function(x) DataList_site[[x]]$SPR_t)
      RelAbund_site <- sapply(1:length(DataList_site), function(x) DataList_site[[x]]$D_t[length(DataList_site[[x]]$D_t)])
      if(write==TRUE) saveRDS(SPR_site, file.path(iterpath, "SPR_site.rds"))

      ## length frequency data at each site
      if(Nyears_comp>1) LF_site <- lapply(1:length(DataList_site), function(x) DataList_site[[x]]$LF)
      if(Nyears_comp==1) LF_site <- lapply(1:length(DataList_site), function(x) as.matrix(DataList_site[[x]]$LF))
      ncols_site <- sapply(1:length(LF_site), function(x) ncol(LF_site[[x]]))
      for(i in 1:length(LF_site)){
        ncol <- ncol(LF_site[[i]])
        if(ncol < max(ncols_site)){
          add <- max(ncols_site) - ncol
          LF_site[[i]] <- cbind(LF_site[[i]], matrix(0, nrow=nrow(LF_site[[i]]), ncol=add))
        }
      } 

      ## length frequency by site
      LF_site_array <- array(NA, dim=c(dim(LF_site[[1]]), length(LF_site)))
      ML_t_site <- matrix(NA, nrow=nrow(LF_site[[1]]), ncol=length(LF_site))
      for(i in 1:length(LF_site)){
        LF_site_array[,,i] <- as.matrix(LF_site[[i]])
        if(Nyears_comp>1) ML_t_site[,i] <- sapply(1:nrow(LF_site_array[,,i]), function(x) sum(LF_site[[i]][x,]*1:ncol(LF_site[[i]]))/sum(LF_site[[i]][x,]))
        if(Nyears_comp==1) ML_t_site[,i] <- sum(LF_site[[i]]*1:length(LF_site[[i]]))/sum(LF_site[[i]])
      }
      rownames(LF_site_array) <- rownames(ML_t_site) <- (Nyears-Nyears_comp+1):Nyears #

      ## length frequency pooled across sites
      LF_pool <- matrix(NA, nrow=nrow(LF_site[[1]]), ncol=max(ncols_site))
      for(i in 1:nrow(LF_site[[1]])){
        for(j in 1:ncol(LF_site[[1]])){
          LF_pool[i,j] <- sum(sapply(1:length(LF_site), function(x) LF_site[[x]][i,j]))
        }
      }
      rownames(LF_pool) <- (Nyears-Nyears_comp+1):Nyears

    if(plotLF==TRUE){
      ## length frequency in the last year at each site
      par(mfrow=c(4,4), mar=c(0,0,0,0), omi=c(1,1,1,1))
      barplot(DataList$LF[nrow(DataList$LF),]/sum(DataList$LF[nrow(DataList$LF),]), axes=F, xlim=c(0,45), ylim=c(0,0.2), col="tomato")
      mtext(side=3, "no spatial process", font=2, line=-3, cex=2)
      axis(2, at=pretty(c(0,0.2)))      
      for(i in 1:length(LF_site)){
        barplot(LF_site[[i]][nrow(LF_site[[i]]),]/sum(LF_site[[i]][nrow(LF_site[[i]]),]), axes=F, xlim=c(0,45), ylim=c(0,0.2))
        mtext(paste0("site ", i), side=3, font=2, line=-3, cex=2)
        if(i %in% 12:15) axis(1, at=pretty(c(0,45)))
        if(i %% 4==0) axis(2, at=pretty(c(0,0.2)))
      }
      mtext(side=1, "Length bin (1 cm)", outer=TRUE, line=3)
      mtext(side=2, "Proportion", outer=TRUE, line=3)
    }
    if(plotLF_compare==TRUE){
      par(mfrow=c(2,1), mar=c(0,0,0,0), omi=c(1,1,1,1))
      barplot(DataList$LF[nrow(DataList$LF),]/sum(DataList$LF[nrow(DataList$LF),]), axes=F, xlim=c(0,45), ylim=c(0,0.2), col="tomato")
      mtext(side=3, "no spatial process", font=2, line=-3, cex=2)
      axis(2, at=pretty(c(0,0.2)))      

      barplot(LF_pool[nrow(LF_pool),]/sum(LF_pool[nrow(LF_pool),]), axes=F, xlim=c(0,45), ylim=c(0,0.2), col="black")
      mtext(side=3, "spatial process pooled", font=2, line=-3, cex=2)
      axis(1, at=pretty(c(0,45)))
      axis(2, at=pretty(c(0,0.2)))
      mtext(side=1, "Length bin (1 cm)", outer=TRUE, line=3)
      mtext(side=2, "Proportion", outer=TRUE, line=3)      
    }


      ## mean length over time pooled across sites
      ML_t_pool <- sapply(1:nrow(LF_pool), function(x) sum(LF_pool[x,]*1:ncol(LF_pool))/sum(LF_pool[x,])) 

      if(plotML==TRUE){
      # ## plot mean length at each site
      # # png("SIM_Mean_length_by_site_year.png", res=200, units="in", width=12, height=9)
      par(mfrow=c(4,4), mar=c(0,0,0,0), omi=c(1,1,1,1))
      plot(ML_t_pool, col="red", type="o", pch=17, lwd=2, ylim=c(0,linf), xaxt="n", yaxt="n")
      axis(2, cex=1.2, las=2)
      mtext(side=3, "pooled", font=2, line=-1.5)
      for(i in 1:length(DataList_site)){
        plot(DataList_site[[i]]$ML_t, col="black", pch=19, lwd=2, ylim=c(0,linf), xaxt="n", yaxt="n")
        lines(ML_t_pool, col="red", pch=17, type="o")
        if(i %in% c(12:15)) axis(1, cex=1.2)
        if(i %in% c(4,8,12)) axis(2, cex=1.2, las=2)
        mtext(paste0("site ", i), side=3, font=2, line=-1.5)
      }
      legend("bottomright", legend=c("site-specific", "pooled"), pch=c(19,17), col=c("black", "red"))
      mtext("Year", outer=TRUE, line=3, side=1)
      mtext("Mean length in catch (cm)", outer=TRUE, line=3,  side=2)
      # # dev.off()
    }
      if(LType==1){
        DataList$LF <- LF_pool
        DataList$ML_t <- ML_t_pool  
      }
      if(LType==0){
        DataList$LF <- LF_site_array
        DataList$ML_t <- ML_t_site
      }


      rm(spatial_sim)
      rm(lh_spatial)
    }
   
    DataList_out <- DataList
    if(nrow(DataList$LF)==1) DataList_out$LF <- t(as.matrix(DataList$LF[,1:length(lh_choose$mids)]))
    if(nrow(DataList$LF)>1) DataList_out$LF <- as.matrix(DataList$LF[,1:length(lh_choose$mids)])
    rownames(DataList_out$LF) <- rownames(DataList$LF)

    ## project the truth forward
    DataList_proj <- with(c(lh_choose, data_avail_list), SimData_LB(Nyears=Nyears, AgeMax=AgeMax, 
      M=M, F1=F1, h=h, S_a=S_a, qcoef=qcoef, Frate=Frate, Fequil=Fequil, 
      SigmaF=SigmaF, Fdynamics=Fdynamics, Rdynamics=Rdynamics, 
      R0=R0, Fmax=Fmax, CVlen=CVlen, mids=mids, highs=highs, 
      lows=lows, W_a=W_a, L_a=L_a, Mat_a=Mat_a, M50=M50, 
      comp_sample=comp_sample, SigmaR=SigmaR, Nyears_comp=Nyears, 
      alt_yrs=FALSE, sample=FALSE, nburn=20, seed=iter, modname=modname)) 
    Inputs <- FormatInput_LB(Nyears=Nyears, DataList=DataList_proj, linf=lh_choose$linf, vbk=lh_choose$vbk, t0=lh_choose$t0, M=lh_choose$M, AgeMax=lh_choose$AgeMax, lbhighs=lh_choose$highs, lbmids=lh_choose$mids, Mat_a=lh_choose$Mat_a, lwa=lh_choose$lwa, lwb=lh_choose$lwb, log_sigma_C=log(lh_choose$SigmaC), log_sigma_I=log(0.001), log_CV_L=log(0.001), F1=DataList$F_t[1], SigmaR=lh_choose$SigmaR, qcoef=lh_choose$qcoef, R0=mean(DataList_out$R_t), S50=lh_choose$S50, model="Rich_LC", Fpen=1, Dpen=0, Dprior=c(0,0), SigRpen=1, SigRprior=c(lh_choose$SigmaR, 0.2), obs_per_yr=rep(1000,Nyears), SigmaF=lh_choose$SigmaF, RecType=0, FType=0, LType=1, h=lh_choose$h, SelexTypeDesc="asymptotic", est_sigma="log_sigma_R", REML=FALSE, site=1, estimate_same=FALSE, start_f=0)
    ParList <- Inputs$Parameters
    obj <- MakeADFun(data=Inputs[["Data"]], parameters=ParList, random=Inputs[["Random"]], map=Inputs[["Map"]],inner.control=list(maxit=1e3), hessian=FALSE, DLL="LIME")  
    Derived <- Calc_derived_quants(obj)

    DataList_out$MSY <- Derived$MSY
    DataList_out$Fmsy <- Derived$Fmsy
    DataList_out$SBBmsy <- Derived$SBBmsy
    DataList_out$SBmsy <- Derived$SBmsy
    DataList_out$F30 <- Derived$F30
    DataList_out$FF30 <- Derived$FF30
    DataList_out$F40 <- Derived$F40
    DataList_out$FF40 <- Derived$FF40
    DataList_out$TBmsy <- Derived$TBmsy

      if(write==TRUE) saveRDS(DataList_out, file.path(iterpath, "True.rds"))
      if(write==FALSE) return(DataList_out)
      rm(DataList)
      rm(DataList_out)
      rm(iterpath)

}

  if(write==TRUE) return(paste0(length(itervec), " iterates of data generated in ", modpath))

}