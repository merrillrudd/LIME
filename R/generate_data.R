#' Generate data for simulation testing
#'
#' \code{generate_data} Generates data from the operating model for use in simulation testing
#'
#' @param modpath directory to save generated data
#' @param data_avail types of data included, must at least include LCX where X is the number of years of length composition data. May also include Catch or Index separated by underscore. For example, LC10, Catch_LC1, Index_Catch_LC20.
#' @param itervec number of iterations of data to generate
#' @param Fdynamics Specify name of pattern of fishing mortality dynamics, Constant, Endogenous, Ramp, Increasing, or None
#' @param Rdynamics Specify name of pattern of recruitment dynamics, Constant, Pulsed, Pulsed_up, or BH
#' @param lh list of life history information to feed to population simulation function, output from create_lh_list
#' @param write write generated dataset? default=TRUE. FALSE helpful for plotting.
#' @param Nyears number of years to simulate
#' @param comp_sample vector with sample sizes of length composition data each year
#' @param rewrite TRUE will re-run OM and observation model. FALSE will skip if it's already written in directory.
#' @param init_depl default=0.4, can specify a different value or "random" to generate random initial depletions per iteration

#' @return print how many iterations were written into the model directory
#' @export
generate_data <- function(modpath, data_avail, itervec, Fdynamics, Rdynamics, lh, write=TRUE, Nyears, comp_sample, rewrite=TRUE, init_depl){

    for(iter in itervec){

    iterpath <- file.path(modpath, iter)
    if(write==TRUE) dir.create(iterpath, showWarnings=FALSE) 
    if(rewrite==FALSE){
      if(file.exists(file.path(iterpath, "True.rds"))) next
    }

    ## find out how many years of length comp data are available 
    if(grepl("LBSPR", data_avail)==FALSE){
        split_name <- unlist(strsplit(data_avail, "_"))
        lc_name <- split_name[grepl("LC", split_name)]
        Nyears_comp <- as.numeric(unlist(strsplit(lc_name, "LC"))[2])
    }
    if(grepl("LBSPR", data_avail)){
        Nyears_comp <- as.numeric(unlist(strsplit(data_avail, "LBSPR"))[2])
    }

    if(init_depl!="random") init_depl_input <- init_depl
    if(init_depl=="random"){
        set.seed(iter+1000)
        init_depl_input <- runif(1,0.05,1)
    }

    ## simulated data with no spatial structure in growth
    DataList <- sim_pop(lh=lh, Nyears=Nyears, Fdynamics=Fdynamics, Rdynamics=Rdynamics, Nyears_comp=Nyears_comp, comp_sample=comp_sample, init_depl=init_depl_input, nburn=50, seed=iter, modname=data_avail)

    # simulated data with spatial structure
    # if(spatial==TRUE){

    #   set.seed(max(itervec)+iter)
    #   ## spatial process in Linf - varies with each iteration
    #   spatial_sim <- spatialgrowth_sim(n_i=20, linf=lh$linf)
    #   if(write==TRUE) saveRDS(spatial_sim, file.path(iterpath, "spatial_sim.rds"))  

    #   ## life history - truth with spatial structure - varies with each iteration
    #   lh_spatial <- lapply(1:nrow(spatial_sim), function(x) choose_lh_list(species=lh_num, param_adjust=c("linf","ML50","binwidth","SigmaR","SigmaF"), val=c(spatial_sim[x,"linf_i"], lh$ML50*(spatial_sim[x,"linf_i"]/lh$linf), lh$binwidth, lh$SigmaR, lh$SigmaF), start_ages=ifelse(lh$AgeMax < length(lh$S_a), 0, 1))) 


    #   ## simulated data with spatial structure in growth
    #   DataList_site <- lapply(1:length(lh_spatial), function(x) with(c(lh_spatial[[x]], data_avail_list), SimData_LB(Nyears=Nyears, AgeMax=AgeMax,
    #       M=M, F1=F1, h=h, S_a=S_a, qcoef=qcoef, Frate=Frate, Fequil=Fequil, 
    #       SigmaF=SigmaF, Fdynamics=Fdynamics, Rdynamics=Rdynamics,
    #       R0=R0, Fmax=Fmax, CVlen=CVlen, mids=mids, highs=highs,
    #       lows=lows, W_a=W_a, L_a=L_a, Mat_a=Mat_a, M50=M50,
    #       comp_sample=comp_sample/nrow(spatial_sim), SigmaR=SigmaR, Nyears_comp=Nyears_comp,
    #       alt_yrs=FALSE, sample=FALSE, nburn=20, seed=iter, modname=modname)))  
    #   SPR_site <- sapply(1:length(DataList_site), function(x) DataList_site[[x]]$SPR_t)
    #   RelAbund_site <- sapply(1:length(DataList_site), function(x) DataList_site[[x]]$D_t[length(DataList_site[[x]]$D_t)])
    #   if(write==TRUE) saveRDS(SPR_site, file.path(iterpath, "SPR_site.rds"))

    #   ## length frequency data at each site
    #   if(Nyears_comp>1) LF_site <- lapply(1:length(DataList_site), function(x) DataList_site[[x]]$LF)
    #   if(Nyears_comp==1) LF_site <- lapply(1:length(DataList_site), function(x) as.matrix(DataList_site[[x]]$LF))
    #   ncols_site <- sapply(1:length(LF_site), function(x) ncol(LF_site[[x]]))
    #   for(i in 1:length(LF_site)){
    #     ncol <- ncol(LF_site[[i]])
    #     if(ncol < max(ncols_site)){
    #       add <- max(ncols_site) - ncol
    #       LF_site[[i]] <- cbind(LF_site[[i]], matrix(0, nrow=nrow(LF_site[[i]]), ncol=add))
    #     }
    #   } 

    #   ## length frequency by site
    #   LF_site_array <- array(NA, dim=c(dim(LF_site[[1]]), length(LF_site)))
    #   ML_t_site <- matrix(NA, nrow=nrow(LF_site[[1]]), ncol=length(LF_site))
    #   for(i in 1:length(LF_site)){
    #     LF_site_array[,,i] <- as.matrix(LF_site[[i]])
    #     if(Nyears_comp>1) ML_t_site[,i] <- sapply(1:nrow(LF_site_array[,,i]), function(x) sum(LF_site[[i]][x,]*1:ncol(LF_site[[i]]))/sum(LF_site[[i]][x,]))
    #     if(Nyears_comp==1) ML_t_site[,i] <- sum(LF_site[[i]]*1:length(LF_site[[i]]))/sum(LF_site[[i]])
    #   }
    #   rownames(LF_site_array) <- rownames(ML_t_site) <- (Nyears-Nyears_comp+1):Nyears #

    #   ## length frequency pooled across sites
    #   LF_pool <- matrix(NA, nrow=nrow(LF_site[[1]]), ncol=max(ncols_site))
    #   for(i in 1:nrow(LF_site[[1]])){
    #     for(j in 1:ncol(LF_site[[1]])){
    #       LF_pool[i,j] <- sum(sapply(1:length(LF_site), function(x) LF_site[[x]][i,j]))
    #     }
    #   }
    #   rownames(LF_pool) <- (Nyears-Nyears_comp+1):Nyears

    # if(plotLF==TRUE){
    #   ## length frequency in the last year at each site
    #   par(mfrow=c(4,4), mar=c(0,0,0,0), omi=c(1,1,1,1))
    #   barplot(DataList$LF[nrow(DataList$LF),]/sum(DataList$LF[nrow(DataList$LF),]), axes=F, xlim=c(0,45), ylim=c(0,0.2), col="tomato")
    #   mtext(side=3, "no spatial process", font=2, line=-3, cex=2)
    #   axis(2, at=pretty(c(0,0.2)))      
    #   for(i in 1:length(LF_site)){
    #     barplot(LF_site[[i]][nrow(LF_site[[i]]),]/sum(LF_site[[i]][nrow(LF_site[[i]]),]), axes=F, xlim=c(0,45), ylim=c(0,0.2))
    #     mtext(paste0("site ", i), side=3, font=2, line=-3, cex=2)
    #     if(i %in% 12:15) axis(1, at=pretty(c(0,45)))
    #     if(i %% 4==0) axis(2, at=pretty(c(0,0.2)))
    #   }
    #   mtext(side=1, "Length bin (1 cm)", outer=TRUE, line=3)
    #   mtext(side=2, "Proportion", outer=TRUE, line=3)
    # }
    # if(plotLF_compare==TRUE){
    #   par(mfrow=c(2,1), mar=c(0,0,0,0), omi=c(1,1,1,1))
    #   barplot(DataList$LF[nrow(DataList$LF),]/sum(DataList$LF[nrow(DataList$LF),]), axes=F, xlim=c(0,45), ylim=c(0,0.2), col="tomato")
    #   mtext(side=3, "no spatial process", font=2, line=-3, cex=2)
    #   axis(2, at=pretty(c(0,0.2)))      

    #   barplot(LF_pool[nrow(LF_pool),]/sum(LF_pool[nrow(LF_pool),]), axes=F, xlim=c(0,45), ylim=c(0,0.2), col="black")
    #   mtext(side=3, "spatial process pooled", font=2, line=-3, cex=2)
    #   axis(1, at=pretty(c(0,45)))
    #   axis(2, at=pretty(c(0,0.2)))
    #   mtext(side=1, "Length bin (1 cm)", outer=TRUE, line=3)
    #   mtext(side=2, "Proportion", outer=TRUE, line=3)      
    # }


    #   ## mean length over time pooled across sites
    #   ML_t_pool <- sapply(1:nrow(LF_pool), function(x) sum(LF_pool[x,]*1:ncol(LF_pool))/sum(LF_pool[x,])) 

    #   if(plotML==TRUE){
    #   # ## plot mean length at each site
    #   # # png("SIM_Mean_length_by_site_year.png", res=200, units="in", width=12, height=9)
    #   par(mfrow=c(4,4), mar=c(0,0,0,0), omi=c(1,1,1,1))
    #   plot(ML_t_pool, col="red", type="o", pch=17, lwd=2, ylim=c(0,linf), xaxt="n", yaxt="n")
    #   axis(2, cex=1.2, las=2)
    #   mtext(side=3, "pooled", font=2, line=-1.5)
    #   for(i in 1:length(DataList_site)){
    #     plot(DataList_site[[i]]$ML_t, col="black", pch=19, lwd=2, ylim=c(0,linf), xaxt="n", yaxt="n")
    #     lines(ML_t_pool, col="red", pch=17, type="o")
    #     if(i %in% c(12:15)) axis(1, cex=1.2)
    #     if(i %in% c(4,8,12)) axis(2, cex=1.2, las=2)
    #     mtext(paste0("site ", i), side=3, font=2, line=-1.5)
    #   }
    #   legend("bottomright", legend=c("site-specific", "pooled"), pch=c(19,17), col=c("black", "red"))
    #   mtext("Year", outer=TRUE, line=3, side=1)
    #   mtext("Mean length in catch (cm)", outer=TRUE, line=3,  side=2)
    #   # # dev.off()
    # }
    #   if(LType==1){
    #     DataList$LF <- LF_pool
    #     DataList$ML_t <- ML_t_pool  
    #   }
    #   if(LType==0){
    #     DataList$LF <- LF_site_array
    #     DataList$ML_t <- ML_t_site
    #   }


    #   rm(spatial_sim)
    #   rm(lh_spatial)
    # }
   
    DataList_out <- DataList
    if(nrow(DataList$LF)==1) DataList_out$LF <- t(as.matrix(DataList$LF[,1:length(lh$mids)]))
    if(nrow(DataList$LF)>1) DataList_out$LF <- as.matrix(DataList$LF[,1:length(lh$mids)])
    rownames(DataList_out$LF) <- rownames(DataList$LF)

    ## project the truth forward
    inits <- create_inputs(lh=lh, input_data=DataList, param=FALSE, val=FALSE)
    Inputs <- format_input(input=inits, data_avail="Index_Catch_LC", Fpen=1, SigRpen=1, SigRprior=c(inits$SigmaR, 0.2), est_sigma="log_sigma_R", REML=FALSE, fix_f=0, f_startval=DataList$F_t, Sel0=0)
    ParList <- Inputs$Parameters

    # dyn.load(paste0(cpp_dir, "\\", dynlib("LIME")))     

    obj <- MakeADFun(data=Inputs[["Data"]], parameters=ParList, random=Inputs[["Random"]], map=Inputs[["Map"]],inner.control=list(maxit=1e3), hessian=FALSE, DLL="LIME")  
    Derived <- calc_derived_quants(obj)

    DataList_out$MSY <- Derived$MSY
    DataList_out$Fmsy <- Derived$Fmsy
    DataList_out$FFmsy <- Derived$FFmsy
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


  if(write==TRUE) return(modpath)

}