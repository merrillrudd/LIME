#' Generate data for simulation testing
#'
#' \code{generate_data} Generates data from the operating model for use in simulation testing
#'
#' @author M.B. Rudd
#' @param modpath directory to save generated data
#' @param itervec number of iterations of data to generate
#' @param Fdynamics Specify name of pattern of fishing mortality dynamics, Constant, Endogenous, Ramp, Increasing, or None
#' @param Rdynamics Specify name of pattern of recruitment dynamics, Constant, Pulsed, Pulsed_up, or BH
#' @param lh list of life history information to feed to population simulation function, output from create_lh_list
#' @param pool if nseasons (in life history list) is greater than one, pool the generated data into annual time steps, or leave at the season-level? default=TRUE, FALSE will generate shorter time step life history info, mean length
#' @param Nyears number of years to simulate
#' @param Nyears_comp number of years to generate length composition data
#' @param comp_sample sample size of length composition data each year
#' @param rewrite TRUE will re-run OM and observation model. FALSE will skip if it's already written in directory.
#' @param mismatch default=FALSE, if TRUE, catch and index overlap with length comp only 1 year
#' @param init_depl default=0.4, can specify a different value or 2 values that indicate range from which to choose them
#' @param derive_quants default=FALSE (takes longer to run), can set to TRUE to output additional derived quantities.
#' @param nburn number of years burn-in for operating model, default=50
#' @importFrom stats runif
#' @importFrom TMB MakeADFun

#' @return print how many iterations were written into the model directory
#' @export
generate_data <- function(modpath, itervec, Fdynamics, Rdynamics, lh, pool=TRUE, Nyears, Nyears_comp, comp_sample, rewrite=TRUE, mismatch=FALSE, init_depl, derive_quants=FALSE, nburn=50){

    if(is.null(modpath) & length(itervec)>1) stop("must specify path to save simulation iterations")
    if(is.null(modpath)) itervec <- 1
    for(iter in itervec){

    if(is.null(modpath)==FALSE){
        iterpath <- file.path(modpath, iter)
        dir.create(iterpath, showWarnings=FALSE)
        if(rewrite==FALSE){
            if(file.exists(file.path(iterpath, "True.rds"))) next
        }
    }

    ## if level of depletion in first year is specified:
    if(length(init_depl)==1){
        init_depl_input <- init_depl
        ## simulated data with no spatial structure in growth
        DataList <- sim_pop(lh=lh, Nyears=Nyears, pool=pool, Fdynamics=Fdynamics, Rdynamics=Rdynamics, Nyears_comp=Nyears_comp, comp_sample=comp_sample, init_depl=init_depl_input, nburn=nburn, seed=iter, mismatch=mismatch)
    }
    ## if we are choosing randomly from a range of initial depletion:
    if(length(init_depl)==2){
        DataList <- NA
        add <- 0
        while(all(is.na(DataList))){
            seed_init <- iter + 1000 + add
            set.seed(seed_init)
            init_depl_input <- runif(1,init_depl[1],init_depl[2])
            ## simulated data with no spatial structure in growth
            DataList <- tryCatch(sim_pop(lh=lh, pool=pool, Nyears=Nyears, Fdynamics=Fdynamics, Rdynamics=Rdynamics, Nyears_comp=Nyears_comp, comp_sample=comp_sample, init_depl=init_depl_input, nburn=50, seed=iter, mismatch=mismatch), error=function(e) NA)
            if(all(is.na(DataList))==FALSE) write(seed_init, file.path(modpath, iter, paste0("init_depl_seed", seed_init,".txt")))
            if(all(is.na(DataList))) add <- add + 1000
        }

    }
    if(length(init_depl)!=1 & length(init_depl)!=2) stop("init_depl must be a single proportion or 2 numbers inidicating minimum or maximum of range")

   
    DataList_out <- DataList
    if(nrow(DataList$LF)==1) DataList_out$LF <- t(as.matrix(DataList$LF[,1:length(lh$mids)]))
    if(nrow(DataList$LF)>1) DataList_out$LF <- as.matrix(DataList$LF[,1:length(lh$mids)])
    rownames(DataList_out$LF) <- rownames(DataList$LF)

    if(derive_quants==TRUE){
        ## project the truth forward
        inits <- create_inputs(lh=lh, input_data=DataList)
        Inputs <- format_input(input=inits, data_avail="Index_Catch_LC", theta_type=0, Fpen=1, SigRpen=1, SigRprior=c(inits$SigmaR, 0.2), est_sigma="log_sigma_R", f_startval=DataList$F_t, LFdist=1, S_l_input=-1, randomR=TRUE)
        ParList <- Inputs$Parameters    

        # dyn.load(paste0(cpp_dir, "\\", dynlib("LIME")))       

        obj <- MakeADFun(data=Inputs[["Data"]], parameters=ParList, random=Inputs[["Random"]], map=Inputs[["Map"]],inner.control=list(maxit=1e3), hessian=FALSE, DLL="LIME")  
        Derived <- calc_derived_quants(Obj=obj, lh=lh) 

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
    }

      if(is.null(modpath)==FALSE) saveRDS(DataList_out, file.path(iterpath, "True.rds"))
      if(is.null(modpath)) return(DataList_out)
      rm(DataList)
      rm(DataList_out)
      rm(iterpath)

}


  if(is.null(modpath)==FALSE) return(modpath)

}