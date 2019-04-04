#' Generate data for simulation testing
#'
#' \code{generate_data} Generates data from the operating model for use in simulation testing
#'
#' @author M.B. Rudd
#' @param modpath directory to save generated data
#' @param itervec number of iterations of data to generate
#' @param Fdynamics Specify name of pattern of fishing mortality dynamics, Constant, Endogenous, Ramp, Increasing, or None. Input number to project forward using a specific F.
#' @param Rdynamics Specify name of pattern of recruitment dynamics, Constant, Pulsed, Pulsed_up, or BH
#' @param lh list of life history information to feed to population simulation function, output from create_lh_list
#' @param pool if nseasons (in life history list) is greater than one, pool the generated data into annual time steps, or leave at the season-level? default=TRUE, FALSE will generate shorter time step life history info, mean length
#' @param Nyears number of years to simulate
#' @param Nyears_comp number of years to generate length composition data
#' @param comp_sample sample size of length composition data each year
#' @param rewrite TRUE will re-run OM and observation model. FALSE will skip if it's already written in directory.
#' @param init_depl default=0.4, can specify a different value or 2 values that indicate range from which to choose them
#' @param derive_quants default=FALSE (takes longer to run), can set to TRUE to output additional derived quantities.
#' @param seed single seed or vector of seeds for each iteration
#' @param mgt_type removals based on F (default) or catch
#' @param fleet_proportions vector specifying the relative size of each fleet in terms of fishing pressure. must have length = nfleets and sum to 1.
#' @param nareas number of areas, default = 1, if greater than 1, must be equal to the number of fleets
#' @importFrom stats runif
#' @importFrom TMB MakeADFun

#' @return print how many iterations were written into the model directory
#' @export
generate_data <- 
    function(modpath, 
            itervec, 
            Fdynamics, 
            Rdynamics, 
            lh, 
            pool=TRUE, 
            Nyears, 
            Nyears_comp, 
            comp_sample, 
            rewrite=TRUE, 
            init_depl, 
            derive_quants=FALSE, 
            seed, 
            mgt_type="F",
            fleet_proportions=1,
            nareas = 1){

    if(is.null(modpath) & length(itervec)>1) stop("must specify path to save simulation iterations")
    if(is.null(modpath)) itervec <- 1
    for(iter in itervec){

        iseed <- seed[iter]
        set.seed(iseed)

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
        DataList <- sim_pop(lh=lh, Nyears=Nyears, pool=pool, Fdynamics=Fdynamics, Rdynamics=Rdynamics, Nyears_comp=Nyears_comp, comp_sample=comp_sample, init_depl=init_depl_input, seed=iseed, mgt_type=mgt_type, fleet_proportions=fleet_proportions, nareas = nareas)
        if(all(is.na(DataList))==FALSE & all(is.null(modpath)==FALSE)) write(iseed, file.path(modpath, iter, paste0("init_depl_seed", iseed,".txt")))
         
    }

  
    ## if we are choosing randomly from a range of initial depletion:
    if(length(init_depl)==2){
        DataList <- NA
        add <- 0
        while(all(is.na(DataList))){
            seed_init <- iseed + add
            init_depl_input <- runif(1,init_depl[1],init_depl[2])
            ## simulated data with no spatial structure in growth
            DataList <- tryCatch(sim_pop(lh=lh, pool=pool, Nyears=Nyears, Fdynamics=Fdynamics, Rdynamics=Rdynamics, Nyears_comp=Nyears_comp, comp_sample=comp_sample, init_depl=init_depl_input, seed=seed_init, fleet_proportions=fleet_proportions, nareas = nareas), error=function(e) NA)
            if(all(is.na(DataList))==FALSE & all(is.null(modpath)==FALSE)) write(seed_init, file.path(modpath, iter, paste0("init_depl_seed", seed_init,".txt")))
            if(all(is.na(DataList))) add <- add + 1000
        }

    }
    if(length(init_depl)!=1 & length(init_depl)!=2) stop("init_depl must be a single proportion or 2 numbers inidicating minimum or maximum of range")

    inits <- create_inputs(lh=lh, input_data=DataList)

    # if(nrow(DataList$LF)==1) DataList_out$LF <- t(as.matrix(DataList$LF[,1:length(lh$mids)]))
    # if(nrow(DataList$LF)>1) DataList_out$LF <- as.matrix(DataList$LF[,1:length(lh$mids)])
    # rownames(DataList_out$LF) <- rownames(DataList$LF)

    if(derive_quants==TRUE){
        ## project the truth forward
        inits$C_ft <- DataList$Cw_ft
        TmbList <- format_input(input=inits, data_avail="Index_Catch_LC", Fpen=1, SigRpen=1, SigRprior=c(0.737,0.3), LFdist=1, C_type=2, est_more=FALSE, fix_more=FALSE, f_startval_ft=DataList$F_ft, rdev_startval_t=DataList$R_t, est_selex_f=FALSE, vals_selex_ft=DataList$S_fl, est_rdev_t=TRUE, mirror=FALSE, est_totalF=FALSE, prop_f=1, est_F_ft=TRUE)
        ParList <- TmbList$Parameters    

        # dyn.load(paste0(cpp_dir, "\\", dynlib("LIME")))       

        obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList, random=TmbList[["Random"]], map=TmbList[["Map"]],inner.control=list(maxit=1e3), hessian=FALSE, DLL="LIME")  
        Derived <- calc_derived_quants(Obj=obj, lh=lh) 

        inits$MSY <- Derived$MSY
        inits$Fmsy <- Derived$Fmsy
        inits$FFmsy <- Derived$FFmsy
        inits$SBBmsy <- Derived$SBBmsy
        inits$SBmsy <- Derived$SBmsy
        inits$F30 <- Derived$F30
        inits$FF30 <- Derived$FF30
        inits$F40 <- Derived$F40
        inits$FF40 <- Derived$FF40
        inits$TBmsy <- Derived$TBmsy
        inits$TBBmsy <- Derived$TBBmsy
    }

      if(is.null(modpath)==FALSE) saveRDS(inits, file.path(iterpath, "True.rds"))
      if(is.null(modpath)) return(inits)
      rm(DataList)
    rm(inits)
      rm(iterpath)

}


  if(is.null(modpath)==FALSE) return(modpath)

}