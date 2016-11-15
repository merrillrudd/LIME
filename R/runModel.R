#' run LIME model - simulation mode
#'
#' \code{runModel} run length-based integrated mixed-effects model with generated data
#'
#' @param modpath model directory
#' @param itervec number of iterations of data to generate
#' @param data_avail other settings for data availability
#' @param est_sigma which variance parameters to estimate, match parameter names
#' @param sensitivity_inputs (in development) named list (parameters) with matrix of 2 rows (low, high) and # of columns for 'life histories' (artifact of testing multiple life histories in the simulation, would only have 1 column for an assessment)
#' @param sensitivity_ESS (in development) will do sensitivity analysis for effective sample size of length composition
#' @param REML default off (FALSE)
#' @param estimate_same TRUE=estimate least-common-denominator parameters for all data availability scenarios, FALSE=estimate parameters specific to data availability
#' @param lh_list list of life history information
#' @param rewrite if results already exist in the directory, should we rewrite them? TRUE or FALSE
#' @param start_f year (in numbers, not actual year) to start estimating fishing mortality (e.g. year 11 out of 20 to get estimates for last 10 years); the value of F in this year will be used as the estimate and SE for all previous years. 0=estimate all years.
#' @param simulation is this a simulation? default TRUE, FALSE means you are using real data (no need for iterations or multiple life history inputs)
#' @param input_data use this to input data for a real-world application (not simulation)
#' @useDynLib LIME

#' @return prints how many iterations were run in model directory
#' 
#' @details need to adjust to run with real data
#' @export
runModel <- function(modpath, itervec, estimate_same=FALSE, REML=FALSE, est_sigma, data_avail, lh_list, sensitivity_inputs=NULL, sensitivity_ESS=NULL, rewrite, start_f, simulation=TRUE, input_data=NULL){

  if(simulation==TRUE){
    lh_num <- ifelse(grepl("LH1", modpath), 1, ifelse(grepl("LH2", modpath), 2, ifelse(grepl("LH3", modpath), 3, ifelse(grepl("LH4", modpath), 4, ifelse(grepl("LH5", modpath), 5, ifelse(grepl("CRSNAP", modpath), "CRSNAP", ifelse(grepl("SIGSUT", modpath), "SIGSUT", ifelse(grepl("HAKE", modpath), "HAKE", stop("No match to life history number")))))))))
    if(lh_num=="CRSNAP") lh_vec_num <- 1
    if(lh_num=="SIGSUT") lh_vec_num <- 2
    if(lh_num=="HAKE") lh_vec_num <- 3
    if(is.numeric(lh_num)) lh_vec_num <- lh_num
    lh_choose <- lh_list[[lh_num]]
  }
  if(simulation==FALSE){
    lh_choose <- lh_list
  }

  if(simulation==FALSE) itervec <- 1

for(iter in itervec){

    if(simulation==TRUE) iterpath <- file.path(modpath, iter)
    if(simulation==FALSE) iterpath <- modpath

    if(rewrite==FALSE){
      if(file.exists(file.path(iterpath, "Derived_quants.rds"))) next
      # if(any(grepl("LBSPR_results", list.files(path=iterpath)))) next
      if(file.exists(file.path(iterpath, "NAs_final_gradient.txt"))) next
      if(file.exists(file.path(iterpath, "high_final_gradient.txt"))) next
    }

    if(rewrite==TRUE){
      if(file.exists(file.path(iterpath, "NAs_final_gradient.txt"))) unlink(file.path(iterpath, "NAs_final_gradient.txt"), TRUE)
      if(file.exists(file.path(iterpath, "high_final_gradient.txt"))) unlink(file.path(iterpath, "high_final_gradient.txt"), TRUE)
    }

    if(simulation==TRUE){
        DataList <- readRDS(file.path(iterpath, "True.rds"))
        modname <- DataList$DataScenario 
    }

    if(simulation==TRUE) if(grepl("MixedEffects", modname)) modname <- strsplit(modname, "_")[[1]][2]
    if(simulation==FALSE) modname <- data_avail

    ## copies life history information with any adjustments for sensitivity analyses
    if(is.null(sensitivity_inputs)){
      param <- FALSE
      val <- FALSE
    }
    if(is.null(sensitivity_inputs)==FALSE){
      param_set <- names(sensitivity_inputs)
      param <- param_set[which(sapply(1:length(param_set), function(x) grepl(paste0("/",param_set[x],"/"), modpath)))]
      val_index <- ifelse(grepl("Low", modpath), 1, ifelse(grepl("High", modpath), 2, stop("Not set up for specified level of sensitivity")))
      if(simulation==TRUE) val <- as.numeric(sensitivity_inputs[[param]][val_index, lh_vec_num])
      if(simulation==FALSE) val <- as.numeric(sensitivity_inputs[[param]][val_index])
    }
    if(simulation==TRUE) inits <- create_inputs(lh_list=lh_choose, data_avail_list=c(data_avail[[modname]],DataList), param=param, val=val)
    if(simulation==FALSE) inits <- create_inputs(lh_list=lh_choose, data_avail_list=input_data, param=param, val=val) 
    Nyears <- inits$Nyears 

    if(simulation==TRUE) obs_per_yr <- inits$obs_per_yr
    if(simulation==FALSE) obs_per_yr <- input_data$obs_per_year

    if(simulation==FALSE){
        DataList <- inits
        saveRDS(DataList, file.path(iterpath, "obsData.rds"))
    }

  if(grepl("LBSPR", modpath)==FALSE){
    
    # if(biascorrect==FALSE) vec <- 1
    # if(biascorrect==TRUE) vec <- 1:2
    Sdreport <- NA
    ParList <- NA  
    df <- NULL

    # for(bb in vec){
      if(all(is.na(Sdreport))) RecDev_biasadj <- rep(0, Nyears)
      if(all(is.na(Sdreport))==FALSE){
          SD <- summary(Sdreport)
          RecDev_biasadj <- 1 - SD[which(rownames(SD)=="Nu_input"), "Std. Error"]^2 / Report$sigma_R^2    
      }
      if(all(is.na(RecDev_biasadj))) RecDev_biasadj <- rep(0, Nyears)
      TmbList <- FormatInput_LB(Nyears=Nyears, DataList=DataList, linf=inits$linf, vbk=inits$vbk, t0=inits$t0, M=inits$M, AgeMax=inits$AgeMax, lbhighs=inits$highs, lbmids=inits$mids, Mat_a=inits$Mat_a, lwa=inits$lwa, lwb=inits$lwb, log_sigma_C=inits$log_sigma_C, log_sigma_I=inits$log_sigma_I, log_CV_L=inits$log_CV_L, F1=inits$F1, SigmaR=inits$SigmaR, qcoef=inits$qcoef, R0=inits$R0, S50=inits$S50, model=as.character(modname), RecDev_biasadj=RecDev_biasadj,SigmaF=inits$SigmaF, Fpen=1, Dpen=0, Dprior=c(0.8, 0.2), SigRpen=1, SigRprior=c(inits$SigmaR, 0.2), obs_per_yr=obs_per_yr, RecType=0, FType=0, LType=1, h=inits$h, SelexTypeDesc="asymptotic", est_sigma=est_sigma, REML=REML, site=1, estimate_same=FALSE, start_f=start_f)
      saveRDS(TmbList, file.path(iterpath, "Inputs.rds")) 

      # if(bb==1) saveRDS(TmbList, file.path(iterpath, "Inputs1.rds")) 
      # if(bb==2) saveRDS(TmbList, file.path(iterpath, "Inputs2.rds"))  

      # dyn.load(paste0(run_exe, "\\", dynlib(version)))     

      if(all(is.na(ParList))) ParList <- TmbList[["Parameters"]]  

      ## create objects to save best results
      # if(bb==1){
        obj_save <- NULL
        jnll <- NULL
        opt_save <- NULL
        opt_save[["final_gradient"]] <- NA
      # }

      ## first run
      obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList, random=TmbList[["Random"]], map=TmbList[["Map"]],inner.control=list(maxit=1e3), hessian=FALSE, DLL="LIME")  
        # if(bb==1) check_id <- Check_Identifiable(obj)[[4]]
        # fix_f <- grep("Bad", check_id[which(check_id[,"Param"]=="log_F_t_input"),3])      
        # good_f <- c(1:Nyears)[which(1:Nyears %in% fix_f == FALSE)] 
        # TmbList$Map[["log_F_t_input"]] = 1:length(TmbList$Parameters[["log_F_t_input"]])
        # TmbList$Map[["log_F_t_input"]][fix_f] <- NA
        # TmbList$Map[["log_F_t_input"]] <- factor(TmbList$Map[["log_F_t_input"]])
        # if(length(fix_f)>0){
        #   TmbList$Data$fix_f <- fix_f
        #   TmbList$Data$fill_f <- good_f[length(good_f)]
        # }
      # if(bb==1) saveRDS(TmbList, file.path(iterpath, "Inputs1.1.rds"))
      # obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList, random=TmbList[["Random"]], map=TmbList[["Map"]],inner.control=list(maxit=1e3), hessian=FALSE, DLL="LIME")  


      ## Settings
      obj$env$inner.control <- c(obj$env$inner.control, "step.tol"=c(1e-8,1e-12)[1], "tol10"=c(1e-6,1e-8)[1], "grad.tol"=c(1e-8,1e-12)[1]) 
        Upr = rep(Inf, length(obj$par))
        Upr[match("log_sigma_R",names(obj$par))] = log(2)
        # Upr[match("logS95", names(obj$par))] = log(inits$AgeMax)
        Upr[match("logS50", names(obj$par))] = log(inits$AgeMax)
        Upr[which(names(obj$par)=="log_F_t_input")] = log(10)
        Upr[match("log_F_sd", names(obj$par))] <- log(2)
        Lwr <- rep(-Inf, length(obj$par))
        Lwr[match("logS50", names(obj$par))] = log(1)

        ## Run optimizer
        opt <- tryCatch( nlminb( start=obj$par, objective=obj$fn, gradient=obj$gr, upper=Upr, lower=Lwr, control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10) ), error=function(e) NA)    
        jnll <- obj$report()$jnll   
        if(all(is.na(opt))==FALSE & is.na(jnll)==FALSE){
          opt[["final_gradient"]] = obj$gr( opt$par ) 
          opt_save <- opt
          obj_save <- obj
          jnll_save <- obj_save$report()$jnll
        }      


        ## loop to try to get opt to run
          for(i in 1:5){
            if(all(is.na(opt))){
              obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList,
                            random=TmbList[["Random"]], map=TmbList[["Map"]], 
                            inner.control=list(maxit=1e3), hessian=FALSE, DLL="LIME")
                    obj$env$inner.control <- c(obj$env$inner.control, "step.tol"=c(1e-8,1e-12)[1], "tol10"=c(1e-6,1e-8)[1], "grad.tol"=c(1e-8,1e-12)[1]) 
                opt <- tryCatch( nlminb( start= obj$env$last.par.best[-obj$env$random] + rnorm(length(obj$par),0,0.2), 
                  objective=obj$fn, gradient=obj$gr, upper=Upr, lower=Lwr, 
                  control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10)), error=function(e) NA)
                jnll <- obj$report()$jnll
            }
            if(all(is.na(opt))==FALSE & is.na(jnll)==FALSE){
              opt[["final_gradient"]] = obj$gr( opt$par )       
              opt_save <- opt
              obj_save <- obj
              jnll_save <- jnll
              break
            }
          }
          

        ## if opt ran: 
        if(all(is.na(opt_save))==FALSE){  

          ## check convergence -- don't let it become NA after it has had a high final gradient
          for(i in 1:2){
            if(all(is.na(opt_save[["final_gradient"]])==FALSE)){
               if(abs(min(opt_save[["final_gradient"]]))>0.01){
                obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList,
                            random=TmbList[["Random"]], map=TmbList[["Map"]], 
                            inner.control=list(maxit=1e3), hessian=FALSE, DLL="LIME")
                      obj$env$inner.control <- c(obj$env$inner.control, "step.tol"=c(1e-8,1e-12)[1], "tol10"=c(1e-6,1e-8)[1], "grad.tol"=c(1e-8,1e-12)[1]) 
                opt <- tryCatch( nlminb( start= obj$env$last.par.best[-obj$env$random] + rnorm(length(obj$par),0,0.2), 
                  objective=obj$fn, gradient=obj$gr, upper=Upr, lower=Lwr, 
                  control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10)), error=function(e) NA)
                jnll <- obj$report()$jnll
              }
            }
              if(all(is.na(opt))==FALSE & is.na(jnll)==FALSE){
                if(is.null(jnll_save)){
                    opt[["final_gradient"]] = obj$gr( opt$par )       
                    opt_save <- opt
                    obj_save <- obj
                    jnll_save <- jnll
                }
                if(is.null(jnll_save)==FALSE){
                    if(jnll<=jnll_save){
                        opt[["final_gradient"]] = obj$gr( opt$par )       
                        opt_save <- opt
                        obj_save <- obj
                        jnll_save <- jnll
                    }
                }
              }
            if(abs(min(opt_save[["final_gradient"]]))<=0.01) break
          }
        }
        if(all(is.na(opt_save))==FALSE)  df <- data.frame(opt_save$final_gradient, names(obj_save$par), opt_save$par)


        ## write error message in directory if opt wouldn't run
        # if(bb==length(vec)){
          if(all(is.null(opt_save))) write("NAs final gradient", file.path(iterpath, "NAs_final_gradient.txt"))
          if(all(is.null(opt_save)==FALSE)) if(abs(min(opt_save[["final_gradient"]]))>0.01) write(opt_save[["final_gradient"]], file.path(iterpath, "high_final_gradient.txt"))
        # }

        ParList <- obj_save$env$parList( x=obj_save$par, par=obj_save$env$last.par.best )
        
        ## Standard errors
        Report = tryCatch( obj_save$report(), error=function(x) NA)
        # if(bb==length(vec)) saveRDS(Report, file.path(iterpath, "Report.rds"))  
        saveRDS(Report, file.path(iterpath, "Report.rds"))  

        Sdreport = tryCatch( sdreport(obj_save, bias.correct=TRUE), error=function(x) NA )
        # if(bb==length(vec)) saveRDS(Sdreport, file.path(iterpath, "Sdreport.rds"))
        saveRDS(Sdreport, file.path(iterpath, "Sdreport.rds"))


        # if(bb==length(vec)){
          Derived = Calc_derived_quants( Obj=obj_save )
          # if(bb==length(vec)) saveRDS(Derived, file.path(iterpath, "Derived_quants.rds"))
          saveRDS(Derived, file.path(iterpath, "Derived_quants.rds"))

        # }
      # } 
        if(iter==1) write.csv(df, file.path(modpath, "df.csv"))  

        rm(Report)
        rm(Sdreport)
        rm(TmbList)
        rm(opt)
        rm(obj)
        rm(df)  
        rm(opt_save)
        rm(obj_save)
  }

}

return(paste0(max(itervec), " iterates run in ", modpath))

}
