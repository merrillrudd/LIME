#' run LIME model - simulation mode
#'
#' \code{run_LIME} run length-based integrated mixed-effects model with generated data
#'
#' @author M.B. Rudd
#' @param modpath model directory
#' @param input tagged list of LIME inputs. Output from create_inputs.
#' @param data_avail types of data included, must at least include LCX where X is the number of years of length composition data. May also include "Catch" or "Index" separated by underscore. For example, "LC10", "Catch_LC1", "Index_Catch_LC20".
#' @param Fpen penalty on fishing mortality 0= off, 1=on
#' @param SigRpen penalty on sigmaR, 0=off, 1=on
#' @param SigRprior vector with prior info for sigmaR penalty, first term is the mean and second term is the standard deviation
#' @param LFdist likelihood distribution for length composition data, default=0 for multinomial, alternate=1 for dirichlet-multinomial
#' @param C_type  default=0, NO catch data available. Copt=1 means the catch is in numbers, Copt2 means the catch is in weight. 
#' @param est_more list of variance parameters to estimate, must match parameter names: log_sigma_R, log_sigma_C, log_sigma_I, log_CV_L, log_sigma_F
#' @param fix_more default=FALSE - parameters are fixed depending on the data available. Can also list vector of parameter names to fix at their starting values (use param_adjust and val_adjust to set these adjustments)
#' @param est_F_ft default=TRUE, otherwise 0 for off and 1 for on in matrix that matches fleets in rows and years in columns
#' @param f_startval_ft default=NULL and F starting values are at 0 for all years. Can also specify vector of F starting values for all years to be modeled (can start at truth for debugging).
#' @param rdev_startval_t default=NULL and Recruitment deviation starting values are at 0 for all years. Can also specify vector of recruitment deviation starting values for all years to be modeled (can start at truth for debugging)
#' @param est_selex_f default=TRUE to estimate selectivity parameters, can set to FALSE for all or multiple fleets
#' @param vals_selex_ft input selectivity-at-length (columns) by fleet (rows) - negative values in the first column indicate to estimate selectivity
#' @param est_rdev_t default=TRUE to estimate recruitment deviations, or specify vector with 0 to turn off deviations in a specific year and 1 to keep them on
#' @param newtonsteps number of extra newton steps to take after optimization; FALSE to turn off
#' @param F_up upper bound of fishing mortality estimate; default=10
#' @param S50_up upper bound of length at 50 percent selectivity; default=NULL
#' @param derive_quants if TRUE, derive MSY-related reference points, default=FALSE
#' @param itervec number of datasets to generate in a simulation study. default=NULL for real stock assessment application. 
#' @param rewrite default=TRUE; if results already exist in the directory, should we rewrite them? TRUE or FALSE
#' @param simulation is this a simulation? default FALSE means you are using real data (can set itervec=NULL)
#' @param mirror vector of parameter names to mirror between fleets
#' @param est_totalF TRUE estimate total F instead of by fleet
#' @param prop_f proportion of catch from each fleet
#' @importFrom TMB MakeADFun sdreport
#' @importFrom TMBhelper fit_tmb
#' @importFrom utils write.csv
#' 

#' @useDynLib LIME

#' @return prints how many iterations were run in model directory
#' 
#' @export
run_LIME <- function(modpath, 
                      input,
                      data_avail, 
                      Fpen=1,
                      SigRpen=1,
                      SigRprior=c(0.737,0.3),
                      LFdist=1,
                      C_type=0,
                      est_more=FALSE,
                      fix_more=FALSE,
                      est_F_ft=TRUE,
                      f_startval_ft=NULL,
                      rdev_startval_t=NULL,
                      est_selex_f=TRUE,
                      vals_selex_ft=-1,
                      est_rdev_t=TRUE,
                      newtonsteps=3,
                      F_up=10,
                      S50_up=NULL,
                      derive_quants=FALSE,
                      itervec=NULL,
                      simulation=FALSE,
                      rewrite=TRUE,
                      mirror=NULL,
                      est_totalF=FALSE,
                      prop_f=1){

                      # Fpen=1
                      # SigRpen=1
                      # SigRprior=c(0.737,0.3)
                      # LFdist=1
                      # est_more=FALSE
                      # fix_more=FALSE
                      # est_F_ft=TRUE
                      # f_startval_ft=NULL
                      # rdev_startval_t=NULL
                      # est_selex_f=TRUE
                      # vals_selex_ft=-1
                      # Rdet=FALSE
                      # newtonsteps=3
                      # F_up=10
                      # S50_up=NULL
                      # derive_quants=FALSE
                      # itervec=NULL
                      # simulation=FALSE
                      # rewrite=TRUE
                      # mirror=NULL
                      # est_totalF=FALSE
                      # prop_f=1

  if(simulation==FALSE) itervec <- 1 
  if(simulation==TRUE & is.null(itervec)) stop("Must specify number of iterations for simulation")    

for(iter in 1:length(itervec)){
    if(simulation==TRUE & is.null(modpath)==FALSE) iterpath <- file.path(modpath, iter)
    if(simulation==TRUE & is.null(modpath)) iterpath <- NULL
    if(simulation==FALSE & is.null(modpath)==FALSE) iterpath <- modpath
    if(simulation==FALSE & is.null(modpath)) iterpath <- NULL

    if(rewrite==FALSE & is.null(modpath)==FALSE){
      if(file.exists(file.path(iterpath, "LIME_output.rds"))){
        output <- readRDS(file.path(iterpath, "LIME_output.rds"))
        next
      }
    }

    if(rewrite==TRUE & is.null(modpath)==FALSE){
      if(file.exists(file.path(iterpath, "LIME_output.rds"))) unlink(file.path(iterpath, "LIME_output.rds"), TRUE)
    }

    if(simulation==TRUE & is.null(modpath)==FALSE){
      stop("Need to edit code for multifleet")
      # sim <- readRDS(file.path(iterpath, "True.rds"))
      # f_inits <- sim$F_ft
      # # f_inits <- NULL
      # if(C_type==0) C_t_input <- NULL
      # if(C_type==1) C_t_input <- sim$Cn_ft
      # if(C_type==2) C_t_input <- sim$Cw_ft
      # if(LFdist==0) obs_input <- sim$obs_per_year
      # if(LFdist==1) obs_input <- rep(0, sim$Nyears)
      # true_nt <- sim$Nyears/sim$nseasons
      # s_all <- as.vector(sapply(1:true_nt, function(x) rep(x,sim$nseasons)))
      # years_i <- s_all[as.numeric(rownames(sim$LF))]
      # input_data <- list("years"=1:sim$Nyears, "LF"=sim$LF, "years_i"=years_i, "I_t"=sim$I_t, "C_t"=C_t_input, "F_t"=f_inits)
    }
      
    Sdreport <- NA
    ParList <- NA  
    df <- NULL
    if("log_sigma_R" %in% fix_more) SigRpen <- 0
    output <- NULL
    output$input <- input
    output$data_avail <- data_avail
    if(grepl("catch", tolower(data_avail)) & C_type == 0) stop("If including catch data,  must specify C_type as 1 for numbers or 2 for biomass.")

    if(all(vals_selex_ft < 0)){
      vals_selex_ft_new <- matrix(-1, nrow=input$nfleets, ncol=length(input$highs))
    }
    if(any(vals_selex_ft >= 0)){
      vals_selex_ft_new <- vals_selex_ft
    }
    if(all(prop_f==1)){
      prop_f_inp <- rep(1/input$nfleets, input$nfleets)
    }
    if(all(prop_f!=1)){
      checksum <- sum(prop_f) == 1
      if(checksum) prop_f_inp <- prop_f
      if(checksum==FALSE) stop("prop_f must sum to 1 and be equal to nfleets, or set prop_f=1 to be equal across fleets")
    }

      TmbList <- format_input(input=input, 
                              data_avail=data_avail, 
                              Fpen=Fpen, 
                              SigRpen=SigRpen, 
                              SigRprior=SigRprior, 
                              LFdist=LFdist, 
                              C_type=C_type, 
                              est_more=est_more, 
                              fix_more=fix_more, 
                              est_F_ft=est_F_ft,
                              f_startval_ft=f_startval_ft, 
                              rdev_startval_t=rdev_startval_t, 
                              est_selex_f=est_selex_f, 
                              vals_selex_ft=vals_selex_ft_new, 
                              est_rdev_t=est_rdev_t, 
                              mirror=mirror,
                              est_totalF=est_totalF,
                              prop_f=prop_f_inp)
                            
                             # input=input 
                             #  data_avail=data_avail 
                             #  Fpen=Fpen 
                             #  SigRpen=SigRpen 
                             #  SigRprior=SigRprior 
                             #  LFdist=LFdist 
                             #  C_type=C_type 
                             #  est_more=est_more 
                             #  fix_more=fix_more 
                             #  est_F_ft=est_F_ft
                             #  f_startval_ft=f_startval_ft 
                             #  rdev_startval_t=rdev_startval_t 
                             #  est_selex_f=est_selex_f 
                             #  vals_selex_ft=vals_selex_ft_new 
                             #  randomR=randomR 
                             #  mirror=mirror 
                             #  est_totalF=est_totalF 
                             #  prop_f=prop_f_inp 

      output$Inputs <- TmbList

      # if(all(is.na(ParList))) ParList <- TmbList[["Parameters"]]  

      ## create objects to save best results
        obj_save <- NULL
        jnll <- NULL
        opt_save <- NULL
        opt_save[["final_gradient"]] <- NA

      ## first run
      obj <- MakeADFun(data=TmbList[["Data"]], parameters=TmbList[["Parameters"]], random=TmbList[["Random"]], map=TmbList[["Map"]], inner.control=list(maxit=1e3), hessian=FALSE, DLL="LIME")

      ## Settings
        Upr = rep(Inf, length(obj$par))
        Upr[match("log_sigma_R",names(obj$par))] = log(2)
        if(is.null(S50_up)==FALSE) Upr[which(names(obj$par)=="log_S50_f")] <- log(S50_up)
        if(is.null(S50_up)) Upr[which(names(obj$par)=="log_S50_f")] <- log(input$linf)
        Upr[which(names(obj$par)=="log_F_ft")] = log(F_up)
        Upr[match("log_sigma_F", names(obj$par))] <- log(2)
        # Upr[which(names(obj$par)=="log_theta")] <- log(10)

        Lwr <- rep(-Inf, length(obj$par))
        Lwr[match("log_CV_L",names(obj$par))] = log(0.001)
        Lwr[match("log_sigma_C",names(obj$par))] = log(0.001)
        Lwr[match("log_sigma_I",names(obj$par))] = log(0.001) 
        Lwr[which(names(obj$par)=="log_S50_f")] = log(1)

        ## Run optimizer
        # opt <- tryCatch( nlminb( start=obj$par, objective=obj$fn, gradient=obj$gr, upper=Upr, lower=Lwr, control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10) ), error=function(e) NA)    
        # opt <- TMBhelper::fit_tmb(obj=obj, upper=Upr, lower=Lwr, newtonsteps=newtonsteps, getsd=FALSE)
        if(is.numeric(newtonsteps)) opt <- tryCatch(TMBhelper::fit_tmb(obj=obj, upper=Upr, lower=Lwr, newtonsteps=newtonsteps, getsd=FALSE), error=function(e) NA)
        if(is.numeric(newtonsteps)==FALSE) opt <- tryCatch(TMBhelper::fit_tmb(obj=obj, upper=Upr, lower=Lwr, loopnum=3, getsd=FALSE), error=function(e) NA)
        opt <- TMBhelper::fit_tmb(obj=obj, upper=Upr, lower=Lwr, loopnum=3, getsd=FALSE)

        # if(is.numeric(newtonsteps)) opt <- TMBhelper::fit_tmb(obj=obj, upper=Upr, lower=Lwr, newtonsteps=newtonsteps, getsd=FALSE)
        # if(is.numeric(newtonsteps)==FALSE) opt <- TMBhelper::fit_tmb(obj=obj, upper=Upr, lower=Lwr, loopnum=3, getsd=FALSE)
        jnll <- obj$report()$jnll   
        if(all(is.na(opt))==FALSE & is.na(jnll)==FALSE){
          opt[["final_gradient"]] = obj$gr( opt$par ) 
          opt_save <- opt
          obj_save <- obj
          jnll_save <- obj_save$report()$jnll
          # ParList <- obj$env$last.par.best
        }      


        ## loop to try to get opt to run
          for(i in 1:5){
            if(all(is.na(opt)) | is.na(jnll) | all(is.na(opt_save[["final_gradient"]]))){
              obj <- MakeADFun(data=TmbList[["Data"]], parameters=TmbList[["Parameters"]],
                            random=TmbList[["Random"]], map=TmbList[["Map"]], 
                            inner.control=list(maxit=1e3), hessian=FALSE, DLL="LIME")
                opt <-  tryCatch(TMBhelper::fit_tmb(obj=obj, start= obj$env$last.par.best[-obj$env$random] + rnorm(length(obj$par),0,0.2), upper=Upr, lower=Lwr, newtonsteps=newtonsteps, getsd=FALSE), error=function(e) NA)
                # opt <-  TMBhelper::fit_tmb(obj=obj, start= obj$env$last.par.best[-obj$env$random] + rnorm(length(obj$par),0,0.2), upper=Upr, lower=Lwr, newtonsteps=newtonsteps, getsd=FALSE)
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
          

        # ## if opt ran: 
        if(all(is.na(opt_save))==FALSE){  

          ## check convergence -- don't let it become NA after it has had a high final gradient
          for(i in 1:5){
            if(all(is.na(opt_save[["final_gradient"]])==FALSE)){
               if(max(abs(opt_save[["final_gradient"]]))>0.001){
                obj <- MakeADFun(data=TmbList[["Data"]], parameters=TmbList[["Parameters"]],
                            random=TmbList[["Random"]], map=TmbList[["Map"]], 
                            inner.control=list(maxit=1e3), hessian=FALSE, DLL="LIME")
                opt <-  tryCatch(TMBhelper::fit_tmb(obj=obj, start= obj$env$last.par.best[-obj$env$random] + rnorm(length(obj$par),0,0.2), upper=Upr, lower=Lwr, newtonsteps=newtonsteps, getsd=FALSE), error=function(e) NA)
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
            if(all(is.na(opt_save[["final_gradient"]]))==FALSE){
              if(max(abs(opt_save[["final_gradient"]]))<=0.001) break
            }
          }
        }
        if(all(is.na(opt_save))==FALSE)  df <- data.frame("gradient"=as.vector(opt_save$final_gradient), "parameter"=names(obj_save$par), "estimate"=opt_save$par, "transformed"=exp(opt_save$par))

        }


        ## Standard errors
        Report = tryCatch( obj_save$report(), error=function(x) NA)
        output$Report <- Report

        # Sdreport <- tryCatch( opt_save[["SD"]], error=function(x) NA)
        if(length(TmbList$Random) > 0) Sdreport = tryCatch(sdreport(obj_save, bias.correct=TRUE), error=function(x) NA )
        if(length(TmbList$Random) == 0) Sdreport <- tryCatch(sdreport(obj_save), error=function(x) NA)
        output$Sdreport <- Sdreport



      ## can calculate derived quants later if you want
      output$obj <- obj_save
      output$opt <- opt_save

      if(derive_quants==TRUE){
        if(all(is.na(Report))) Derived <- "Model NA"
        if(all(is.na(Report))==FALSE){
                    # MSY calculations
          Fmsy <- optimize(calc_msy, ages=input$ages, M=input$M, R0=exp(Report$beta), W_a=input$W_a, S_fa=Report$S_fa, lower=0, upper=10, maximum=TRUE)$maximum
          FFmsy <- Report$F_t/Fmsy
          msy <- calc_msy(F=Fmsy, ages=input$ages, M=input$M, R0=exp(Report$beta), W_a=input$W_a, S_fa=Report$S_fa)
          Bmsy <- sum(calc_equil_abund(ages=input$ages, M=input$M, F=Fmsy, S_fa=Report$S_fa, R0=exp(Report$beta)) * input$W_a)
          SBmsy <- sum(calc_equil_abund(ages=input$ages, M=input$M, F=Fmsy, S_fa=Report$S_fa, R0=exp(Report$beta)) * input$W_a * input$Mat_a)        

          SBBmsy <- Report$SB_t/SBmsy
          BBmsy <- Report$TB_t/Bmsy


          F30 <- tryCatch(uniroot(calc_ref, lower=0, upper=50, ages=input$ages, Mat_a=input$Mat_a, W_a=input$W_a, M=input$M, S_fa=Report$S_fa, ref=0.3)$root, error=function(e) NA) * input$nseasons
          F40 <- tryCatch(uniroot(calc_ref, lower=0, upper=50, ages=input$ages, Mat_a=input$Mat_a, W_a=input$W_a, M=input$M, S_fa=Report$S_fa, ref=0.4)$root, error=function(e) NA) * input$nseasons
          FF30 <- FF40 <- NULL
          if(is.na(F30)==FALSE) FF30 <- Report$F_t/F30
          if(is.na(F40)==FALSE) FF40 <- Report$F_t/F40

          Derived <- NULL
          Derived$Fmsy <- Fmsy
          Derived$FFmsy <- FFmsy
          Derived$msy <- msy
          Derived$Bmsy <- Bmsy
          Derived$BBmsy <- BBmsy
          Derived$SBmsy <- SBmsy
          Derived$SBBmsy <- SBBmsy
          Derived$F30 <- F30
          Derived$F40 <- F40
          Derived$FF30 <- FF30
          Derived$FF40 <- FF40
          output$Derived <- Derived
        }
      }

        output$df <- df

        rm(Report)
        rm(Sdreport)
        rm(TmbList)
        rm(opt)
        rm(obj)
        rm(df)  
        rm(opt_save)
        rm(obj_save)
  }

  if(is.null(modpath)==FALSE){
    if(rewrite==TRUE | file.exists(file.path(iterpath, "LIME_output.rds"))==FALSE){
      saveRDS(output, file.path(iterpath, "LIME_output.rds"))
    }
    if(rewrite==FALSE){
      output <- readRDS(file.path(iterpath, "LIME_output.rds"))
    }
  }

  return(output)


}
