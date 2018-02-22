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
#' @param f_startval_ft default=NULL and F starting values are at 0 for all years. Can also specify vector of F starting values for all years to be modeled (can start at truth for debugging).
#' @param rdev_startval_t default=NULL and Recruitment deviation starting values are at 0 for all years. Can also specify vector of recruitment deviation starting values for all years to be modeled (can start at truth for debugging)
#' @param est_selex_f default=TRUE to estimate selectivity parameters, can set to FALSE for all or multiple fleets
#' @param vals_selex_ft input selectivity-at-length (columns) by fleet (rows) - negative values in the first column indicate to estimate selectivity
#' @param randomR default = TRUE, estimate recruitment as a random effect; if FALSE, turn off random effect on recruitment (do not derive deviations)
#' @param newtonsteps number of extra newton steps to take after optimization; FALSE to turn off
#' @param F_up upper bound of fishing mortality estimate; default=10
#' @param S50_up upper bound of length at 50 percent selectivity; default=NULL
#' @param derive_quants if TRUE, derive MSY-related reference points, default=FALSE
#' @param itervec number of datasets to generate in a simulation study. default=NULL for real stock assessment application. 
#' @param rewrite default=TRUE; if results already exist in the directory, should we rewrite them? TRUE or FALSE
#' @param simulation is this a simulation? default FALSE means you are using real data (can set itervec=NULL)
#' @param mirror vector of parameter names to mirror between fleets
#' @importFrom TMB MakeADFun sdreport
#' @importFrom TMBhelper Optimize
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
                      f_startval_ft=NULL,
                      rdev_startval_t=NULL,
                      est_selex_f=TRUE,
                      vals_selex_ft=-1,
                      randomR=TRUE,
                      newtonsteps=FALSE,
                      F_up=10,
                      S50_up=NULL,
                      derive_quants=FALSE,
                      itervec=NULL,
                      simulation=FALSE,
                      rewrite=TRUE,
                      mirror=NULL){

      # dyn.load(paste0(cpp_dir, "\\", dynlib("LIME")))

  if(simulation==FALSE) itervec <- 1 
  if(simulation==TRUE & is.null(itervec)) stop("Must specify number of iterations for simulation")    

for(iter in 1:length(itervec)){
    if(simulation==TRUE & is.null(modpath)==FALSE) iterpath <- file.path(modpath, iter)
    if(simulation==TRUE & is.null(modpath)) iterpath <- NULL
    if(simulation==FALSE & is.null(modpath)==FALSE) iterpath <- modpath
    if(simulation==FALSE & is.null(modpath)) iterpath <- NULL

    if(rewrite==FALSE & is.null(modpath)==FALSE){
      if(file.exists(file.path(iterpath, "Sdreport.rds"))) next
      if(file.exists(file.path(iterpath, "NAs_final_gradient.txt"))) next
      if(file.exists(file.path(iterpath, "high_final_gradient.txt"))) next
    }

    if(rewrite==TRUE & is.null(modpath)==FALSE){
      if(file.exists(file.path(iterpath, "NAs_final_gradient.txt"))) unlink(file.path(iterpath, "NAs_final_gradient.txt"), TRUE)
      if(file.exists(file.path(iterpath, "high_final_gradient.txt"))) unlink(file.path(iterpath, "high_final_gradient.txt"), TRUE)
      if(file.exists(file.path(iterpath, "model_NA.txt"))) unlink(file.path(iterpath, "model_NA.txt"), TRUE)
    }

    if(simulation==TRUE & is.null(modpath)==FALSE){
      stop("Need to edit code for multifleet")
      # sim <- readRDS(file.path(iterpath, "True.rds"))
      # f_inits <- sim$F_t
      # f_inits <- NULL
      # if(C_type==0) C_t_input <- NULL
      # if(C_type==1) C_t_input <- sim$Cn_t
      # if(C_type==2) C_t_input <- sim$Cw_t
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
    # if(inits$SigmaR > 0.05) SigRpen <- 0
    # if(inits$SigmaR <= 0.05) SigRpen <- 1
    if(is.null(modpath)) output <- NULL

    if(all(vals_selex_ft < 0)){
      vals_selex_ft_new <- matrix(-1, nrow=input$nfleets, ncol=length(input$highs))
    }
    if(any(vals_selex_ft >= 0)){
      vals_selex_ft_new <- vals_selex_ft
    }

      TmbList <- format_input(input=input, data_avail=data_avail, Fpen=Fpen, SigRpen=SigRpen, SigRprior=SigRprior, LFdist=LFdist, C_type=C_type, est_more=est_more, fix_more=fix_more, f_startval_ft=f_startval_ft, rdev_startval_t=rdev_startval_t, est_selex_f=est_selex_f, vals_selex_ft=vals_selex_ft_new, randomR=randomR, mirror=mirror)

      if(is.null(modpath)==FALSE) saveRDS(TmbList, file.path(iterpath, "Inputs.rds")) 
      if(is.null(modpath)) output$Inputs <- TmbList

      if(all(is.na(ParList))) ParList <- TmbList[["Parameters"]]  

      ## create objects to save best results
        obj_save <- NULL
        jnll <- NULL
        opt_save <- NULL
        opt_save[["final_gradient"]] <- NA

      ## first run
      obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList, random=TmbList[["Random"]], map=TmbList[["Map"]], inner.control=list(maxit=1e3), hessian=FALSE, DLL="LIME")

        # check_id <- Check_Identifiable(obj)
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
        Upr = rep(Inf, length(obj$par))
        Upr[match("log_sigma_R",names(obj$par))] = log(2)
        if(is.null(S50_up)==FALSE) Upr[which(names(obj$par)=="log_S50_f")] <- log(S50_up)
        if(is.null(S50_up)) Upr[which(names(obj$par)=="log_S50_f")] <- log(input$linf)
        Upr[which(names(obj$par)=="log_F_ft")] = log(F_up)
        Upr[match("log_sigma_F", names(obj$par))] <- log(2)

        Lwr <- rep(-Inf, length(obj$par))
        Lwr[match("log_CV_L",names(obj$par))] = log(0.001)
        Lwr[match("log_sigma_C",names(obj$par))] = log(0.001)
        Lwr[match("log_sigma_I",names(obj$par))] = log(0.001) 
        Lwr[which(names(obj$par)=="log_S50_f")] = log(1)
        Lwr[which(names(obj$par)=="log_F_ft")] <- log(0.001)

        ## Run optimizer
        # opt <- tryCatch( nlminb( start=obj$par, objective=obj$fn, gradient=obj$gr, upper=Upr, lower=Lwr, control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10) ), error=function(e) NA)    
        if(is.numeric(newtonsteps)) opt <- tryCatch(TMBhelper::Optimize(obj=obj, upper=Upr, lower=Lwr, newtonsteps=newtonsteps, getsd=FALSE), error=function(e) NA)
        if(is.numeric(newtonsteps)==FALSE) opt <- tryCatch(TMBhelper::Optimize(obj=obj, upper=Upr, lower=Lwr, loopnum=3, getsd=FALSE), error=function(e) NA)
        jnll <- obj$report()$jnll   
        if(all(is.na(opt))==FALSE & is.na(jnll)==FALSE){
          opt[["final_gradient"]] = obj$gr( opt$par ) 
          opt_save <- opt
          obj_save <- obj
          jnll_save <- obj_save$report()$jnll
          ParList <- list("log_F_ft"=log(obj_save$report()$F_ft), 
                          "log_q_f"=log(obj_save$report()$q_f), 
                          "beta"=obj_save$report()$beta,
                          "log_sigma_R"=log(obj_save$report()$sigma_R),
                          "log_S50_f"=log(obj_save$report()$S50),
                          "log_Sdelta_f"=log(obj_save$report()$S95 - obj_save$report()$S50),
                          "log_sigma_F"=log(obj_save$report()$sigma_F),
                          "log_sigma_C"=log(obj_save$report()$sigma_C),
                          "log_sigma_I"=log(obj_save$report()$sigma_I),
                          "log_CV_L"=log(obj_save$report()$CV_L),
                          "log_theta"=log(obj_save$report()$theta),
                          "Nu_input"=rep(0,length(TmbList$Parameters$Nu_input)))
        }      


        ## loop to try to get opt to run
          for(i in 1:5){
            if(all(is.na(opt)) | is.na(jnll) | all(is.na(opt_save[["final_gradient"]]))){
              obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList,
                            random=TmbList[["Random"]], map=TmbList[["Map"]], 
                            inner.control=list(maxit=1e3), hessian=FALSE, DLL="LIME")
                opt <-  tryCatch(TMBhelper::Optimize(obj=obj, start= obj$env$last.par.best[-obj$env$random] + rnorm(length(obj$par),0,0.2), upper=Upr, lower=Lwr, newtonsteps=newtonsteps, getsd=FALSE), error=function(e) NA)
                jnll <- obj$report()$jnll
            }
            if(all(is.na(opt))==FALSE & is.na(jnll)==FALSE){
              opt[["final_gradient"]] = obj$gr( opt$par )       
              opt_save <- opt
              obj_save <- obj
              jnll_save <- jnll
              ParList <- list("log_F_ft"=log(obj_save$report()$F_ft), 
                          "log_q_f"=log(obj_save$report()$q_f), 
                          "beta"=obj_save$report()$beta,
                          "log_sigma_R"=log(obj_save$report()$sigma_R),
                          "log_S50_f"=log(obj_save$report()$S50),
                          "log_Sdelta_f"=log(obj_save$report()$S95 - obj_save$report()$S50),
                          "log_sigma_F"=log(obj_save$report()$sigma_F),
                          "log_sigma_C"=log(obj_save$report()$sigma_C),
                          "log_sigma_I"=log(obj_save$report()$sigma_I),
                          "log_CV_L"=log(obj_save$report()$CV_L),
                          "log_theta"=log(obj_save$report()$theta),
                          "Nu_input"=rep(0,length(TmbList$Parameters$Nu_input)))
              break
            }
          }
          

        ## if opt ran: 
        if(all(is.na(opt_save))==FALSE){  

          ## check convergence -- don't let it become NA after it has had a high final gradient
          for(i in 1:5){
            if(all(is.na(opt_save[["final_gradient"]])==FALSE)){
               if(max(abs(opt_save[["final_gradient"]]))>0.001){
                obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList,
                            random=TmbList[["Random"]], map=TmbList[["Map"]], 
                            inner.control=list(maxit=1e3), hessian=FALSE, DLL="LIME")
                opt <-  tryCatch(TMBhelper::Optimize(obj=obj, start= ParList, upper=Upr, lower=Lwr, newtonsteps=newtonsteps, getsd=FALSE), error=function(e) NA)
                jnll <- obj$report()$jnll
              }
            }
              if(all(is.na(opt))==FALSE & is.na(jnll)==FALSE){
                if(is.null(jnll_save)){
                    opt[["final_gradient"]] = obj$gr( opt$par )       
                    opt_save <- opt
                    obj_save <- obj
                    jnll_save <- jnll
                    ParList <- list("log_F_ft"=log(obj_save$report()$F_ft), 
                          "log_q_f"=log(obj_save$report()$q_f), 
                          "beta"=obj_save$report()$beta,
                          "log_sigma_R"=log(obj_save$report()$sigma_R),
                          "log_S50_f"=log(obj_save$report()$S50),
                          "log_Sdelta_f"=log(obj_save$report()$S95 - obj_save$report()$S50),
                          "log_sigma_F"=log(obj_save$report()$sigma_F),
                          "log_sigma_C"=log(obj_save$report()$sigma_C),
                          "log_sigma_I"=log(obj_save$report()$sigma_I),
                          "log_CV_L"=log(obj_save$report()$CV_L),
                          "log_theta"=log(obj_save$report()$theta),
                          "Nu_input"=rep(0,length(TmbList$Parameters$Nu_input)))    
                }
                if(is.null(jnll_save)==FALSE){
                    if(jnll<=jnll_save){
                        opt[["final_gradient"]] = obj$gr( opt$par )       
                        opt_save <- opt
                        obj_save <- obj
                        jnll_save <- jnll
                    ParList <- list("log_F_ft"=log(obj_save$report()$F_ft), 
                          "log_q_f"=log(obj_save$report()$q_f), 
                          "beta"=obj_save$report()$beta,
                          "log_sigma_R"=log(obj_save$report()$sigma_R),
                          "log_S50_f"=log(obj_save$report()$S50),
                          "log_Sdelta_f"=log(obj_save$report()$S95 - obj_save$report()$S50),
                          "log_sigma_F"=log(obj_save$report()$sigma_F),
                          "log_sigma_C"=log(obj_save$report()$sigma_C),
                          "log_sigma_I"=log(obj_save$report()$sigma_I),
                          "log_CV_L"=log(obj_save$report()$CV_L),
                          "log_theta"=log(obj_save$report()$theta),
                          "Nu_input"=rep(0,length(TmbList$Parameters$Nu_input)))                     }
                }
              }
            if(all(is.na(opt_save[["final_gradient"]]))==FALSE){
              if(max(abs(opt_save[["final_gradient"]]))<=0.001) break
            }
          }
        }
        if(all(is.na(opt_save))==FALSE)  df <- data.frame(opt_save$final_gradient, names(obj_save$par), opt_save$par, exp(opt_save$par))


        ## write error message in directory if opt wouldn't run
          if(is.null(modpath)) output$issue <- NULL
          if(all(is.null(opt_save)==FALSE)){
            if(all(is.na(opt_save[["final_gradient"]]))==FALSE){
              if(max(abs(opt_save[["final_gradient"]]))>0.001){
                if(is.null(modpath)==FALSE) write(opt_save[["final_gradient"]], file.path(iterpath, "high_final_gradient.txt"))
                }
              if(is.null(modpath)) output$issue <- c(output$issue, "high_final_gradient")
              }
          }
          if(all(is.na(opt_save)) & is.null(modpath)==FALSE) write("model_NA", file.path(iterpath, "model_NA.txt"))
          if(all(is.na(opt_save)) & is.null(modpath)) output$issue <- c(output$issue, "model_NA")
        
        ## Standard errors
        Report = tryCatch( obj_save$report(), error=function(x) NA)
        if(is.null(modpath)==FALSE) saveRDS(Report, file.path(iterpath, "Report.rds"))  
        if(is.null(modpath)) output$Report <- Report

        Sdreport = tryCatch( sdreport(obj_save, bias.correct=TRUE), error=function(x) NA )
        if(is.null(modpath)==FALSE) saveRDS(Sdreport, file.path(iterpath, "Sdreport.rds"))
        if(is.null(modpath)) output$Sdreport <- Sdreport



# FUN <- function(InputMat, log=TRUE, rel=FALSE){
#           index <- which(is.na(InputMat[,2])==FALSE)
#           if(log==TRUE) return(c( exp(InputMat[index,1]-1.96*InputMat[index,2]), rev(exp(InputMat[index,1]+1.96*InputMat[index,2]))))
#           if(log==FALSE) return(c( InputMat[index,1]-1.96*InputMat[index,2], rev(InputMat[index,1]+1.96*InputMat[index,2])))
# } 
#   plot(Report$F_t, lwd=2, col="blue", ylim=c(0, max(Report$F_t)*1.5), type="l")
#   polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_t"),], log=TRUE), x=c(which(is.na(summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_t"),2])==FALSE), rev(which(is.na(summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)

      ## can calculate derived quants later if you want
      if(is.null(modpath)==FALSE) saveRDS(obj_save, file.path(iterpath, "TMB_obj.rds"))
      if(is.null(modpath)==FALSE) saveRDS(opt_save, file.path(iterpath, "TMB_opt.rds"))
      if(is.null(modpath)) output$obj <- obj_save
      if(is.null(modpath)) output$opt <- opt_save

      if(derive_quants==TRUE){
          Derived = calc_derived_quants( Obj=obj_save, lh=lh )
          if(is.null(modpath)==FALSE) saveRDS(Derived, file.path(iterpath, "Derived_quants.rds"))
          if(is.null(modpath)) output$Derived <- Derived
      }

        if(is.null(modpath)==FALSE) saveRDS(df, file.path(iterpath, "check_convergence.rds"))

        if(iter==1 & is.null(modpath)==FALSE) write.csv(df, file.path(modpath, "df.csv"))  
        if(iter==1 & is.null(modpath)) output$df <- df

        rm(Report)
        rm(Sdreport)
        rm(TmbList)
        rm(opt)
        rm(obj)
        rm(df)  
        rm(opt_save)
        rm(obj_save)
  }


if(is.null(modpath)==FALSE) return(paste0(max(itervec), " iterates run in ", modpath))
if(is.null(modpath)) return(output)

}
