#' Predictive stacking
#'
#' \code{runstack} run predicted stacking for life history parameters using LIME
#'
#' @author M.B. Rudd
#' @param savedir directory to save results
#' @param lh life history list, with elements contained as output from LIME::create_lh_list
#' @param nodes matrix of nodes where each column is a different parameter and each row is a value from a distribution
#' @param param parameters (column names for nodes)
#' @param mean mean of each parameter value 
#' @param cov covariance matrix across parameters
#' @param modname model name to save in directory
#' @param data_avail data available to model, default = "LC", can adjust to "Catch_LC", "Index_LC", or "Index_Catch_LC"
#' @param max_gradient maximum final gradient criterion, default = 0.001
#' @param rewrite rewrite results?
#' @param input_data for LIME use with real data (not simulated)
#' @param simulation is this a simulation? Default = FALSE
#' @param iter iteration of generated data
#' @param seed set seed
#' @param Fscenario fishing mortality scenario to generate data

#' @useDynLib LIME

#' @return prints how many iterations were run in model directory
#' 
#' @export
runstack <- function(savedir, lh, nodes, param, mean, cov, modname, data_avail="LC", max_gradient=0.001, rewrite=TRUE, input_data=NULL, simulation=FALSE, iter=NULL, seed=NULL, Fscenario=NULL){

	## check inputs and find directories
	if(simulation == TRUE){
		if(is.null(iter)) stop("Must specify iteration number when simulation==TRUE")
		if(is.null(seed)) stop("Must specify seed when simulation==TRUE")
		if(is.null(Fscenario)) stop("Must specify Fscenario == harvestdyn or equil when simulation==TRUE")
		if(all(is.null(input_data))==FALSE) warning("Ignoring input_data for simulation")

		modpath <- file.path(savedir, modname)
		dir.create(modpath, showWarnings=FALSE)

		iterpath <- file.path(modpath, iter)
		dir.create(iterpath, showWarnings=FALSE)
	}
	if(simulation == FALSE){
		if(all(is.null(input_data))) stop("Must specify input data list with at least years and Length frequency matrix (array or list for multi-fleet)")
		if(is.null(iter)==FALSE | is.null(seed)==FALSE | is.null(Fscenario)==FALSE) warning("Ignoring iter, seed, and Fscenario when simulation==FALSE")
		iterpath <- file.path(savedir, modname)
		dir.create(iterpath, showWarnings=FALSE)
	}


	## delete files if rewriting
	if(rewrite==TRUE){
		files <- list.files(iterpath)
		ignore <- sapply(1:length(files), function(x) unlink(file.path(iterpath, files[x]), TRUE))
	}

	## simulation only -- generate data and test run with true values
	if(simulation==TRUE){
		set.seed(seed)
		
		###################
		## generate data
		###################
			## true values OR based on FishLife OR will be ignored if not generating data OR starting values
			vbk_choose <- ifelse("K" %in% param, rlnorm(1, mean=mean["K"], sd=sqrt(cov["K","K"])), exp(mean["K"]))
			M_choose <- ifelse("M" %in% param, rlnorm(1, mean=mean["M"], sd=sqrt(cov["M","M"])), exp(mean["M"]))
			Linf_choose <- ifelse("Loo" %in% param, rlnorm(1, mean=mean["Loo"], sd=sqrt(cov["Loo","Loo"])), exp(mean["Loo"]))
			if(Fscenario=="equil"){
				SigmaF_inp <- 0.01
				SigmaR_inp <- 0.01
				rho_inp <- 0
				Fdynamics_inp <- "Constant"
			}
			if(Fscenario=="harvestdyn" | Fscenario==FALSE){
				SigmaF_inp <- 0.2
				SigmaR_inp <- 0.737
				rho_inp <- 0.4
				Fdynamics_inp <- "Endogenous"
			}
			plist <- with(lh, create_lh_list(linf=Linf_choose, vbk=vbk_choose, t0=t0,
									lwa=lwa, lwb=lwb,
									M=M_choose,
									M50=ML50, maturity_input="length",
									S50=SL50, S95=SL95, selex_input="length",
									SigmaF=SigmaF_inp, SigmaR=SigmaR_inp, rho=rho_inp,
									AgeMax=AgeMax,
									binwidth=binwidth,
									Fequil=1.1,
									theta=10))


			if(rewrite==TRUE | file.exists(file.path(iterpath, "True.rds"))==FALSE){
				## use seed + 1000 to generate data
					data <- generate_data(modpath=savedir, itervec=iter, 
									Fdynamics=Fdynamics_inp, Rdynamics="Constant", 
									lh=plist, 
									Nyears=20, Nyears_comp=20, comp_sample=200,
									init_depl=c(0.10,0.90), 
									seed=rep(seed+1000,iter),
									rewrite=TRUE)
			}
			data <- readRDS(file.path(iterpath, "True.rds"))
			input_data <- list("years"=data$years, "LF"=data$LF)

		## run at true values
		if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("res_IterTrue.rds")))==FALSE){	

				## input file and run model
				input <- create_inputs(lh=plist, input_data=input_data)
				out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=3)

				## check_convergence
				isNA <- all(is.null(out$df))
				if(isNA==FALSE){
					gradient <- out$opt$max_gradient <= max_gradient
					pdHess <- out$Sdreport$pdHess
				}

				if(isNA==FALSE & (gradient == FALSE | pdHess == FALSE)){
					out <- get_converged(results=out)
				}

				## flag non-convergence or NAs
				if(all(is.null(out$df))) write("model NA", file.path(iterpath, "modelNA_IterTrue.txt"))
				if(all(is.null(out$df))==FALSE){
					gradient <- out$opt$max_gradient <= max_gradient
					pdHess <- out$Sdreport$pdHess
					if(gradient==FALSE) write("highgradient", file.path(iterpath,"highgradient_IterTrue.txt"))
					if(pdHess==FALSE) write("Hessian not positive definite", file.path(iterpath, "pdHess_IterTrue.txt"))
					## save results if converged
					if(gradient == TRUE & pdHess == TRUE) saveRDS(out, file.path(iterpath, paste0("res_IterTrue.rds")))	
				}
			}
	}

	## run at means from FishLife for ensemble parameters
	if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("res_FishLifeMeans.rds")))==FALSE){	

			## life history inputs
			vbk_inp <- ifelse("K" %in% param, exp(mean["K"]), lh$vbk)
			M_inp <- ifelse("M" %in% param, exp(mean["M"]), lh$M)
			linf_inp <- ifelse("Loo" %in% param, exp(mean["Loo"]), lh$linf)
			lhinp <- with(lh, 
					create_lh_list(linf=linf_inp, vbk=vbk_inp, t0=t0,
									lwa=lwa, lwb=lwb,
									M=M_inp,
									M50=ML50, maturity_input="length",
									S50=SL50, S95=SL95, selex_input="length",
									SigmaF=SigmaF, SigmaR=SigmaR,
									AgeMax=AgeMax,
									binwidth=binwidth,
									theta=10))		

			## input file and run model
			input <- create_inputs(lh=lhinp, input_data=input_data)
			out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=3)	

				## check_convergence
				isNA <- all(is.null(out$df))
				if(isNA==FALSE){
					gradient <- out$opt$max_gradient <= max_gradient
					pdHess <- out$Sdreport$pdHess
				}

				if(isNA==FALSE & (gradient == FALSE | pdHess == FALSE)){
					out <- get_converged(results=out)
				}

				## flag non-convergence or NAs
				if(all(is.null(out$df))) write("model NA", file.path(iterpath, "modelNA_FishLifeMeans.txt"))
				if(all(is.null(out$df))==FALSE){
					gradient <- out$opt$max_gradient <= max_gradient
					pdHess <- out$Sdreport$pdHess
					if(gradient==FALSE) write("highgradient", file.path(iterpath,"highgradient_FishLifeMeans.txt"))
					if(pdHess==FALSE) write("Hessian not positive definite", file.path(iterpath, "pdHess_FishLifeMeans.txt"))
					## save results if converged
					if(gradient == TRUE & pdHess == TRUE) saveRDS(out, file.path(iterpath, paste0("res_FishLifeMeans.rds")))	
				}
	}	
	

	## predictive stacking
	if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("res_stacking.rds")))==FALSE){
		res <- lapply(1:nrow(nodes), function(x){

			## life history inputs -- nodes
			vbk_inp <- ifelse("K" %in% param, exp(nodes[x,"K"]), exp(mean["K"]))
			M_inp <- ifelse("M" %in% param, exp(nodes[x,"M"]), exp(mean["M"]))
			linf_inp <- ifelse("Loo" %in% param, exp(nodes[x,"Loo"]), exp(mean["Loo"]))
			lhinp <- with(lh, 
				 		create_lh_list(linf=linf_inp, vbk=vbk_inp, t0=t0,
										lwa=lwa, lwb=lwb,
										M=M_inp,
										M50=ML50, maturity_input="length",
										S50=SL50, S95=SL95, selex_input="length",
										SigmaF=SigmaF, SigmaR=SigmaR,
										AgeMax=AgeMax,
										binwidth=binwidth,
										theta=10))			

			# plot_LCfits(LFlist=list("LF"=input_data$LF[,,1]))

			## input files and run model
			input <- create_inputs(lh=lhinp, input_data=input_data)
			out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=3)		

					  # Fpen=1
       #                SigRpen=1
       #                SigRprior=c(0.737,0.3)
       #                LFdist=1
       #                C_type=0
       #                est_more=FALSE
       #                fix_more=FALSE
       #                f_startval_ft=NULL
       #                rdev_startval_t=NULL
       #                est_selex_f=TRUE
       #                vals_selex_ft=-1
       #                randomR=TRUE
       #                F_up=10
       #                S50_up=NULL
       #                derive_quants=FALSE
       #                itervec=NULL
       #                simulation=FALSE
       #                mirror=NULL
       #                est_totalF=FALSE
       #                prop_f=1

				## check_convergence
				isNA <- all(is.null(out$df))
				if(isNA==FALSE){
					gradient <- out$opt$max_gradient <= max_gradient
					pdHess <- out$Sdreport$pdHess
				}	

				if(isNA==FALSE & (gradient == FALSE | pdHess == FALSE)){
					out <- get_converged(results=out)
				}

				## flag non-convergence or NAs
				if(all(is.null(out$df))) write("model NA", file.path(iterpath, paste0("modelNA_node_", x, ".txt")))
				if(all(is.null(out$df))==FALSE){
					gradient <- out$opt$max_gradient<=0.001
					pdHess <- out$Sdreport$pdHess
					if(gradient==FALSE) write("highgradient", file.path(iterpath,paste0("highgradient_node_", x, ".txt")))
					if(pdHess==FALSE) write("Hessian not positive definite", file.path(iterpath, paste0("pdHess_node_", x, ".txt")))
					## save results if converged
					if(gradient == TRUE & pdHess == TRUE) saveRDS(out, file.path(iterpath, paste0("res_node_", x, ".rds")))	
				}
					
				return(out)
		})
			saveRDS(res, file.path(iterpath, paste0("res_stacking.rds")))
			files <- list.files(path=file.path(iterpath))
			remove <- files[grepl(paste0("res_node"), files)]
			ignore <- sapply(1:length(remove), function(x) unlink(file.path(iterpath, remove[x]), TRUE))
	}
	

	return(paste0("Ran iter ", iter, " in ", savedir))
}