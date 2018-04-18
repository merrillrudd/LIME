#' Predictive stacking
#'
#' \code{runstack} run predicted stacking for life history parameters using LIME
#'
#' @author M.B. Rudd
#' @param savedir directory to save results
#' @param lh life history list, with elements contained as output from LIME::create_lh_list
#' @param nodes matrix of nodes where each column is a different parameter and each row is a value from a distribution
#' @param param parameters (column names for nodes)
#' @param mean means of each parameter value,
#' @param cov covariance matrix across parameters
#' @param modname model name to save in directory
#' @param data_avail data available to model, default = "LC", can adjust to "Catch_LC", "Index_LC", or "Index_Catch_LC"
#' @param max_gradient maximum final gradient criterion, default = 0.001
#' @param C_type default = 0 (no catch data), 1=catch in numbers, 2= catch in biomass
#' @param LFdist default = 1 dirichlet-multinomial, 0=multinomial
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
#' 
runstack <- function(savedir, 
					lh, 
					nodes, 
					param, 
					mean, 
					cov, 
					modname, 
					data_avail="LC", 
					max_gradient=0.001, 
					C_type=0, 
					LFdist=1,
					rewrite=TRUE, 
					input_data=NULL, 
					simulation=FALSE, 
					iter=NULL, 
					seed=NULL, 
					Fscenario=NULL){

	## check inputs and find directories
	if(simulation == TRUE){
		if(is.null(iter)) stop("Must specify iteration number when simulation==TRUE")
		if(is.null(seed)) stop("Must specify seed when simulation==TRUE")
		if(is.null(Fscenario)) stop("Must specify Fscenario == harvestdyn or equil when simulation==TRUE")
		if(all(is.null(input_data))==FALSE) warning("Ignoring input_data for simulation")

		iterpath <- file.path(savedir, iter)
		dir.create(iterpath, showWarnings=FALSE)
	}
	if(simulation == FALSE){
		if(all(is.null(input_data))) stop("Must specify input data list with at least years and Length frequency matrix (array or list for multi-fleet)")
		if(is.null(iter)==FALSE | is.null(seed)==FALSE | is.null(Fscenario)==FALSE) warning("Ignoring iter, seed, and Fscenario when simulation==FALSE")
		iterpath <- file.path(savedir)
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
				SigmaF_inp <- 0.1
				SigmaR_inp <- 0.737
				rho_inp <- 0.4
				Fdynamics_inp <- "Endogenous"
			}
			plist <- with(lh, create_lh_list(linf=Linf_choose, vbk=vbk_choose, t0=t0,
									lwa=lwa, lwb=lwb,
									M=M_choose,
									M50=M50, maturity_input="age",
									S50=S50, S95=S95, selex_input="age",
									SigmaF=SigmaF_inp, SigmaR=SigmaR_inp, rho=rho_inp,
									AgeMax=AgeMax,
									binwidth=binwidth))

			# p <- ggplot(plist$df %>% filter(By=="Age")) +
			# 	geom_line(aes(x=X, y=Value, color=Fleet), lwd=2) + 
			# 	facet_grid(Variable~., scale="free_y") +
			# 	xlab("Age") +
			# 	guides(color=FALSE)
			# ggsave(file.path(iterpath, "LH_info.png"), p)

		if(is.list(input_data)==FALSE & is.data.frame(input_data)==FALSE){
			if(rewrite==TRUE | file.exists(file.path(iterpath, "True.rds"))==FALSE){
				## use seed + 1000 to generate data
				data <- generate_data(modpath=savedir, itervec=iter, 
								Fdynamics=Fdynamics_inp, Rdynamics="AR", 
								lh=plist, 
								Nyears=20, Nyears_comp=20, comp_sample=200,
								init_depl=c(0.10,0.90), 
								seed=rep(seed+1000,iter),
								rewrite=TRUE)
				# LFlist <- NULL
				# for(f in 1:plist$nfleets){
				# 	LFlist[[f]] <- data$LF[,,f]
				# }
				# png(file.path(iterpath, "LF_data.png"), height=8, width=10, res=200, units="in")
				# plot_LCfits(LFlist=LFlist, ylim=c(0,0.15))	
				# dev.off()
			}
			data <- readRDS(file.path(iterpath, "True.rds"))
			input_data <- list("years"=data$years, "LF"=data$LF)	
		}


		################
		## run models
		################
		## run at means from FishLife
		if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("res_FishLifeMeans.rds")))==FALSE){	

			## remove any flags or results files if re-running
			if(file.exists(file.path(iterpath,"nonconvergence_FishLifeMeans.txt"))) unlink(file.path(iterpath, "nonconvergence_FishLifeMeans.txt"), TRUE)
			if(file.exists(file.path(iterpath,"modelNA_FishLifeMeans.txt"))) unlink(file.path(iterpath, "modelNA_FishLifeMeans.txt"), TRUE)
			if(file.exists(file.path(iterpath,"res__FishLifeMeans.txt"))) unlink(file.path(iterpath, "res__FishLifeMeans.txt"), TRUE)	


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
				out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=3, C_type=C_type, LFdist=LFdist)

				## check_convergence
				isNA <- all(is.null(out$df))
				if(isNA==TRUE){
					gradient <- FALSE
					pdHess <- FALSE
				}
				if(isNA==FALSE){
					gradient <- out$opt$max_gradient <= max_gradient
					pdHess <- out$Sdreport$pdHess
				}	

				if(all(is.null(out$df))==FALSE & (gradient == FALSE | pdHess == FALSE)){
					out <- get_converged(results=out, saveFlagsDir=iterpath, saveFlagsName="IterTrue")
				}


				## flag non-convergence or NAs
				if(all(is.null(out$df))){
					write("model NA", file.path(iterpath, "modelNA_IterTrue.txt"))
				}
				if(all(is.null(out$df))==FALSE){
					gradient <- out$opt$max_gradient <= max_gradient
					pdHess <- out$Sdreport$pdHess
					if(gradient==FALSE){
						write("highgradient", file.path(iterpath,"highgradient_IterTrue.txt"))
					}
					if(pdHess==FALSE){
						write("Hessian not positive definite", file.path(iterpath, "pdHess_IterTrue.txt"))
					}
					## save results if converged
					if(gradient == TRUE & pdHess == TRUE) saveRDS(out, file.path(iterpath, paste0("res_IterTrue.rds")))	
				}
			}
	}

	## run at means from FishLife for ensemble parameters
	if(rewrite==TRUE | file.exists(file.path(iterpath, paste0(modname, "_res_FishLifeMeans.rds")))==FALSE){	

			## life history inputs
			vbk_inp <- ifelse("K" %in% param, exp(mean["K"]), lh$vbk)
			M_inp <- ifelse("M" %in% param, exp(mean["M"]), lh$M)
			linf_inp <- ifelse("Loo" %in% param, exp(mean["Loo"]), lh$linf)

			lhinp <- with(plist, 
					create_lh_list(linf=linf_inp, vbk=vbk_inp, t0=t0,
									lwa=lwa, lwb=lwb,
									M=M_inp,
									M50=M50, maturity_input="age",
									S50=S50, S95=S95, selex_input="age",
									SigmaF=SigmaF, SigmaR=SigmaR,
									AgeMax=AgeMax,
									binwidth=binwidth))		

			## input file and run model
			# input <- create_inputs(lh=lhinp, input_data=input_data)
			out <- run_LIME(modpath=NULL, lh=lhinp, input_data=input_data, est_sigma="log_sigma_R", data_avail="LC", rewrite=TRUE, newtonsteps=3)	

			## flag non-convergence or NAs
			if(all(is.null(out$df))) write("model NA", file.path(iterpath, "modelNA_FishLifeMeans.txt"))
			if(all(is.null(out$df))==FALSE){
				if(max(abs(out$df[,1]))>0.001) write("nonconvergence", file.path(iterpath,"nonconvergence_FishLifeMeans.txt"))
				## save results if converged
				if(max(abs(out$df[,1]))<=0.001) saveRDS(out, file.path(iterpath, paste0("res_FishLifeMeans.rds")))	
			}
			 		## if model doesn't converge:
			 		try <- 0
					while(try <= 3 & file.exists(file.path(iterpath, paste0("nonconvergence_FishLifeMeans.txt"))) | file.exists(file.path(iterpath, paste0("modelNA_FishLifeMeans.txt")))){		
						try <- try + 1
						## change starting values for F to those estimated in previous, nonconverged run
						if(all(is.null(out$df))==FALSE) lhinp$F_t <- out$Report$F_t
						if(all(is.null(out$df))) lhinp$F_t <- rnorm(length(input_data$years), mean=1, sd=0.2)
						out <- run_LIME(modpath=NULL, lh=lhinp, input=input_data, est_sigma="log_sigma_R", data_avail="LC", rewrite=TRUE, newtonsteps=3)	
						if(all(is.null(out$df))) write("modelNA", file.path(iterpath, paste0("modelNA_FishLifeMeans.txt")))

						if(max(abs(out$df[,1]))>0.001){
							write("nonconvergence", file.path(iterpath, paste0("nonconvergence_FishLifeMeans.txt")))
							if(file.exists(file.path(iterpath, paste0("modelNA_FishLifeMeans.txt")))) unlink(file.path(iterpath, paste0("modelNA_FishLifeMeans.txt")), TRUE)
						}

						if(all(is.null(out$df))==FALSE){
							if(max(abs(out$df[,1]))<=0.001){
								remove <- unlink(file.path(iterpath, paste0("nonconvergence_FishLifeMeans.txt")))
								saveRDS(out, file.path(iterpath, paste0("nonconvergence_FishLifeMeans.txt")))
							}
						}
						if(try == 3){
							if(all(is.null(out$df))==FALSE){
								if(max(abs(out$df[,1]))<=0.01){
									remove <- unlink(file.path(iterpath, paste0("nonconvergence_FishLifeMeans.txt")), TRUE)
									saveRDS(out, file.path(iterpath, paste0("res_FishLifeMeans.rds")))	
									write("convergence threshold 0.01", file.path(iterpath, paste0("minimal_convergence_FishLifeMeans.txt")))
								}
							}
							saveRDS(out, file.path(iterpath, paste0("res_FishLifeMeans_NC.rds")))
						}
					}
		}	

		## run at true values
		if(is.null(iter)==FALSE){
		if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("res_IterTrue.rds")))==FALSE){	
			

			## remove any flags or results files if re-running
			if(file.exists(file.path(iterpath,"nonconvergence_IterTrue.txt"))) unlink(file.path(iterpath, "nonconvergence_IterTrue.txt"), TRUE)
			if(file.exists(file.path(iterpath,"modelNA_IterTrue.txt"))) unlink(file.path(iterpath, "modelNA_IterTrue.txt"), TRUE)
			if(file.exists(file.path(iterpath,"res__IterTrue.txt"))) unlink(file.path(iterpath, "res__IterTrue.txt"), TRUE)	

			## input file and run model
			# input <- create_inputs(lh=plist, input_data=input_data)
			# input$SigmaF <- 0.1
			out <- run_LIME(modpath=NULL, lh=plist, input=input_data, est_sigma="log_sigma_R", data_avail="LC", rewrite=TRUE, newtonsteps=3)	

			## flag non-convergence or NAs
			if(max(abs(out$df[,1]))>0.001) write("nonconvergence", file.path(iterpath,"nonconvergence_IterTrue.txt"))
			if(all(is.null(out$df))) write("model NA", file.path(iterpath, "modelNA_IterTrue.txt"))

			## save results if converged
			if(max(abs(out$df[,1]))<=0.001) saveRDS(out, file.path(iterpath, paste0("res_IterTrue.rds")))	

			 	# 	## if model doesn't converge:
			 	# 	try <- 0
					# while(try <= 3 & file.exists(file.path(iterpath, paste0("nonconvergence_IterTrue.txt")))){		
					# 	try <- try + 1
					# 	## change starting values for F to those estimated in previous, nonconverged run
					# 	out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=3, f_startval_ft=out$Report$F_ft)	
					# 	if(max(abs(out$df[,1]))>0.001) write("nonconvergence", file.path(iterpath, paste0("nonconvergence_IterTrue.txt")))

					# 	if(all(is.null(out$df))) write("modelNA", file.path(iterpath, paste0("modelNA_IterTrue.txt")))
					# 	if(max(abs(out$df[,1]))<=0.001){
					# 		remove <- unlink(file.path(iterpath, paste0("nonconvergence_IterTrue.txt")))
					# 		saveRDS(out, file.path(iterpath, paste0("nonconvergence_IterTrue.txt")))
					# 	}
					# 	if(try == 3){
					# 		if(max(abs(out$df[,1]))<=0.01){
					# 			remove <- unlink(file.path(iterpath, paste0("nonconvergence_IterTrue.txt")), TRUE)
					# 			saveRDS(out, file.path(iterpath, paste0("res_IterTrue.rds")))	
					# 			write("convergence threshold 0.01", file.path(iterpath, paste0("minimal_convergence_IterTrue.txt")))
					# 		}
					# 	}
					# }
		}	
		}
	

		## predictive stacking
		if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("res_", modname, ".rds")))==FALSE){
			res <- lapply(1:nrow(nodes), function(x){

					## life history inputs -- nodes
					vbk_inp <- ifelse("K" %in% param, exp(nodes[x,"K"]), exp(mean["K"]))
					M_inp <- ifelse("M" %in% param, exp(nodes[x,"M"]), exp(mean["M"]))
					linf_inp <- ifelse("Loo" %in% param, exp(nodes[x,"Loo"]), exp(mean["Loo"]))
			 		lhinp <- with(plist, 
			 				 create_lh_list(linf=linf_inp, vbk=vbk_inp, t0=t0,
										lwa=lwa, lwb=lwb,
										M=M_inp,
										M50=M50, maturity_input="age",
										S50=S50, S95=S95, selex_input="age",
										SigmaF=0.1, SigmaR=SigmaR,
										AgeMax=AgeMax,
										binwidth=binwidth))		

			 		## input files and run model
					# input <- create_inputs(lh=lhinp, input_data=input_data)
					out <- run_LIME(modpath=NULL, lh=lhinp, input=input_data, est_sigma="log_sigma_R", data_avail="LC", rewrite=TRUE, newtonsteps=3)	

					## flag non-convergence or NAs
			 		if(max(abs(out$df[,1]))>0.001) write("nonconvergence", file.path(iterpath, paste0("nonconvergence_", modname, "node", x, ".txt")))
			 		if(all(is.null(out$df))) write("modelNA", file.path(iterpath, paste0("modelNA_", modname, "node", x, ".txt")))

			 		## save results if converged
			 		if(max(abs(out$df[,1]))<=0.001) saveRDS(out, file.path(iterpath, paste0("res_", modname, "node", x, ".txt")))

			 	# 	## if model doesn't converge:
			 	# 	try <- 0
					# while(try <= 3 & file.exists(file.path(iterpath, paste0("nonconvergence_", modname, "node", x, ".txt")))){		
					# 	try <- try + 1
					# 	## change starting values for F to those estimated in previous, nonconverged run
					# 	out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=3, f_startval_ft=out$Report$F_ft)	
					# 	if(max(abs(out$df[,1]))>0.001) write("nonconvergence", file.path(iterpath, paste0("nonconvergence_", modname, "node", x, ".txt")))

					# 	if(all(is.null(out$df))) write("modelNA", file.path(iterpath, paste0("modelNA_", modname, "node", x, ".txt")))
					# 	if(max(abs(out$df[,1]))<=0.001){
					# 		remove <- unlink(file.path(iterpath, paste0("nonconvergence_", modname, "node", x, ".txt")))
					# 		saveRDS(out, file.path(iterpath, paste0("res_", modname, "node", x, ".txt")))
					# 	}
					# 	if(try == 3){
					# 		if(max(abs(out$df[,1]))<=0.01){
					# 			remove <- unlink(file.path(iterpath, paste0("nonconvergence_", modname, "node", x, ".txt")), TRUE)
					# 			saveRDS(out, file.path(iterpath, paste0("res_", modname, "node", x, ".txt")))
					# 			write("convergence threshold 0.01", file.path(iterpath, paste0("minimal_convergence_", modname, "node", x, ".txt")))
					# 		}
					# 	}
					# }
				return(out)

			})
			saveRDS(res, file.path(iterpath, paste0("res_", modname, ".rds")))
			files <- list.files(path=file.path(iterpath))
			remove <- files[grepl(paste0("res_", modname, "node"), files)]
			ignore <- sapply(1:length(remove), function(x) unlink(file.path(iterpath, remove[x]), TRUE))
		}
	

	return(paste0("Ran iter ", iter, " in ", savedir))

}
