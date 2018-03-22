#' Predictive stacking
#'
#' \code{runstack} run predicted stacking for life history parameters using LIME
#'
#' @author M.B. Rudd
#' @param savedir directory to save results
#' @param iter iteration of generated data
#' @param seed set seed
#' @param lh life history list, with elements contained as output from LIME::create_lh_list
#' @param nodes matrix of nodes where each column is a different parameter and each row is a value from a distribution
#' @param param parameters (column names for nodes)
#' @param mean mean of each parameter value 
#' @param cov covariance matrix across parameters
#' @param modname model name to save in directory
#' @param input_data for LIME use with real data (not simulated)
#' @param Fscenario fishing mortality scenario to generate data
#' @param rewrite rewrite results?
#' @param binwidth default=1cm, can change

#' @useDynLib LIME

#' @return prints how many iterations were run in model directory
#' 
#' @export
runstack <- function(savedir, iter, seed, lh, nodes, param, mean, cov, modname, input_data, Fscenario, rewrite, binwidth=1){

	if(is.null(iter)==FALSE) iterpath <- file.path(savedir, iter)
	if(is.null(iter)) iterpath <- file.path(savedir, modname)
	dir.create(iterpath, showWarnings=FALSE)

		## delete files if rewriting
		if(rewrite==TRUE){
			files <- list.files(iterpath)
			ignore <- sapply(1:length(files), function(x) unlink(file.path(iterpath, files[x]), TRUE))
		}

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
									M50=M50, maturity_input="age",
									S50=S50, S95=S95, selex_input="age",
									SigmaF=SigmaF_inp, SigmaR=SigmaR_inp, rho=rho_inp,
									AgeMax=AgeMax,
									binwidth=binwidth,
									Fequil=1.1,
									theta=10))

		if(is.list(input_data)==FALSE & is.data.frame(input_data)==FALSE){
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
		}


		################
		## run models
		################
		## run at means from FishLife
		if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("res_FishLifeMeans.rds")))==FALSE){	

			## remove any flags or results files if re-running
			if(file.exists(file.path(iterpath,"pdHess_FishLifeMeans.txt"))) unlink(file.path(iterpath, "pdHess_FishLifeMeans.txt"), TRUE)
			if(file.exists(file.path(iterpath,"highgradient_FishLifeMeans.txt"))) unlink(file.path(iterpath, "highgradient_FishLifeMeans.txt"), TRUE)
			if(file.exists(file.path(iterpath,"modelNA_FishLifeMeans.txt"))) unlink(file.path(iterpath, "modelNA_FishLifeMeans.txt"), TRUE)
			if(file.exists(file.path(iterpath,"res__FishLifeMeans.txt"))) unlink(file.path(iterpath, "res__FishLifeMeans.txt"), TRUE)	

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
									binwidth=binwidth,
									theta=theta))		

			## input file and run model
			input <- create_inputs(lh=lhinp, input_data=input_data)
			out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=3)	

			## flag non-convergence or NAs
			if(all(is.null(out$df))) write("model NA", file.path(iterpath, "modelNA_FishLifeMeans.txt"))
			if(all(is.null(out$df))==FALSE){
				gradient <- out$opt$max_gradient<=0.001
				pdHess <- out$Sdreport$pdHess
				if(gradient==FALSE) write("highgradient", file.path(iterpath,"highgradient_FishLifeMeans.txt"))
				if(pdHess==FALSE) write("Hessian not positive definite", file.path(iterpath, "pdHess_FishLifeMeans.txt"))
				## save results if converged
				if(gradient == TRUE & pdHess == TRUE) saveRDS(out, file.path(iterpath, paste0("res_FishLifeMeans.rds")))	
			}

					## check and rerun in case of nonconvergence
					if(all(is.null(out$df))==FALSE & (gradient == FALSE | pdHess == FALSE)){
						## first check that theta is not estimated extremely high
						## often a problem that theta is estimated very large, and high final gradient is on selectivity
						## more important to estimate selectivity and fix theta at a high number
						if(out$Report$theta > 50){
							input$theta <- 50
							out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=3, fix_more="log_theta")
							
							## flag non-convergence or NAs
							if(all(is.null(out$df))) write("model NA", file.path(iterpath, paste0("modelNA_FishLifeMeans.txt")))
							if(all(is.null(out$df))==FALSE){
								gradient <- out$opt$max_gradient<=0.001
								pdHess <- out$Sdreport$pdHess
								if(gradient==FALSE) write("highgradient", file.path(iterpath,paste0("highgradient_FishLifeMeans.txt")))
								if(pdHess==FALSE) write("Hessian not positive definite", file.path(iterpath, paste0("pdHess_FishLifeMeans.txt")))
								## save results if converged
								if(gradient == TRUE & pdHess == TRUE){
									## remove flags if they exist
									unlink(file.path(iterpath,paste0( "pdHess_FishLifeMeans.txt")), TRUE)
									unlink(file.path(iterpath,paste0( "highgradient_FishLifeMeans.txt")), TRUE)
									saveRDS(out, file.path(iterpath, paste0("res_FishLifeMeans.rds")))
									write("fixed theta high", file.path(iterpath, "fixed_theta_high_FishLifeMeans.txt"))
								}	
							}						
						}

						if(gradient==FALSE){
							## fix parameter with high final gradient
							find_param <- as.character(out$df[,2][which(abs(out$df[,1])>=0.001)])
							out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=3, fix_more=find_param)

							## flag non-convergence or NAs
							if(all(is.null(out$df))) write("model NA", file.path(iterpath, paste0("modelNA_FishLifeMeans.txt")))
							if(all(is.null(out$df))==FALSE){
								gradient <- out$opt$max_gradient<=0.001
								pdHess <- out$Sdreport$pdHess
								if(gradient==FALSE) write("highgradient", file.path(iterpath,paste0("highgradient_FishLifeMeans.txt")))
								if(pdHess==FALSE) write("Hessian not positive definite", file.path(iterpath, paste0("pdHess_FishLifeMeans.txt")))
								## save results if converged
								if(gradient == TRUE & pdHess == TRUE){
									## remove flags if they exist
									unlink(file.path(iterpath,paste0( "pdHess_FishLifeMeans.txt")), TRUE)
									unlink(file.path(iterpath,paste0( "highgradient_FishLifeMeans.txt")), TRUE)
									saveRDS(out, file.path(iterpath, paste0("res_FishLifeMeans.rds")))
									ignore <- sapply(1:length(find_param), function(x) write("fixed high gradient parameter", file.path(iterpath, paste0("fixed_", find_param[x], "_FishLifeMeans", ".txt"))))
								}	
							}
						}
					}
		}	

		## run at true values
		if(is.null(iter)==FALSE){
			if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("res_IterTrue.rds")))==FALSE){	
		
				## remove any flags or results files if re-running
				if(file.exists(file.path(iterpath,"pdHess_IterTrue.txt"))) unlink(file.path(iterpath, "pdHess_IterTrue.txt"), TRUE)
				if(file.exists(file.path(iterpath,"highgradient_IterTrue.txt"))) unlink(file.path(iterpath, "highgradient_IterTrue.txt"), TRUE)
				if(file.exists(file.path(iterpath,"modelNA_IterTrue.txt"))) unlink(file.path(iterpath, "modelNA_IterTrue.txt"), TRUE)
				if(file.exists(file.path(iterpath,"res__IterTrue.txt"))) unlink(file.path(iterpath, "res__IterTrue.txt"), TRUE)		

				## input file and run model
				input <- create_inputs(lh=plist, input_data=input_data)
				# input$SigmaF <- 0.1
				out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=3)		

				## flag non-convergence or NAs
				if(all(is.null(out$df))) write("model NA", file.path(iterpath, "modelNA_IterTrue.txt"))
				if(all(is.null(out$df))==FALSE){
					gradient <- out$opt$max_gradient<=0.001
					pdHess <- out$Sdreport$pdHess
					if(gradient==FALSE) write("highgradient", file.path(iterpath,"highgradient_IterTrue.txt"))
					if(pdHess==FALSE) write("Hessian not positive definite", file.path(iterpath, "pdHess_IterTrue.txt"))
					## save results if converged
					if(gradient == TRUE & pdHess == TRUE) saveRDS(out, file.path(iterpath, paste0("res_IterTrue.rds")))	
				}

					## check and rerun in case of nonconvergence
					if(all(is.null(out$df))==FALSE & (gradient == FALSE | pdHess == FALSE)){
						## first check that theta is not estimated extremely high
						## often a problem that theta is estimated very large, and high final gradient is on selectivity
						## more important to estimate selectivity and fix theta at a high number
						if(out$Report$theta > 50){
							input$theta <- 50
							out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=3, fix_more="log_theta")
							
							## flag non-convergence or NAs
							if(all(is.null(out$df))) write("model NA", file.path(iterpath, paste0("modelNA_IterTrue.txt")))
							if(all(is.null(out$df))==FALSE){
								gradient <- out$opt$max_gradient<=0.001
								pdHess <- out$Sdreport$pdHess
								if(gradient==FALSE) write("highgradient", file.path(iterpath,paste0("highgradient_IterTrue.txt")))
								if(pdHess==FALSE) write("Hessian not positive definite", file.path(iterpath, paste0("pdHess_IterTrue.txt")))
								## save results if converged
								if(gradient == TRUE & pdHess == TRUE){
									## remove flags if they exist
									unlink(file.path(iterpath,paste0( "pdHess_IterTrue.txt")), TRUE)
									unlink(file.path(iterpath,paste0( "highgradient_IterTrue.txt")), TRUE)
									saveRDS(out, file.path(iterpath, paste0("res_IterTrue.rds")))
									write("fixed theta high", file.path(iterpath, "fixed_theta_high_IterTrue.txt"))
								}
							}						
						}

						if(gradient==FALSE){
							## fix parameter with high final gradient
							find_param <- as.character(out$df[,2][which(abs(out$df[,1])>=0.001)])
							out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=3, fix_more=find_param)

							## flag non-convergence or NAs
							if(all(is.null(out$df))) write("model NA", file.path(iterpath, paste0("modelNA_IterTrue.txt")))
							if(all(is.null(out$df))==FALSE){
								gradient <- out$opt$max_gradient<=0.001
								pdHess <- out$Sdreport$pdHess
								if(gradient==FALSE) write("highgradient", file.path(iterpath,paste0("highgradient_IterTrue.txt")))
								if(pdHess==FALSE) write("Hessian not positive definite", file.path(iterpath, paste0("pdHess_IterTrue.txt")))
								## save results if converged
								if(gradient == TRUE & pdHess == TRUE){
									## remove flags if they exist
									unlink(file.path(iterpath,paste0( "pdHess_IterTrue.txt")), TRUE)
									unlink(file.path(iterpath,paste0( "highgradient_IterTrue.txt")), TRUE)
									saveRDS(out, file.path(iterpath, paste0("res_IterTrue.rds")))
									ignore <- sapply(1:length(find_param), function(x) write("fixed high gradient parameter", file.path(iterpath, paste0("fixed_", find_param[x], "_IterTrue", ".txt"))))
								}	
							}
						}
					}
			}
		}
	

		## predictive stacking
		if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("res_", modname, ".rds")))==FALSE){
			res <- lapply(1:nrow(nodes), function(x){

			## remove any flags or results files if re-running
			if(file.exists(file.path(iterpath,paste0("pdHess_node_", x, ".txt")))) unlink(file.path(iterpath, paste0("pdHess_node_", x, ".txt")), TRUE)
			if(file.exists(file.path(iterpath,paste0("highgradient_node_", x, ".txt")))) unlink(file.path(iterpath,paste0( "highgradient_node_", x, ".txt")), TRUE)
			if(file.exists(file.path(iterpath,paste0("modelNA_node_", x, ".txt")))) unlink(file.path(iterpath,paste0( "modelNA_node_", x, ".txt")), TRUE)
			if(file.exists(file.path(iterpath,paste0("res__node_", x, ".txt")))) unlink(file.path(iterpath,paste0( "res__node_", x, ".txt")), TRUE)	

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
										SigmaF=SigmaF, SigmaR=SigmaR,
										AgeMax=AgeMax,
										binwidth=binwidth,
										theta=theta))		

			 		## input files and run model
					input <- create_inputs(lh=lhinp, input_data=input_data)
					out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=3)	

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

					## check and rerun in case of nonconvergence
					if(all(is.null(out$df))==FALSE & (gradient == FALSE | pdHess == FALSE)){
						## first check that theta is not estimated extremely high
						## often a problem that theta is estimated very large, and high final gradient is on selectivity
						## more important to estimate selectivity and fix theta at a high number
						if(out$Report$theta > 50){
							input$theta <- 50
							out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=3, fix_more="log_theta")
							
							## flag non-convergence or NAs
							if(all(is.null(out$df))) write("model NA", file.path(iterpath, paste0("modelNA_node_", x, ".txt")))
							if(all(is.null(out$df))==FALSE){
								gradient <- out$opt$max_gradient<=0.001
								pdHess <- out$Sdreport$pdHess
								if(gradient==FALSE) write("highgradient", file.path(iterpath,paste0("highgradient_node_", x, ".txt")))
								if(pdHess==FALSE) write("Hessian not positive definite", file.path(iterpath, paste0("pdHess_node_", x, ".txt")))
								## save results if converged
								if(gradient == TRUE & pdHess == TRUE){
									## remove flags if they exist
									unlink(file.path(iterpath,paste0( "pdHess_node_", x, ".txt")), TRUE)
									unlink(file.path(iterpath,paste0( "highgradient_node_", x, ".txt")), TRUE)
									saveRDS(out, file.path(iterpath, paste0("res_node_", x, ".rds")))
									write("fixed theta high", file.path(iterpath, paste0("fixed_theta_high_node_", x, ".txt")))
								}	
							}						
						}

						if(gradient==FALSE){
							## fix parameter with high final gradient
							find_param <- as.character(out$df[,2][which(abs(out$df[,1])>=0.001)])
							out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=3, fix_more=find_param)

							## flag non-convergence or NAs
							if(all(is.null(out$df))) write("model NA", file.path(iterpath, paste0("modelNA_node_", x, ".txt")))
							if(all(is.null(out$df))==FALSE){
								gradient <- out$opt$max_gradient<=0.001
								pdHess <- out$Sdreport$pdHess
								if(gradient==FALSE) write("highgradient", file.path(iterpath,paste0("highgradient_node_", x, ".txt")))
								if(pdHess==FALSE) write("Hessian not positive definite", file.path(iterpath, paste0("pdHess_node_", x, ".txt")))
								## save results if converged
								if(gradient == TRUE & pdHess == TRUE){
									## remove flags if they exist
									unlink(file.path(iterpath,paste0( "pdHess_node_", x, ".txt")), TRUE)
									unlink(file.path(iterpath,paste0( "highgradient_node_", x, ".txt")), TRUE)
									saveRDS(out, file.path(iterpath, paste0("res_node_", x, ".rds")))
									ignore <- sapply(1:length(find_param), function(x) write("fixed high gradient parameter", file.path(iterpath, paste0("fixed_", find_param[x], "_node_", x, ".txt"))))
								}	
							}
						}
					}
				return(out)
			})
			saveRDS(res, file.path(iterpath, paste0("res_", modname, ".rds")))
			files <- list.files(path=file.path(iterpath))
			remove <- files[grepl(paste0("res_node"), files)]
			ignore <- sapply(1:length(remove), function(x) unlink(file.path(iterpath, remove[x]), TRUE))
		}
	

	return(paste0("Ran iter ", iter, " in ", savedir))
}