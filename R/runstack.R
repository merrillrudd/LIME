#' Predictive stacking
#'
#' \code{runstack} run predicted stacking for life history parameters using LIME
#'
#' @author M.B. Rudd
#' @param savedir directory to save results
#' @param nodes matrix of nodes where each column is a different parameter and each row is a value from a distribution
#' @param param parameters (column names for nodes)
#' @param mean means of each parameter value,
#' @param cov covariance matrix across parameters
#' @param taxon taxonomic level for distributions
#' @param dim description of dimensional name
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
#' @param model default="LIME", alternate = "LBSPR"
#' @param sim_model default="LIME", alterate = "LBSPR"


#' @useDynLib LIME

#' @return prints how many iterations were run in model directory
#' 
#' @export
runstack <- function(savedir, 
					nodes, 
					param, 
					mean, 
					cov, 
					taxon,
					dim, 
					data_avail="LC", 
					max_gradient=0.001, 
					C_type=0, 
					LFdist=1,
					rewrite=TRUE, 
					input_data=NULL, 
					simulation=FALSE, 
					iter=NULL, 
					seed=NULL, 
					Fscenario=NULL,
					Nyears=NULL,
					model="LIME",
					sim_model="LIME"){

	## check inputs and find directories
	if(simulation == TRUE){
		if(is.null(iter)) stop("Must specify iteration number when simulation==TRUE")
		if(is.null(seed)) stop("Must specify seed when simulation==TRUE")
		if(is.null(Fscenario)) stop("Must specify Fscenario == harvestdyn or equil when simulation==TRUE")
		if(is.null(Nyears)) stop("Must specify number of years when simulation==TRUE")
		if(all(is.null(input_data))==FALSE) warning("Ignoring input_data for simulation")

		iterpath <- file.path(savedir, iter)
		dir.create(iterpath, showWarnings=FALSE)
	}
	if(simulation == FALSE){
		if(all(is.null(input_data))) stop("Must specify input data list with at least years and Length frequency matrix (array or list for multi-fleet)")
		if(is.null(iter)==FALSE | is.null(seed)==FALSE | is.null(Fscenario)==FALSE | is.null(Nyears)==FALSE) warning("Ignoring iter, seed, and Fscenario when simulation==FALSE")
		iterpath <- file.path(savedir)
		dir.create(iterpath, showWarnings=FALSE)
	}


	## delete files if rewriting
	if(rewrite==TRUE){
		files <- list.files(iterpath)
		ignore <- sapply(1:length(files), function(x) unlink(file.path(iterpath, files[x]), TRUE))
	}

	## simulation only -- generate data and test run with true values
	## should only be run at species level
	if(simulation==TRUE){
		set.seed(seed)
		
		###################
		## generate data
		###################
			lh_inp <- rmvnorm(1, mean=mean[c("Loo", "K","M","Lm")], sigma=cov[which(rownames(cov) %in% c("Loo", "K","M","Lm")), which(colnames(cov) %in% c("Loo", "K","M","Lm"))])

			if(Fscenario=="equil"){
				SigmaF_inp <- 0.001
				SigmaR_inp <- 0.001
				rho_inp <- 0
				# CVlen_inp <- 0.07
				Fdynamics_inp <- "Constant"
			}
			if(Fscenario=="harvestdyn" | Fscenario==FALSE){
				SigmaF_inp <- 0.1
				SigmaR_inp <- 0.737
				rho_inp <- 0.4
				# CVlen_inp <- lh$CVlen
				Fdynamics_inp <- "Endogenous"
			}
			plist <- create_lh_list(linf=exp(lh_inp[,"Loo"]), vbk=exp(lh_inp[,"K"]),
									lwa=0.01, lwb=3.04,
									M=exp(lh_inp[,"M"]),
									M50=exp(lh_inp[,"Lm"]), maturity_input="length",
									S50=exp(lh_inp[,"Lm"]), S95=min(exp(lh_inp[,"Loo"])*0.95, exp(lh_inp[,"Lm"])*1.2), selex_input="length",
									SigmaF=SigmaF_inp, SigmaR=SigmaR_inp, rho=rho_inp,
									# AgeMax=Amax_choose,
									binwidth=1,
									theta=10)

			if(rewrite==TRUE | file.exists(file.path(iterpath, "True.rds"))==FALSE){
				## use seed + 1000 to generate data
				if(sim_model=="LIME"){
					data <- generate_data(modpath=savedir, itervec=iter, 
									Fdynamics=Fdynamics_inp, Rdynamics="Constant", 
									lh=plist, 
									Nyears=Nyears, Nyears_comp=Nyears, comp_sample=200,
									init_depl=c(0.1,0.9), 
									seed=rep(seed+1000,iter),
									rewrite=TRUE)
				}
				if(sim_model=="LBSPR"){
					LB_pars <- new("LB_pars")
					LB_pars@MK <- plist$M/plist$vbk
					LB_pars@Linf <- plist$linf
					LB_pars@L50 <- plist$ML50
					LB_pars@L95 <- plist$ML95
					LB_pars@Walpha <- plist$lwa
					LB_pars@Wbeta <- plist$lwb
					LB_pars@BinWidth <- plist$binwidth	
					LB_pars@SL50 <- plist$SL50
					LB_pars@SL95 <- plist$SL95
					LB_pars@R0 <- plist$R0
					LB_pars@BinMin <- 0
					LB_pars@BinMax <- length(plist$mids)
					LB_pars@Steepness <- ifelse(plist$h==1, 0.99, plist$h)

            		init_depl_input <- runif(1,0.2,0.9)
            		LB_pars@SPR <- init_depl_input

            		sim <- LBSPRsim(LB_pars)
            		simlf <- rmultinom(Nyears, size=200, prob=sim@pLCatch)

            		data <- list()
            		data$LF <- t(simlf)
            		colnames(data$LF) <- sim@LMids
            		rownames(data$LF) <- 1:Nyears
            		data$mids <- sim@LMids
            		data$SPR <- sim@SPR
            		data$FM <- sim@FM
            		data$D_t <- sim@SSB/sim@SSB0
            		data$years <- 1:Nyears
            		data$SL50 <- sim@SL50
            		data$SL95 <- sim@SL95
            		data$S_fl <- matrix((1 /(1 + exp(-log(19)*(sim@LMids-sim@SL50)/(sim@SL95-sim@SL50)))), nrow=1)
            		saveRDS(data, file.path(iterpath, "True.rds"))
				}

			}
			data <- readRDS(file.path(iterpath, "True.rds"))
			input_data <- list("years"=data$years, "LF"=data$LF)

		## run at true values
		if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("res_IterTrue_", model, ".rds")))==FALSE){	

			input <- create_inputs(lh=plist, input_data=input_data)

			if(model=="LIME"){
				## input file and run model
				if(is.vector(input$LF[,,1])) Rdet <- TRUE
				if(is.vector(input$LF[,,1])==FALSE) Rdet <- FALSE
				out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=3, C_type=C_type, LFdist=LFdist, Rdet=Rdet)

				## check_convergence
				isNA <- all(is.null(out$df))
				if(isNA==TRUE){
					out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=FALSE, C_type=C_type, LFdist=LFdist, Rdet=Rdet)
				}
				isNA <- all(is.null(out$df))
				if(isNA==FALSE){
					gradient <- out$opt$max_gradient <= max_gradient
					pdHess <- out$Sdreport$pdHess
					if(pdHess==FALSE){
						input$theta <- 50
						out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=3, C_type=C_type, LFdist=LFdist, Rdet=Rdet, fix_more="log_theta")
					}
					isNA <- all(is.null(out$df))
					if(isNA==TRUE){
						out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=FALSE, C_type=C_type, LFdist=LFdist, Rdet=Rdet, fix_more="log_theta")
					}
					isNA <- all(is.null(out$df))
					if(isNA==FALSE){
						gradient <- out$opt$max_gradient <= max_gradient
						pdHess <- out$Sdreport$pdHess
					}	
				}	


				# if(all(is.null(out$df))==FALSE & (gradient == FALSE | pdHess == FALSE)){
				# 	out <- get_converged(results=out, saveFlagsDir=iterpath, saveFlagsName="IterTrue")
				# }


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
					if(gradient == TRUE & pdHess == TRUE) saveRDS(out, file.path(iterpath, paste0("res_IterTrue_LIME.rds")))	
				}
			}
			if(model=="LBSPR"){
				LB_lengths <- new("LB_lengths")
				LB_lengths@LMids <- input$mids
				LB_lengths@LData <- as.matrix(input$LF[,,1], ncol=Nyears)
				LB_lengths@Years <- as.numeric(rownames(input$LF))
				LB_lengths@NYears <- Nyears
				LB_lengths@L_units <- "cm"

					##----------------------------------------------------------------
					## Step 2: Specify biological inputs and parameter starting values
					##----------------------------------------------------------------
				LB_pars <- new("LB_pars")
				LB_pars@MK <- input$M/input$vbk
				LB_pars@Linf <- input$linf
				LB_pars@L50 <- input$ML50
				LB_pars@L95 <- input$ML95
				LB_pars@Walpha <- input$lwa
				LB_pars@Wbeta <- input$lwb
				LB_pars@BinWidth <- input$binwidth	
				LB_pars@SL50 <- input$SL50
				LB_pars@SL95 <- input$SL95
				LB_pars@R0 <- input$R0
				LB_pars@Steepness <- ifelse(input$h==1, 0.99, input$h)
				LB_pars@L_units <- "cm"
				LB_pars@BinMin <- 0
				LB_pars@BinMax <- input$linf * 1.3

				lbspr_res <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths)
				saveRDS(lbspr_res, file.path(iterpath, paste0("res_IterTrue_LBSPR.rds")))	

			}
		}
	}

	## run at means from FishLife for ensemble parameters
	if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("res_Means_", model,"_", taxon, ".rds")))==FALSE){	

			## life history inputs
			lhinp <- create_lh_list(linf=exp(mean["Loo"]), vbk=exp(mean["K"]),
									lwa=0.01, lwb=3.04,
									M=exp(mean["M"]),
									M50=exp(mean["Lm"]), maturity_input="length",
									S50=exp(mean["Lm"]), S95=min(exp(mean["Loo"])*0.95, exp(mean["Lm"])*1.2), selex_input="length",
									SigmaF=0.1, SigmaR=0.737,
									# AgeMax=exp(mean["tmax"]),
									binwidth=1,
									theta=10)	

		if(simulation==TRUE){
			data <- readRDS(file.path(iterpath, "True.rds"))
			input_data <- list("years"=data$years, "LF"=data$LF)
		}


			input <- create_inputs(lh=lhinp, input_data=input_data)

		if(model=="LIME"){
				## input file and run model
				if(is.vector(input$LF[,,1])) Rdet <- TRUE
				if(is.vector(input$LF[,,1])==FALSE) Rdet <- FALSE
				out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=FALSE, C_type=C_type, LFdist=LFdist, Rdet=Rdet)	

				## check_convergence
				isNA <- all(is.null(out$df))
				if(isNA==TRUE){
					out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=FALSE, C_type=C_type, LFdist=LFdist, Rdet=Rdet)
				}
				isNA <- all(is.null(out$df))
				if(isNA==FALSE){
					gradient <- out$opt$max_gradient <= max_gradient
					pdHess <- out$Sdreport$pdHess
					if(pdHess==FALSE){
						input$theta <- 50
						out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=3, C_type=C_type, LFdist=LFdist, Rdet=Rdet, fix_more="log_theta")
					}
					isNA <- all(is.null(out$df))
					if(isNA==TRUE){
						out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=FALSE, C_type=C_type, LFdist=LFdist, Rdet=Rdet, fix_more="log_theta")
					}
					isNA <- all(is.null(out$df))
					if(isNA==FALSE){
						gradient <- out$opt$max_gradient <= max_gradient
						pdHess <- out$Sdreport$pdHess
					}	
				}

				# if(all(is.null(out$df))==FALSE & (gradient == FALSE | pdHess == FALSE)){
				# 	out <- get_converged(results=out, saveFlagsDir=iterpath, saveFlagsName=paste0(modname, "_Means"))
				# }

				## flag non-convergence or NAs
				if(all(is.null(out$df))){
					write("model NA", file.path(iterpath, paste0("modelNA_Means_", taxon,".txt")))
				}
				if(all(is.null(out$df))==FALSE){
					gradient <- out$opt$max_gradient <= max_gradient
					pdHess <- out$Sdreport$pdHess
					if(gradient==FALSE){
						write("highgradient", file.path(iterpath, paste0("highgradient_Means_", taxon,".txt")))
					}
					if(pdHess==FALSE){
						write("Hessian not positive definite", file.path(iterpath, paste0("pdHess_Means_",taxon,".txt")))
					}
					## save results if converged
					if(gradient == TRUE & pdHess == TRUE) saveRDS(out, file.path(iterpath, paste0("res_Means_LIME_",taxon,".rds")))	
				}
		}
		if(model=="LBSPR"){
				LB_lengths <- new("LB_lengths")
				LB_lengths@LMids <- input$mids
				LB_lengths@LData <- as.matrix(input$LF[,,1], ncol=Nyears)
				LB_lengths@Years <- as.numeric(rownames(input$LF))
				LB_lengths@NYears <- Nyears
				LB_lengths@L_units <- "cm"

					##----------------------------------------------------------------
					## Step 2: Specify biological inputs and parameter starting values
					##----------------------------------------------------------------
				LB_pars <- new("LB_pars")
				LB_pars@MK <- input$M/input$vbk
				LB_pars@Linf <- input$linf
				LB_pars@L50 <- input$ML50
				LB_pars@L95 <- input$ML95
				LB_pars@Walpha <- input$lwa
				LB_pars@Wbeta <- input$lwb
				LB_pars@BinWidth <- input$binwidth	
				LB_pars@SL50 <- input$SL50
				LB_pars@SL95 <- input$SL95
				LB_pars@R0 <- input$R0
				LB_pars@Steepness <- ifelse(input$h==1, 0.99, input$h)
				LB_pars@BinMin <- 0
				LB_pars@BinMax <- input$linf * 1.3

				lbspr_res <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths)
				saveRDS(lbspr_res, file.path(iterpath, paste0("res_Means_LBSPR_",taxon,".rds")))	
		}
	}	
	

	## predictive stacking
	if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("res_stacking_", model, "_", taxon, "_", dim, ".rds")))==FALSE){
		res <- lapply(1:nrow(nodes), function(x){
		# for(x in 1:nrow(nodes)){
			## life history inputs -- nodes
			vbk_inp <- ifelse(any(param=="K"), exp(nodes[x,"K"]), exp(mean["K"]))
			M_inp <- ifelse(any(param=="M"), exp(nodes[x,"M"]), exp(mean["M"]))
			linf_inp <- ifelse(any(param=="Loo"), exp(nodes[x,"Loo"]), exp(mean["Loo"]))
			Lmat_inp <- ifelse(any(param=="Lm"), exp(nodes[x,"Lm"]), exp(mean["Lm"]))	
			Amax_inp <- ifelse(any(param=="tmax"), exp(nodes[x,"tmax"]), exp(mean["tmax"]))
			lhinp <- create_lh_list(linf=linf_inp, vbk=vbk_inp,
									lwa=0.01, lwb=3.04,
										M=M_inp,
										M50=Lmat_inp, maturity_input="length",
										S50=Lmat_inp, S95=min(linf_inp*0.95, Lmat_inp*1.2), selex_input="length",
										SigmaF=0.1, SigmaR=0.737,
										# AgeMax=Amax_inp,
										binwidth=1,
										theta=10)			

		if(simulation==TRUE){
			data <- readRDS(file.path(iterpath, "True.rds"))
			input_data <- list("years"=data$years, "LF"=data$LF)
		}


			input <- create_inputs(lh=lhinp, input_data=input_data)

		if(model=="LIME"){
				if(is.vector(input$LF[,,1])) Rdet <- TRUE
				if(is.vector(input$LF[,,1])==FALSE) Rdet <- FALSE
				out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=FALSE, C_type=C_type, LFdist=LFdist, Rdet=Rdet)	

				## check_convergence
				isNA <- all(is.null(out$df))
				if(isNA==TRUE){
					out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=FALSE, C_type=C_type, LFdist=LFdist, Rdet=Rdet)
				}
				isNA <- all(is.null(out$df))
				if(isNA==FALSE){
					gradient <- out$opt$max_gradient <= max_gradient
					pdHess <- out$Sdreport$pdHess
					if(pdHess==FALSE){
						input$theta <- 50
						out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=3, C_type=C_type, LFdist=LFdist, Rdet=Rdet, fix_more="log_theta")
					}
					isNA <- all(is.null(out$df))
					if(isNA==TRUE){
						out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=FALSE, C_type=C_type, LFdist=LFdist, Rdet=Rdet, fix_more="log_theta")
					}
					isNA <- all(is.null(out$df))
					if(isNA==FALSE){
						gradient <- out$opt$max_gradient <= max_gradient
						pdHess <- out$Sdreport$pdHess
					}	
				}	

				# if(all(is.null(out$df))==FALSE & (gradient == FALSE | pdHess == FALSE)){
				# 	out <- get_converged(results=out, saveFlagsDir=iterpath, saveFlagsName=paste0(modname, "_node_", x))
				# }

				## flag non-convergence or NAs
				if(all(is.null(out$df))){
					write("model NA", file.path(iterpath, paste0("modelNA_node_", x, "_", taxon, "_", dim, ".txt")))
				}
				if(all(is.null(out$df))==FALSE){
					gradient <- out$opt$max_gradient <= max_gradient
					pdHess <- out$Sdreport$pdHess
					if(gradient==FALSE){
						write("highgradient", file.path(iterpath,paste0("highgradient_node_", x, "_", taxon, "_", dim, ".txt")))
					}
					if(pdHess==FALSE){
						write("Hessian not positive definite", file.path(iterpath, paste0("pdHess_node_", x, "_", taxon, "_", dim, ".txt")))
					}
					if(gradient == TRUE & pdHess == TRUE) saveRDS(out, file.path(iterpath, paste0("_res_node_", x, "_", taxon, "_", dim, ".rds")))	
				}
		}
		if(model=="LBSPR"){
				LB_lengths <- new("LB_lengths")
				LB_lengths@LMids <- input$mids
				LB_lengths@LData <- as.matrix(input$LF[,,1], ncol=Nyears)
				LB_lengths@Years <- as.numeric(rownames(input$LF))
				LB_lengths@NYears <- Nyears
				LB_lengths@L_units <- "cm"

					##----------------------------------------------------------------
					## Step 2: Specify biological inputs and parameter starting values
					##----------------------------------------------------------------
				LB_pars <- new("LB_pars")
				LB_pars@MK <- input$M/input$vbk
				LB_pars@Linf <- input$linf
				LB_pars@L50 <- input$ML50
				LB_pars@L95 <- ifelse(is.na(input$ML95), input$linf, input$ML95)
				LB_pars@Walpha <- input$lwa
				LB_pars@Wbeta <- input$lwb
				LB_pars@BinWidth <- input$binwidth	
				LB_pars@SL50 <- input$SL50
				LB_pars@SL95 <- input$SL95
				LB_pars@R0 <- input$R0
				LB_pars@Steepness <- ifelse(input$h==1, 0.99, input$h)
				LB_pars@BinMin <- 0
				LB_pars@BinMax <- input$linf * 1.3

				lbspr_res <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths)
				out <- lbspr_res			
		}
					
				return(out)
		})
			saveRDS(res, file.path(iterpath, paste0("res_stacking_", model, "_", taxon, "_", dim, ".rds")))
			files <- list.files(path=file.path(iterpath))
			remove <- files[grepl(paste0("res_node"), files)]
			ignore <- sapply(1:length(remove), function(x) unlink(file.path(iterpath, remove[x]), TRUE))
	}


	return(paste0("Ran iter ", iter, " in ", savedir))
}
