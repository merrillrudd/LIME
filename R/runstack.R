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
#' @param model default="LIME", alternate = "LBSPR"
#' @param sim_model default="LIME", alterate = "LBSPR"


#' @useDynLib LIME

#' @return prints how many iterations were run in model directory
#' 
#' @export
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
					Fscenario=NULL,
					model="LIME",
					sim_model="LIME"){

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
				SigmaF_inp <- 0.001
				SigmaR_inp <- 0.001
				rho_inp <- 0
				CVlen_inp <- 0.07
				Fdynamics_inp <- "Constant"
			}
			if(Fscenario=="harvestdyn" | Fscenario==FALSE){
				SigmaF_inp <- lh$SigmaF
				SigmaR_inp <- lh$SigmaR
				rho_inp <- lh$rho
				CVlen_inp <- lh$CVlen_inp
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
									Frate=Frate,
									theta=10,
									h=h,
									CVlen=CVlen_inp,
									nfleets=nfleets))


			if(rewrite==TRUE | file.exists(file.path(iterpath, "True.rds"))==FALSE){
				## use seed + 1000 to generate data
				if(sim_model=="LIME"){
					data <- generate_data(modpath=savedir, itervec=iter, 
									Fdynamics=Fdynamics_inp, Rdynamics="Constant", 
									lh=plist, 
									Nyears=10, Nyears_comp=10, comp_sample=200,
									init_depl=c(0.3,0.9), 
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
            		simlf <- rmultinom(10, size=200, prob=sim@pLCatch)

            		data <- list()
            		data$LF <- t(simlf)
            		colnames(data$LF) <- sim@LMids
            		rownames(data$LF) <- 1:10
            		data$mids <- sim@LMids
            		data$SPR <- sim@SPR
            		data$FM <- sim@FM
            		data$D_t <- sim@SSB/sim@SSB0
            		data$years <- 1:10
            		data$SL50 <- sim@SL50
            		data$SL95 <- sim@SL95
            		data$S_fl <- matrix((1 /(1 + exp(-log(19)*(sim@LMids-sim@SL50)/(sim@SL95-sim@SL50)))), nrow=1)
            		saveRDS(data, file.path(iterpath, "True.rds"))
				}

			}
			data <- readRDS(file.path(iterpath, "True.rds"))
			input_data <- list("years"=data$years, "LF"=data$LF)

		## run at true values
		if(rewrite==TRUE | file.exists(file.path(iterpath, paste0(modname, "_res_IterTrue_", model, ".rds")))==FALSE){	

			input <- create_inputs(lh=plist, input_data=input_data)

			if(model=="LIME"){
				## input file and run model
				if(nrow(input_data$LF)==1) Rdet <- TRUE
				if(nrow(input_data$LF)>1) Rdet <- FALSE
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
					if(gradient == TRUE & pdHess == TRUE) saveRDS(out, file.path(iterpath, paste0(modname, "_res_IterTrue_LIME.rds")))	
				}
			}
			if(model=="LBSPR"){
				LB_lengths <- new("LB_lengths")
				LB_lengths@LMids <- plist$mids
				# LB_lengths@LData <- as.matrix(input$LF[nrow(input$LF),,1], ncol=1)
				LB_lengths@LData <- t(input$LF[,,1])
				LB_lengths@Years <- as.numeric(rownames(input$LF))
				# LB_lengths@Years <- as.numeric(rownames(input$LF)[length(rownames(input$LF))])
				LB_lengths@NYears <- nrow(input$LF[,,1])	
				# LB_lengths@NYears <- 1	
				LB_lengths@L_units <- "cm"

					##----------------------------------------------------------------
					## Step 2: Specify biological inputs and parameter starting values
					##----------------------------------------------------------------
				LB_pars <- new("LB_pars")
				LB_pars@MK <- plist$M/plist$vbk
				LB_pars@Linf <- plist$linf
				LB_pars@L50 <- plist$ML50
				LB_pars@L95 <- plist$ML95
				LB_pars@Walpha <- plist$lwa
				LB_pars@Wbeta <- plist$lwb
				LB_pars@BinWidth <- plist$binwidth	
				LB_pars@SL50 <- data$SL50
				LB_pars@SL95 <- data$SL95
				LB_pars@R0 <- plist$R0
				LB_pars@Steepness <- ifelse(plist$h==1, 0.99, plist$h)
				LB_pars@L_units <- "cm"

				lbspr_res <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths, Control=list("GTG"))
				saveRDS(lbspr_res, file.path(iterpath, paste0(modname, "_res_IterTrue_LBSPR.rds")))	

			}
		}
	}

	## run at means from FishLife for ensemble parameters
	if(rewrite==TRUE | file.exists(file.path(iterpath, paste0(modname, "_res_FishLifeMeans_", model,".rds")))==FALSE){	

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
									theta=10,
									h=h,
									CVlen=CVlen,
									nfleets=nfleets))		

		if(model=="LIME"){
			## input file and run model
			input <- create_inputs(lh=lhinp, input_data=input_data)
			out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=FALSE, C_type=C_type, LFdist=LFdist)	

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
					out <- get_converged(results=out, saveFlagsDir=iterpath, saveFlagsName=paste0(modname, "_FishLifeMeans"))
				}

				## flag non-convergence or NAs
				if(all(is.null(out$df))){
					write("model NA", file.path(iterpath, paste0(modname, "_modelNA_FishLifeMeans.txt")))
				}
				if(all(is.null(out$df))==FALSE){
					gradient <- out$opt$max_gradient <= max_gradient
					pdHess <- out$Sdreport$pdHess
					if(gradient==FALSE){
						write("highgradient", file.path(iterpath, paste0(modname, "_highgradient_FishLifeMeans.txt")))
					}
					if(pdHess==FALSE){
						write("Hessian not positive definite", file.path(iterpath, paste0(modname, "_pdHess_FishLifeMeans.txt")))
					}
					## save results if converged
					if(gradient == TRUE & pdHess == TRUE) saveRDS(out, file.path(iterpath, paste0(modname, "_res_FishLifeMeans_LIME.rds")))	
				}
		}
		if(model=="LBSPR"){
				LB_lengths <- new("LB_lengths")
				LB_lengths@LMids <- plist$mids
				# LB_lengths@LData <- as.matrix(input$LF[nrow(input$LF),,1], ncol=1)
				LB_lengths@LData <- t(input$LF[,,1])
				LB_lengths@Years <- as.numeric(rownames(input$LF))
				# LB_lengths@Years <- as.numeric(rownames(input$LF)[length(rownames(input$LF))])
				LB_lengths@NYears <- nrow(input$LF[,,1])	
				# LB_lengths@NYears <- 1	
				LB_lengths@L_units <- "cm"

					##----------------------------------------------------------------
					## Step 2: Specify biological inputs and parameter starting values
					##----------------------------------------------------------------
				LB_pars <- new("LB_pars")
				LB_pars@MK <- lhinp$M/lhinp$vbk
				LB_pars@Linf <- lhinp$linf
				LB_pars@L50 <- lhinp$ML50
				LB_pars@L95 <- lhinp$ML95
				LB_pars@Walpha <- lhinp$lwa
				LB_pars@Wbeta <- lhinp$lwb
				LB_pars@BinWidth <- lhinp$binwidth	
				LB_pars@SL50 <- data$SL50
				LB_pars@SL95 <- data$SL95
				LB_pars@R0 <- lhinp$R0
				LB_pars@Steepness <- ifelse(lhinp$h==1, 0.99, lhinp$h)


				lbspr_res <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths)
				saveRDS(lbspr_res, file.path(iterpath, paste0(modname, "_res_FishLifeMeans_LBSPR.rds")))	
		}
	}	
	

	## predictive stacking
	if(rewrite==TRUE | file.exists(file.path(iterpath, paste0(modname, "_res_stacking_", model, ".rds")))==FALSE){
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
										theta=10,
										h=h,
										CVlen=CVlen,
										nfleets=nfleets))			

		if(model=="LIME"){
			## input files and run model
			input <- create_inputs(lh=lhinp, input_data=input_data)
			out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=3, C_type=C_type, LFdist=LFdist)		

				## check_convergence
				isNA <- all(is.null(out$df))
				if(isNA==TRUE){
					## before entering loop, check:
					out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, rewrite=TRUE, newtonsteps=FALSE, C_type=C_type, LFdist=LFdist)
					isNA <- all(is.null(out$df))
				}
				if(isNA==FALSE){
					gradient <- out$opt$max_gradient <= max_gradient
					pdHess <- out$Sdreport$pdHess
				}	

				if(all(is.null(out$df))==FALSE & (gradient == FALSE | pdHess == FALSE)){
					out <- get_converged(results=out, saveFlagsDir=iterpath, saveFlagsName=paste0(modname, "_node_", x))
				}

				## flag non-convergence or NAs
				if(all(is.null(out$df))){
					write("model NA", file.path(iterpath, paste0(modname, "_modelNA_node_", x, ".txt")))
				}
				if(all(is.null(out$df))==FALSE){
					gradient <- out$opt$max_gradient <= max_gradient
					pdHess <- out$Sdreport$pdHess
					if(gradient==FALSE){
						write("highgradient", file.path(iterpath,paste0(modname, "_highgradient_node_", x, ".txt")))
					}
					if(pdHess==FALSE){
						write("Hessian not positive definite", file.path(iterpath, paste0(modname, "_pdHess_node_", x, ".txt")))
					}
					if(gradient == TRUE & pdHess == TRUE) saveRDS(out, file.path(iterpath, paste0(modname, "_res_node_", x, ".rds")))	
				}
		}
		if(model=="LBSPR"){
				LB_lengths <- new("LB_lengths")
				LB_lengths@LMids <- plist$mids
				# LB_lengths@LData <- as.matrix(input$LF[nrow(input$LF),,1], ncol=1)
				LB_lengths@LData <- t(input$LF[,,1])
				LB_lengths@Years <- as.numeric(rownames(input$LF))
				# LB_lengths@Years <- as.numeric(rownames(input$LF)[length(rownames(input$LF))])
				LB_lengths@NYears <- nrow(input$LF[,,1])	
				# LB_lengths@NYears <- 1	
				LB_lengths@L_units <- "cm"

					##----------------------------------------------------------------
					## Step 2: Specify biological inputs and parameter starting values
					##----------------------------------------------------------------
				LB_pars <- new("LB_pars")
				LB_pars@MK <- lhinp$M/lhinp$vbk
				LB_pars@Linf <- lhinp$linf
				LB_pars@L50 <- lhinp$ML50
				LB_pars@L95 <- lhinp$ML95
				LB_pars@Walpha <- lhinp$lwa
				LB_pars@Wbeta <- lhinp$lwb
				LB_pars@BinWidth <- lhinp$binwidth	
				LB_pars@SL50 <- data$SL50
				LB_pars@SL95 <- data$SL95
				LB_pars@R0 <- lhinp$R0
				LB_pars@Steepness <- ifelse(lhinp$h==1, 0.99, lhinp$h)


				lbspr_res <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths)
				out <- lbspr_res			
		}
					
				return(out)
		})
			saveRDS(res, file.path(iterpath, paste0(modname, "_res_stacking_", model, ".rds")))
			files <- list.files(path=file.path(iterpath))
			remove <- files[grepl(paste0(modname,"_res_node"), files)]
			ignore <- sapply(1:length(remove), function(x) unlink(file.path(iterpath, remove[x]), TRUE))
	}
	

	return(paste0("Ran iter ", iter, " in ", savedir))
}
