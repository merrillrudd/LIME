#' Predictive stacking
#'
#' \code{runstack} run predicted stacking for life history parameters using LIME
#'
#' @author M.B. Rudd
#' @param savedir directory to save results
#' @param iter iteration of generated data
#' @param seed set seed
#' @param tmax maximum age
#' @param nodes matrix of nodes where each column is a different parameter and each row is a value from a distribution
#' @param param parameters (column names for nodes)
#' @param mean mean of each parameter value 
#' @param cov covariance matrix across parameters
#' @param modname model name to save in directory
#' @param input_data for LIME use with real data (not simulated)
#' @param Fscenario fishing mortality scenario to generate data
#' @param rewrite rewrite results?

#' @useDynLib LIME

#' @return prints how many iterations were run in model directory
#' 
#' @export
runstack <- function(savedir, iter, seed, tmax, nodes, param, mean, cov, modname, input_data, Fscenario, rewrite){

	iterpath <- file.path(savedir, iter)
	dir.create(iterpath, showWarnings=FALSE)

	## generate data
		## true values
	set.seed(seed) ## use this seed to draw true value
		if(rewrite==TRUE){
			files <- list.files(iterpath)
			ignore <- sapply(1:length(files), function(x) unlink(file.path(iterpath, files[x]), TRUE))
		}
		vbk_choose <- ifelse("K" %in% param, rlnorm(1, mean=mean["K"], sd=sqrt(cov["K","K"])), exp(mean["K"]))
		M_choose <- ifelse("M" %in% param, rlnorm(1, mean=mean["M"], sd=sqrt(cov["M","M"])), exp(mean["M"]))
		Linf_choose <- ifelse("Loo" %in% param, rlnorm(1, mean=mean["Loo"], sd=sqrt(cov["Loo","Loo"])), exp(mean["Loo"]))
		if(Fscenario=="equil"){
			SigmaF_inp <- 0.01
			SigmaR_inp <- 0.01
			rho_inp <- 0
			Fdynamics_inp <- "Constant"
		}
		if(Fscenario=="harvestdyn"){
			SigmaF_inp <- 0.2
			SigmaR_inp <- 0.6
			rho_inp <- 0.4
			Fdynamics_inp <- "Endogenous"
		}
		plist <- create_lh_list(linf=Linf_choose, vbk=vbk_choose, t0=-1.77,
								lwa=0.053, lwb=2.706,
								M=M_choose,
								M50=16.9, maturity_input="length",
								S50=16.9, S95=25, selex_input="length",
								SigmaF=SigmaF_inp, SigmaR=SigmaR_inp, rho=rho_inp,
								AgeMax=tmax)

	if(all(input_data==FALSE)){
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

		# LFlist <- NULL
		# for(f in 1:plist$nfleets){
		# 	LFlist[[f]] <- input_data$LF[,,f]
		# }
		# plot_LCfits(LFlist=LFlist, ylim=c(0,0.15))

	## run at means from FishLife
	if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("res_FishLifeMeans.rds")))==FALSE){

		if(file.exists(file.path(iterpath,"nonconvergence_FishLifeMeans.txt"))) unlink(file.path(iterpath, "nonconvergence_FishLifeMeans.txt"), TRUE)
		if(file.exists(file.path(iterpath,"modelNA_FishLifeMeans.txt"))) unlink(file.path(iterpath, "modelNA_FishLifeMeans.txt"), TRUE)
		if(file.exists(file.path(iterpath,"res__FishLifeMeans.txt"))) unlink(file.path(iterpath, "res__FishLifeMeans.txt"), TRUE)

				vbk_inp <- exp(mean["K"])
				M_inp <- exp(mean["M"])
				linf_inp <- exp(mean["Loo"])
		 		lhinp <- with(plist, 
		 				 create_lh_list(linf=linf_inp, vbk=vbk_inp, t0=t0,
									lwa=lwa, lwb=lwb,
									M=M_inp,
									M50=M50, maturity_input="length",
									S50=SL50, S95=SL95, selex_input="length",
									SigmaF=0.1, SigmaR=SigmaR,
									AgeMax=AgeMax))		
		 		input <- create_inputs(lh=lhinp, input_data=input_data)
		 		out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=3)	
		 		if(max(abs(out$df[,1]))>0.001) write("nonconvergence", file.path(iterpath,"nonconvergence_FishLifeMeans.txt"))
		 		if(all(is.null(out$df))) write("model NA", file.path(iterpath, "modelNA_FishLifeMeans.txt"))
		 		if(max(abs(out$df[,1]))<=0.001) saveRDS(out, file.path(iterpath, paste0("res_FishLifeMeans.rds")))

		 		# plot_output(Inputs=out$Inputs, Report=out$Report, Sdreport=out$Sdreport, plot=c("Fish","Rec","SPR","SB"), True=data)
	}

	## run at true values
	if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("res_IterTrue.rds")))==FALSE){	
		
		if(file.exists(file.path(iterpath,"nonconvergence_IterTrue.txt"))) unlink(file.path(iterpath, "nonconvergence_IterTrue.txt"), TRUE)
		if(file.exists(file.path(iterpath,"modelNA_IterTrue.txt"))) unlink(file.path(iterpath, "modelNA_IterTrue.txt"), TRUE)
		if(file.exists(file.path(iterpath,"res__IterTrue.txt"))) unlink(file.path(iterpath, "res__IterTrue.txt"), TRUE)

		input <- create_inputs(lh=plist, input_data=input_data)
		input$SigmaF <- 0.1
		out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=3)	
		if(max(abs(out$df[,1]))>0.001) write("nonconvergence", file.path(iterpath,"nonconvergence_IterTrue.txt"))
		if(all(is.null(out$df))) write("model NA", file.path(iterpath, "modelNA_IterTrue.txt"))
		if(max(abs(out$df[,1]))<=0.001) saveRDS(out, file.path(iterpath, paste0("res_IterTrue.rds")))

		# plot_output(Inputs=out$Inputs, Report=out$Report, Sdreport=out$Sdreport, plot=c("Fish","Rec","SPR","SB"), True=data)

	}


	if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("res_", modname, ".rds")))==FALSE){
		res <- lapply(1:nrow(nodes), function(x){
				vbk_inp <- ifelse("K" %in% param, exp(nodes[x,"K"]), exp(mean["K"]))
				M_inp <- ifelse("M" %in% param, exp(nodes[x,"M"]), exp(mean["M"]))
				linf_inp <- ifelse("Loo" %in% param, exp(nodes[x,"Loo"]), exp(mean["Loo"]))
		 		lhinp <- with(plist, 
		 				 create_lh_list(linf=linf_inp, vbk=vbk_inp, t0=t0,
									lwa=lwa, lwb=lwb,
									M=M_inp,
									M50=M50, maturity_input="length",
									S50=SL50, S95=SL95, selex_input="length",
									SigmaF=0.1, SigmaR=SigmaR,
									AgeMax=AgeMax))		
				input <- create_inputs(lh=lhinp, input_data=input_data)
				out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=3)	
		 		if(max(abs(out$df[,1]))>0.001) write("nonconvergence", file.path(iterpath, paste0("nonconvergence_", modname, "node", x, ".txt")))
		 		if(all(is.null(out$df))) write("modelNA", file.path(iterpath, paste0("modelNA_", modname, "node", x, ".txt")))
		 		if(max(abs(out$df[,1]))<=0.001) saveRDS(out, file.path(iterpath, paste0("res_", modname, "node", x, ".txt")))
			return(out)
		})
		saveRDS(res, file.path(iterpath, paste0("res_", modname, ".rds")))
		files <- list.files(path=file.path(iterpath))
		remove <- files[grepl(paste0("res_", modname, "node"), files)]
		ignore <- sapply(1:length(remove), function(x) unlink(file.path(iterpath, remove[x]), TRUE))
	}
	return(paste0("Ran iter ", iter, " in ", savedir))

}
