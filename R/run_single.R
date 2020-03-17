run_single <- function(modpath, iter, lh, model, name){

	iterpath <- file.path(modpath, iter)
	if(file.exists(file.path(iterpath, "True.rds"))==FALSE) stop("Generated data not available in directory")
	true <- readRDS(file.path(iterpath, "True.rds"))

	input_data <- list("years"=true$years, "LF"=true$LF)
	input <- create_inputs(lh=lh, input_data=input_data)

	LFobs <- input$LF[,,1]
	obs <- which(rowSums(LFobs)>0)
	LFobs <- LFobs[obs,]
	Yobs <- input$years[obs]

	if(model=="LBSPR"){
		LB_lengths <- new("LB_lengths")
		LB_lengths@LMids <- input$mids
		LB_lengths@LData <- t(LFobs)
		LB_lengths@Years <- Yobs
		LB_lengths@NYears <- length(Yobs)
		LB_lengths@L_units <- "cm"

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
		saveRDS(lbspr_res, file.path(iterpath, paste0("res_", name, "_LBSPR.rds")))
	}

	if(model=="LIME"){
		input$SigmaR <- 0.737
		input$SigmaF <- 0.1
		out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=3, C_type=0, LFdist=1)

		## check_convergence
		isNA <- all(is.null(out$df))
		if(isNA==TRUE){
			out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=FALSE, C_type=0, LFdist=1)
		}
		isNA <- all(is.null(out$df))
		if(isNA==FALSE){
			gradient <- out$opt$max_gradient <= 0.001
			pdHess <- out$Sdreport$pdHess
			if(pdHess==FALSE){
				input$theta <- 50
				out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=3, C_type=0, LFdist=1, fix_more="log_theta")
			}
			if(gradient==FALSE){
				index <- out$df[which(out$df[,1]>0.001),2]
				if("log_S50_f" %in% index){
					input$SL50 <- out$Report$S50
					out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=FALSE, C_type=0, LFdist=1,fix_more="log_S50_f")
				}
			}
			isNA <- all(is.null(out$df))
			if(isNA==TRUE){
				out <- run_LIME(modpath=NULL, input=input, data_avail="LC", rewrite=TRUE, newtonsteps=FALSE, C_type=0, LFdist=1,fix_more="log_theta")
			}
			isNA <- all(is.null(out$df))
			if(isNA==FALSE){
				gradient <- out$opt$max_gradient <= 0.001
				pdHess <- out$Sdreport$pdHess
			}	
		}
		## flag non-convergence or NAs
		if(all(is.null(out$df))){
			write("model NA", file.path(iterpath, paste0("modelNA_", name, "_LIME.txt")))
		}
		if(all(is.null(out$df))==FALSE){
			gradient <- out$opt$max_gradient <= 0.001
			pdHess <- out$Sdreport$pdHess
			if(gradient==FALSE){
				write("highgradient", file.path(iterpath, paste0("highgradient_", name, "_LIME.txt")))
			}
			if(pdHess==FALSE){
				write("Hessian not positive definite", file.path(iterpath, paste0("pdHess_", name, "_LIME.txt")))
			}
			## save results if converged
			if(gradient == TRUE & pdHess == TRUE) saveRDS(out, file.path(iterpath, paste0("res_", name, "_LIME.rds")))	
		}	
	}
}