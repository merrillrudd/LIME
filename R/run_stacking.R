#' Run stacking algorithm
#'
#' \code{run_stacking} runs stacking algorithm for length-based stock assessment methods given quadrature nodes
#' @author M.B. Rudd
#' @param modpath, the directory to save results
#' @param iter, the iteration for a simulation study, default=NULL for an assessment
#' @param lh, the life history list output from LIME::create_lh_list
#' @param model, "LIME" or "LBSPR"
#' @param nodes, data frame of nodes mapped to life history parameters, number of rows = number of nodes, number of columns = life history parameter K, Loo, M, Lm
#' @param dim, "1D" to label 1-d stacking, "2D" to label 2-d stacking
#' @param input_data, default=NULL for simulation study, data list following rules by LIME for regular assessment
#' @importFrom stats uniroot optimize

#' @return List, a tagged list of potentially useful benchmarks
#' @export
run_stacking <- function(modpath, iter=NULL, lh, model, nodes, dim, input_data=NULL){

	stop("Deprecated: please use function bioensembles::run_stacking.")

	if(all(is.null(iter))==FALSE){
		iterpath <- file.path(modpath, iter)
		if(file.exists(file.path(iterpath, "True.rds"))==FALSE) stop("Generated data not available in directory")
		true <- readRDS(file.path(iterpath, "True.rds"))
		input_data <- list("years"=true$years, "LF"=true$LF)
	}
	if(all(is.null(iter))){
		iterpath <- modpath
		if(all(is.null(input_data))) stop("iter=NULL, not a simulation study, requires input data")
	}

		byNode <- lapply(1:nrow(nodes), function(nn){	

			lh_new <- with(lh, create_lh_list(vbk=exp(nodes[nn,"K"]), linf=exp(nodes[nn,"Loo"]),
											M=exp(nodes[nn,"M"]),
											M50=exp(nodes[nn,"Lm"]), maturity_input="length",
											lwa=lwa, lwb=lwb,
											S50=exp(nodes[nn,"Lm"]), S95=min(exp(nodes[nn,"Loo"])*0.95, exp(nodes[nn,"Lm"])*1.2), selex_input="length",
											SigmaF=0.1, SigmaR=0.737, rho=rho, binwidth=binwidth, theta=theta))	

			input <- create_inputs(lh=lh_new, input_data=input_data)		

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

				out <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths)
				saveRDS(out, file.path(iterpath, paste0("res_node_LBSPR.rds")))
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
					write("model NA", file.path(iterpath, paste0("modelNA_node_", nn, "_LIME_", dim, ".txt")))
				}
				if(all(is.null(out$df))==FALSE){
					gradient <- out$opt$max_gradient <= 0.001
					pdHess <- out$Sdreport$pdHess
					if(gradient==FALSE){
						write("highgradient", file.path(iterpath, paste0("highgradient_node_", nn, "_LIME_", dim, ".txt")))
					}
					if(pdHess==FALSE){
						write("Hessian not positive definite", file.path(iterpath, paste0("pdHess_node_", nn, "_LIME_", dim, ".txt")))
					}
					## save results if converged
					if(gradient == TRUE & pdHess == TRUE) saveRDS(out, file.path(iterpath, paste0("res_node_", nn, "_LIME_", dim, ".rds")))
				}	
			}		

			return(out)
		})
		saveRDS(byNode, file.path(iterpath, paste0("res_stacking_", model, "_", dim, ".rds")))
		files <- list.files(path=file.path(iterpath))
		remove <- files[grepl(paste0("res_node"), files)]
		ignore <- sapply(1:length(remove), function(x) unlink(file.path(iterpath, remove[x]), TRUE))


}