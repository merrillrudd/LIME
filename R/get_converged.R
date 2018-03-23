#' Check convergence and re-run for common issues
#'
#' \code{get_converged} re-run LIME by pre-identified adjustments if not converged
#'
#' @author M.B. Rudd
#' @param results list of LIME results from `run_LIME`
#' @param  max_gradient maximum gradient, default=0.001

#' @useDynLib LIME

#' @return prints how many iterations were run in model directory
#' 
#' @export
get_converged <- function(results, max_gradient=0.001){

		## differentiate results input from results output
		out <- results
		if(all(is.null(out$df))) stop("Model from results list is NA, check model inputs and structure.")

		## model inputs for re-running that won't change between runs
		input <- out$input
		data_avail <- out$data_avail
		C_type <- out$Inputs$Data$C_type

		## identify convergence problems
		gradient <- out$opt$max_gradient <= max_gradient
		pdHess <- out$Sdreport$pdHess

					## check and rerun in case of nonconvergence, try maximum 3 times
					try <- 0
					while(try < 3 & all(is.null(out$df))==FALSE & (gradient == FALSE | pdHess == FALSE)){
						## first check that theta is not estimated extremely high
						## often a problem that theta is estimated very large, and high final gradient is on selectivity
						## more important to estimate selectivity and fix theta at a high number
						try <- try + 1
						print(try)
						if(out$Report$theta > 50){
							input$theta <- 50
							out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, C_type=C_type, rewrite=TRUE, newtonsteps=3, fix_more="log_theta")
							
							## check_convergence
							isNA <- all(is.null(out$df))
							if(isNA==FALSE){
								gradient <- out$opt$max_gradient <= max_gradient
								pdHess <- out$Sdreport$pdHess
							}						
						}

						if(pdHess==FALSE){
							find_param <- unique(rownames(summary(out$Sdreport))[which(is.na(summary(out$Sdreport)[,2]))])
							find_param_est <- find_param[which(find_param %in% names(out$opt$par))]
							if("log_sigma_R" %in% find_param_est){
								input$SigmaR <- out$Report$sigma_R
								out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, C_type=C_type, rewrite=TRUE, newtonsteps=3, fix_more="log_sigma_R")
								
								## check_convergence
								isNA <- all(is.null(out$df))
								if(isNA==FALSE){
									gradient <- out$opt$max_gradient <= max_gradient
									pdHess <- out$Sdreport$pdHess
								}	
							}
						}

						if(pdHess==FALSE){
							find_param <- unique(rownames(summary(out$Sdreport))[which(is.na(summary(out$Sdreport)[,2]))])
							find_param_est <- find_param[which(find_param %in% names(out$opt$par))]
							if("log_S50_f" %in% find_param_est){
								input$SL50 <- out$Report$S50
								input$SL95 <- out$Report$S95 * 1.3
								out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, C_type=C_type, rewrite=TRUE, newtonsteps=3)

								## check_convergence
								isNA <- all(is.null(out$df))
								if(isNA==FALSE){
									gradient <- out$opt$max_gradient <= max_gradient
									pdHess <- out$Sdreport$pdHess
								}	
							}
						}

						if(pdHess==FALSE){
							find_param <- unique(rownames(summary(out$Sdreport))[which(is.na(summary(out$Sdreport)[,2]))])
							find_param_est <- find_param[which(find_param %in% names(out$opt$par))]
							if("log_F_ft" %in% find_param_est){
								out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, C_type=C_type, rewrite=TRUE, newtonsteps=3, f_startval_ft=matrix(mean(out$Report$F_ft), nrow=nrow(out$Report$F_ft), ncol=ncol(out$Report$F_ft)))
							
								## check_convergence
								isNA <- all(is.null(out$df))
								if(isNA==FALSE){
									gradient <- out$opt$max_gradient <= max_gradient
									pdHess <- out$Sdreport$pdHess
								}	
							}							
						}

						if(gradient==FALSE){
							## fix parameter with high final gradient
							find_param <- as.character(out$df[,2][which(abs(out$df[,1])>=0.001)])
							out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, C_type=C_type, rewrite=TRUE, newtonsteps=3, fix_more=find_param)

								## check_convergence
								isNA <- all(is.null(out$df))
								if(isNA==FALSE){
									gradient <- out$opt$max_gradient <= max_gradient
									pdHess <- out$Sdreport$pdHess
								}	
						}
					}

	outs <- NULL
	outs$out <- out
	outs$gradient <- gradient
	outs$pdHess <- pdHess
	outs$isNA <- isNA
	return(outs)
}