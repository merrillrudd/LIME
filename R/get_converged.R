#' Check convergence and re-run for common issues
#'
#' \code{get_converged} re-run LIME by pre-identified adjustments if not converged
#'
#' @author M.B. Rudd
#' @param results list of LIME results from `run_LIME`
#' @param  max_gradient maximum gradient, default=0.001
#' @param saveFlagsDir directory to save flags for what was changed to get model converged
#' @param saveFlagsName name of model or node to save flags for what was changed or get model converged

#' @useDynLib LIME

#' @return prints how many iterations were run in model directory
#' 
#' @export
get_converged <- function(results, max_gradient=0.001, saveFlagsDir=FALSE, saveFlagsName=FALSE){

		## differentiate results input from results output
		out <- results
		out_save <- results
		if(all(is.null(out$df))) stop("Model from results list is NA, check model inputs and structure.")

		## model inputs for re-running that won't change between runs
		input <- out$input
		data_avail <- out$data_avail
		C_type <- out$Inputs$Data$C_type

		## identify convergence problems
		gradient <- out$opt$max_gradient <= max_gradient
		pdHess <- out$Sdreport$pdHess

					## check and rerun in case of nonconvergence, try to address multiple possible issues and rerun 2 times
					try <- 0
					fix_more <- FALSE
					est_selex_f <- TRUE
					while(try < 3 & all(is.null(out_save$df))==FALSE & (gradient == FALSE | pdHess == FALSE)){
						## first check that theta is not estimated extremely high
						## often a problem that theta is estimated very large, and high final gradient is on selectivity
						## more important to estimate selectivity and fix theta at a high number
						try <- try + 1
						print(try)
						if(out$Report$theta > 50){
							input$theta <- 50
							if(all(fix_more != FALSE)) fix_more <- c(fix_more, "log_theta")
							if(all(fix_more == FALSE)) fix_more <- "log_theta"
							out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, C_type=C_type, rewrite=TRUE, newtonsteps=3, fix_more=unique(fix_more), est_selex_f=est_selex_f)							
							## check_convergence
							isNA <- all(is.null(out$df))
							if(isNA) out <- out_save
							if(isNA==FALSE){
								out_save <- out
								gradient <- out$opt$max_gradient <= max_gradient
								pdHess <- out$Sdreport$pdHess
							}						
						}

						if(pdHess==FALSE){
							find_param <- unique(rownames(summary(out$Sdreport))[which(is.na(summary(out$Sdreport)[,2]))])
							find_param_est <- find_param[which(find_param %in% names(out$opt$par))]
							if("log_sigma_R" %in% find_param_est){
								input$SigmaR <- min(2, out$Report$sigma_R)
								if(all(fix_more == FALSE)) fix_more <- "log_sigma_R"
								if(all(fix_more != FALSE)) fix_more <- c(fix_more, "log_sigma_R")

								out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, C_type=C_type, rewrite=TRUE, newtonsteps=3, fix_more=unique(fix_more), est_selex_f=est_selex_f)
								
								## check_convergence
								isNA <- all(is.null(out$df))
								if(isNA) out <- out_save
								if(isNA==FALSE){
									out_save <- out
									gradient <- out$opt$max_gradient <= max_gradient
									pdHess <- out$Sdreport$pdHess
								}	
							}
						}

						if(pdHess==FALSE){
							find_param <- unique(rownames(summary(out$Sdreport))[which(is.na(summary(out$Sdreport)[,2]))])
							find_param_est <- find_param[which(find_param %in% names(out$opt$par))]
							if("log_S50_f" %in% find_param_est){
								if(try == 1){
									input$SL50 <- out$Report$S50
									input$SL95 <- out$Report$S95 * 1.3									
								}
								if(try > 1) est_selex_f <- FALSE
								out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, C_type=C_type, rewrite=TRUE, newtonsteps=3, fix_more=unique(fix_more), est_selex_f=est_selex_f)

								## check_convergence
								isNA <- all(is.null(out$df))
								if(isNA) out <- out_save
								if(isNA==FALSE){
									out_save <- out
									gradient <- out$opt$max_gradient <= max_gradient
									pdHess <- out$Sdreport$pdHess
								}	
							}
						}

						if(pdHess==FALSE){
							find_param <- unique(rownames(summary(out$Sdreport))[which(is.na(summary(out$Sdreport)[,2]))])
							find_param_est <- find_param[which(find_param %in% names(out$opt$par))]
							if("log_F_ft" %in% find_param_est){

									input$SL50 <- out$Report$S50
									input$SL95 <- out$Report$S95
									est_selex_f <- FALSE

								out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, C_type=C_type, rewrite=TRUE, newtonsteps=3, fix_more=unique(fix_more), est_selex_f=est_selex_f, f_startval_ft=matrix(mean(out$Report$F_ft), nrow=nrow(out$Report$F_ft), ncol=ncol(out$Report$F_ft)))
							
								## check_convergence
								isNA <- all(is.null(out$df))
								if(isNA) out <- out_save
								if(isNA==FALSE){
									out_save <- out
									gradient <- out$opt$max_gradient <= max_gradient
									pdHess <- out$Sdreport$pdHess
								}	
							}							
						}

						if(pdHess==FALSE){
							find_param <- unique(rownames(summary(out$Sdreport))[which(is.na(summary(out$Sdreport)[,2]))])
							find_param_est <- find_param[which(find_param %in% names(out$opt$par))]
							if("log_F_ft" %in% find_param_est){

								if(try==1) input$SigmaF <- 0.1
								if(try>1) input$SigmaF <- 0.05

								out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, C_type=C_type, rewrite=TRUE, newtonsteps=3, fix_more=unique(fix_more), est_selex_f=est_selex_f, f_startval_ft=matrix(mean(out$Report$F_ft), nrow=nrow(out$Report$F_ft), ncol=ncol(out$Report$F_ft)))
							
								## check_convergence
								isNA <- all(is.null(out$df))
								if(isNA) out <- out_save
								if(isNA==FALSE){
									out_save <- out
									gradient <- out$opt$max_gradient <= max_gradient
									pdHess <- out$Sdreport$pdHess
								}
								## switch it back in case models need to be run again
								input$SigmaF <- 0.2	
							}							
						}

						if(gradient==FALSE){
							## fix parameter with high final gradient
							if(all(fix_more == FALSE)) fix_more <- as.character(out$df[,2][which(abs(out$df[,1])>=0.001)])
							if(all(fix_more != FALSE)) fix_more <- c(fix_more, as.character(out$df[,2][which(abs(out$df[,1])>=0.001)]))
							fix_more <- fix_more[-which(grepl("F_ft", fix_more))]
							out <- run_LIME(modpath=NULL, input=input, data_avail=data_avail, C_type=C_type, rewrite=TRUE, newtonsteps=3, fix_more=unique(fix_more), est_selex_f=est_selex_f, f_startval_ft=matrix(mean(out$Report$F_ft), nrow=nrow(out$Report$F_ft), ncol=ncol(out$Report$F_ft)))

								## check_convergence
								isNA <- all(is.null(out$df))
								if(isNA) out <- out_save
								if(isNA==FALSE){
									out_save <- out
									gradient <- out$opt$max_gradient <= max_gradient
									pdHess <- out$Sdreport$pdHess
								}	
						}
					}

		## save flags if model converged
		if(isNA==FALSE & (gradient == TRUE & pdHess == TRUE)){
		  if(saveFlagsDir!=FALSE){
			if(all(fix_more!=FALSE)) saveRDS(fix_more, file.path(saveFlagsDir, paste0(saveFlagsName, "_fix_more.rds")))
			if(est_selex_f==FALSE) saveRDS(out$Report$S_fl, file.path(saveFlagsDir, paste0(saveFlagsName, "_fixed_selectivity.rds")))
			if(out$Inputs$Parameters$log_sigma_F != results$Inputs$Parameters$log_sigma_F) saveRDS(out$Report$sigma_F, file.path(saveFlagsDir, paste0(saveFlagsName, "_adjustSigmaF.rds")))
		  }
		  if(saveFlagsDir==FALSE){
		 	out$flags <- NULL
			if(all(fix_more!=FALSE)) out$flags <- c(out$flags, "fix_more")
			if(est_selex_f==FALSE) out$flags <- c(out$flags, "fixed_selectivity")
			if(out$Inputs$Parameters$log_sigma_F != results$Inputs$Parameters$log_sigma_F) out$flags <- c(out$flags, "adjustSigmaF")
		  }
		}

	return(out)
}