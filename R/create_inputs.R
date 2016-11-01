#' Create input parameters for TMB model
#'
#' \code{create_inputs} Gets list of parameter inputs into the proper format
#'
#' @param param parameter name to adjust (sensitivity analysis)
#' @param val value of parameter name to adjust (sensitivity analysis)
#' @param lh_list tagged list of life history/starting value information
#' @param data_avail_list artifact from sensitivity analysis, adjusts some information on sample size, years, etc. for different data-availability scenarios

#' @return List, a tagged list of potentially useful benchmarks
#' @details Specifically used to merge life history information with other model settings; flexibility to change parameter inputs for sensitivity analysis without changing the baseline life history information that was used to generate data in a simulation study, or carefully compiled for real=life application
#' @export
create_inputs <- function(param, val, lh_list, data_avail_list){
    
        ## copy life history
        dat_input <- c(lh_list, data_avail_list)

        ## have the log ready in the input file for some variance parameterss
        dat_input$log_sigma_C <- log(lh_list$SigmaC)
        dat_input$log_sigma_I <- log(lh_list$SigmaI)

        ## change input values for sensitivity analysis
        if(param!=FALSE){
            for(pp in 1:length(param)){
                dat_input[[param[pp]]] <- val[which(param==param[pp])]
            }
        }
        
        dat_input$log_CV_L <- log(dat_input$CVlen)

        ## make sure length bins from life history and observed data match
        obs_lb <- as.numeric(names(which(rev(colSums(dat_input$LF))>0)[1]))
        if(obs_lb <= max(dat_input$highs)){
            dat_input$LF <- dat_input$LF[,1:max(dat_input$highs)]
            dat_input$LFprop <- dat_input$LFprop[,1:max(dat_input$highs)]
        }
        if(obs_lb > max(dat_input$highs)){
            dat_input$LF <- dat_input$LF[,1:obs_lb]
            dat_input$LFprop <- dat_input$LFprop[,1:obs_lb]
            dat_input$highs <- c(dat_input$highs,seq(from=(max(dat_input$highs)+dat_input$binwidth), to=obs_lb, by=dat_input$binwidth))
            dat_input$mids <- c(dat_input$mids, seq(from=(max(dat_input$mids)+dat_input$binwidth), to=obs_lb, by=dat_input$binwidth))
            dat_input$lows <- c(dat_input$lows, seq(from=(max(dat_input$lows)+dat_input$binwidth), to=(obs_lb-dat_input$binwidth), by=dat_input$binwidth))
        }

    return(dat_input)
}