#' Create input parameters for TMB model
#'
#' \code{create_inputs} Gets list of parameter inputs into the proper format
#'
#' @param param parameter name to adjust (sensitivity analysis), default=FALSE
#' @param val value of parameter name to adjust (sensitivity analysis), default=FALSE
#' @param lh tagged list of life history/starting value information
#' @param input_data tagged list of data inputs. Required: years = vector of years (true years or indices); LF = matrix of length frequency (years along rows and length bins along columns), obs_per_year = vector of sample size per year. Optional: I_t = vector of abundance index, named with years; C_t = vector of catch, named with years. 

#' @return List, a tagged list of potentially useful benchmarks
#' @details Specifically used to merge life history information with other model settings; flexibility to change parameter inputs for sensitivity analysis without changing the baseline life history information that was used to generate data in a simulation study, or carefully compiled for real=life application
#' @export
create_inputs <- function(param=FALSE, val=FALSE, lh, input_data){
    
        ## copy life history
        dat_input <- c(lh, input_data)

        ## change input values for sensitivity analysis
        if(param[1]!=FALSE){
            for(pp in 1:length(param)){
                dat_input[[param[pp]]] <- val[which(param==param[pp])]
            }
        }

        ## have the log ready in the input file for some variance parameterss
        dat_input$log_sigma_C <- log(dat_input$SigmaC)
        dat_input$log_sigma_I <- log(dat_input$SigmaI)
        
        dat_input$log_CV_L <- log(dat_input$CVlen)

        ## make sure length bins from life history and observed data match
        if(ncol(dat_input$LF)!=length(dat_input$highs)){
          if(is.null(colnames(dat_input$LF))==FALSE) obs_lb <- as.numeric(names(which(rev(colSums(dat_input$LF))>0)[1]))
          if(is.null(colnames(dat_input$LF))) obs_lb <- max(seq(dat_input$binwidth,length=ncol(dat_input$LF), by=dat_input$binwidth))
            if(obs_lb < max(dat_input$highs)){
                new <- matrix(0, nrow=nrow(dat_input$LF), ncol=length(seq(obs_lb+dat_input$binwidth, max(dat_input$highs), by=dat_input$binwidth)))
                rownames(new) <- rownames(dat_input$LF)
                colnames(new) <- seq(obs_lb+dat_input$binwidth, max(dat_input$highs), by=dat_input$binwidth)
                if(nrow(dat_input$LF)==1) LF_new <- cbind(t(as.matrix(dat_input$LF[,1:which(colnames(dat_input$LF)==obs_lb)])), new)
                if(nrow(dat_input$LF)>1) LF_new <- cbind(dat_input$LF[,1:which(colnames(dat_input$LF)==obs_lb)], new)
                dat_input$LF <- as.matrix(LF_new)
            }
            if(obs_lb >= max(dat_input$highs)){
                max_lb <- max(from=seq(dat_input$binwidth, to=ncol(dat_input$LF), by=dat_input$binwidth))
                test_lb <- max(seq(from=(obs_lb + dat_input$binwidth), length=5, by=dat_input$binwidth))
                change_lb <- min(test_lb, max_lb)
                if(is.null(colnames(dat_input$LF))==FALSE) index_lb <- which(colnames(dat_input$LF)==change_lb)
                if(is.null(colnames(dat_input$LF))) index_lb <- which(seq(dat_input$binwidth, length=ncol(dat_input$LF), by=dat_input$binwidth)==change_lb)
                dat_input$LF <- dat_input$LF[,1:index_lb]
                dat_input$highs <- seq(dat_input$binwidth, change_lb, by=dat_input$binwidth)
                dat_input$mids <- seq(dat_input$binwidth/2, change_lb-dat_input$binwidth/2, by=dat_input$binwidth)
                dat_input$lows <- dat_input$highs - dat_input$binwidth
            }
        }

        LF <- dat_input$LF

        years_i <- seq_along(dat_input$years)
        years_o <- which(dat_input$years %in% as.numeric(rownames(LF)))
        years_oi <- which(years_i %in% as.numeric(rownames(LF)))
        if(length(years_oi)>0) years_name <- years_oi
        if(length(years_o)>0) years_name <- years_o
        rownames(LF) <- years_name
        dat_input$LF <- LF

        if(is.null(dat_input$Nyears)) dat_input$Nyears <- length(dat_input$years)


    return(dat_input)
}