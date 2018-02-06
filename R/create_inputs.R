#' Create input parameters for TMB model
#'
#' \code{create_inputs} Gets list of parameter inputs into the proper format
#'
#' @author M.B. Rudd
#' @param lh tagged list of life history/starting value information
#' @param input_data tagged list of data inputs. Required: years = vector of years (true years or indices); LF = matrix of length frequency (years along rows and length bins along columns), obs_per_year = vector of sample size per year. Optional: I_t = vector of abundance index, named with years; C_t = vector of catch, named with years. 

#' @return List, a tagged list of potentially useful benchmarks
#' @details Specifically used to merge life history information with other model settings; flexibility to change parameter inputs for sensitivity analysis without changing the baseline life history information that was used to generate data in a simulation study, or carefully compiled for real=life application
#' @export
create_inputs <- function(lh, input_data){
    
        ## copy life history
        dat_input <- c(lh, input_data)

        ## have the log ready in the input file for some variance parameterss
        dat_input$log_sigma_C <- log(dat_input$SigmaC)
        dat_input$log_sigma_I <- log(dat_input$SigmaI)
        
        dat_input$log_CV_L <- log(dat_input$CVlen)

        length_raw <- dat_input$LF
        bw <- dat_input$binwidth

        if(is.array(dat_input$LF)){
            length_raw <- dat_input$LF
            bins_dim <- seq(bw, by=bw, length=dim(length_raw)[2])
            max_bin <- max(c(max(dat_input$highs), 
                            max(bins_dim), 
                            sapply(1:dat_input$nfleets, function(x) as.numeric(colnames(length_raw[,,x])[max(which(colSums(length_raw[,,x])>0))]))))
            highs <- seq(bw, max_bin, by=bw)
            mids <- seq(bw/2, max(highs), by=bw)
            lows <- highs - bw
            time <- as.numeric(unique(c(sapply(1:length(length_raw), function(x) unique(rownames(length_raw))))))
            if(max_bin > max(bins_dim)){
                LF <- array(NA, dim=c(length(time), length(highs), dat_input$nfleets))            
                rownames(LF) <- time
                colnames(LF) <- highs
                for(f in 1:dat_input$nfleets){
                    LFsub <- length_raw[,,f]
                    LF[which(rownames(LFsub) %in% rownames(LF)),,f] <- LFsub
                }
            }
            if(max_bin <= max(bins_dim)) LF <- length_raw
            dat_input$LF <- LF
        }

        ## make sure length bins from life history and observed data match
        if(is.data.frame(dat_input$LF)){
            max_bin <- max(c(max(dat_input$highs), ceiling(max(length_raw$Value)*1.25)))
            highs <- seq(bw, max_bin, by=bw)
            mids <- seq(bw/2, max(highs), by=bw)
            lows <- highs - bw
            time <- unique(length_raw$X)[order(unique(length_raw$X))]
            LF <- array(NA, dim=c(length(time), length(highs), dat_input$nfleets))
            for(f in 1:dat_input$nfleets){
                lfind <- length_raw %>% filter(Fleet==f)
                lfreq <- t(sapply(1:length(time), function(x){
                    sub <- lfind %>% filter(X==time[x])
                    if(nrow(sub)>0){
                        out <- sapply(1:length(highs), function(y){
                            sub2 <- sub$Value[which(sub$Value > highs[y]-bw & sub$Value <= highs[y])]
                            return(length(sub2))
                        })
                    }
                    if(nrow(sub)==0) out <- rep(0, length(highs))
                    return(out)
                })) 

                LF[,,f] <- lfreq
            }
            rownames(LF) <- time
            colnames(LF) <- highs
            dat_input$LF <- LF
        }

        if(is.list(dat_input$LF)){
            for(f in 1:length(length_raw)){
                colnames(length_raw[[x]]) <- seq(bw, by=bw, length=ncol(length_raw[[x]]))
            }
            max_bin <- max(c(max(dat_input$highs), 
                            sapply(1:length(length_raw), function(x) max(as.numeric(colnames(length_raw[[x]])))), 
                            sapply(1:length(length_raw), function(x) as.numeric(colnames(length_raw[[x]])[max(which(colSums(length_raw[[x]])>0))]))))
            highs <- seq(bw, max_bin, by=bw)
            mids <- seq(bw/2, max(highs), by=bw)
            lows <- highs - bw
            time <- as.numeric(unique(c(sapply(1:length(length_raw), function(x) unique(rownames(length_raw[[x]]))))))
            LF <- array(NA, dim=c(length(time), length(highs), dat_input$nfleets))
            rownames(LF) <- time
            colnames(LF) <- highs
            for(f in 1:dat_input$nfleets){
                LFsub <- length_raw[[f]]
                LF[which(rownames(LFsub) %in% rownames(LF)),,f] <- LFsub
            }
            dat_input$LF <- LF
        }
        dat_input$highs <- highs
        dat_input$mids <- mids
        dat_input$lows <- lows

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