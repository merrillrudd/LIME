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
            bins_dim <- as.numeric(colnames(length_raw))
            # if(bins_dim[1] != bw){
            #     ## add zeros to the beginning of the length comps
            #     bins_pre <- rev(seq(from=min(bins)-bw, to=bw, by=-bw))
            #     add <- matrix(0, nrow=nrow(lf), ncol=length(bins_pre))
            #     colnames(add) <- bins_pre           

            #     length_raw <- cbind(add, length_raw)      
            #     bins_dim <- as.numeric(colnames(length_raw))          
            # }

            if(is.matrix(length_raw)){
                new <- array(NA, dim=c(dim(length_raw),1))
                new[,,1] <- length_raw
                rownames(new) <- rownames(length_raw)
                colnames(new) <- colnames(length_raw)
                length_raw <- new
            }
            if(dim(length_raw)[1] > 1) max_bin <- max(c(max(bins_dim)+(bw/2), 
                                                        sapply(1:dat_input$nfleets, function(x) as.numeric(bins_dim[max(which(colSums(length_raw[,,x])>0))]) + (bw/2))))
            if(dim(length_raw)[1] == 1) max_bin <- max(c(max(bins_dim)+(bw/2),
                                                        sapply(1:dat_input$nfleets, function(x) as.numeric(bins_dim[max(which(length_raw[,,x]>0))]) + (bw/2))))
        
            mids <- seq(bins_dim[1], max_bin - (bw/2), by=bw)
            highs <- seq(bins_dim[1] + (bw/2), max_bin, by=bw)
            lows <- highs - bw
            time <- dat_input$years
            if(max_bin > (max(bins_dim) + (bw/2))){
                LF <- array(0, dim=c(length(time), length(highs), dat_input$nfleets))            
                rownames(LF) <- time
                colnames(LF) <- mids
                for(f in 1:dat_input$nfleets){
                    LFsub <- matrix(length_raw[,,f], nrow=nrow(length_raw), ncol=ncol(length_raw))
                    colnames(LFsub) <- bins_dim
                    rownames(LFsub) <- rownames(length_raw)
                    LF[which(rownames(LF) %in% rownames(LFsub)),which(colnames(LF) %in% colnames(LFsub)),f] <- LFsub
                }
            }
            if(max_bin <= (max(bins_dim) + (bw/2))) LF <- length_raw
            if(dim(LF)[1] != length(time)){
                LF <- array(0, dim=c(length(time), length(highs), dat_input$nfleets))            
                rownames(LF) <- time
                colnames(LF) <- mids
                for(f in 1:dat_input$nfleets){
                    LFsub <- matrix(length_raw[,,f], nrow=nrow(length_raw), ncol=ncol(length_raw))
                    colnames(LFsub) <- bins_dim
                    rownames(LFsub) <- rownames(length_raw)
                    LF[which(rownames(LF) %in% rownames(LFsub)),which(colnames(LF) %in% colnames(LFsub)),f] <- LFsub
                }                
            }
            dat_input$LF <- LF
        }

        ## make sure length bins from life history and observed data match
        if(is.data.frame(dat_input$LF)){
            stop("Input data as length frequency matrix with years along rows and length bins along columns. For multiple fleets, input array where 3rd dimension is by fleet, or list with matrices of equal dimension and zeros for no observations. ")
            # max_bin <- max(c(ceiling(max(length_raw$Value)*1.25)))
            # mids <- seq(bins_dim[1], max_bin - (bw/2), by=bw)
            # highs <- seq(bins_dim[1] + (bw/2), max_bin, by=bw)
            # lows <- highs - bw
            # time <- dat_input$years
            # check_dim <- dim(dat_input$LF)
            # if(check_dim[2]!=length(highs)) stop("Is your matrix cast as a data frame? Data frame input must be in long form, with column names 'Fleet', 'Time', and 'Value'")
            # LF <- array(NA, dim=c(length(time), length(highs), dat_input$nfleets))
            # for(f in 1:dat_input$nfleets){
            #     lfind <- length_raw %>% filter(Fleet==f)
            #     lfreq <- t(sapply(1:length(time), function(x){
            #         sub <- lfind %>% filter(Time==time[x])
            #         if(nrow(sub)>0){
            #             out <- sapply(1:length(highs), function(y){
            #                 sub2 <- sub$Value[which(sub$Value > highs[y]-bw & sub$Value <= highs[y])]
            #                 return(length(sub2))
            #             })
            #         }
            #         if(nrow(sub)==0) out <- rep(0, length(highs))
            #         return(out)
            #     })) 

            #     LF[,,f] <- lfreq
            # }
            # rownames(LF) <- time
            # colnames(LF) <- highs
            # dat_input$LF <- LF
        }

        if(is.list(dat_input$LF)){
            length_raw <- dat_input$LF
            bins_dim <- as.numeric(colnames(length_raw[[1]]))
            # if(bins_dim[1] != bw){
            #     ## add zeros to the beginning of the length comps
            #     bins_pre <- rev(seq(from=min(bins)-bw, to=bw, by=-bw))
            #     add <- matrix(0, nrow=nrow(lf), ncol=length(bins_pre))
            #     colnames(add) <- bins_pre           

            #     for(f in 1:length(length_raw)){
            #         length_raw[[f]] <- cbind(add, length_raw[[f]])      
            #     }
            #     bins_dim <- as.numeric(colnames(length_raw[[1]]))          
            # }
            max_bin <- max(c(sapply(1:length(length_raw), function(x) max(as.numeric(colnames(length_raw[[x]]))) + (bw/2)), 
                            sapply(1:length(length_raw), function(x) as.numeric(colnames(length_raw[[x]])[max(which(colSums(length_raw[[x]])>0))]) + (bw/2))))
            mids <- seq(bins_dim[1], max_bin - (bw/2), by=bw)
            highs <- seq(bins_dim[1] + (bw/2), max_bin, by=bw)            
            lows <- highs - bw
            time <- dat_input$years
            LF <- array(0, dim=c(length(time), length(highs), dat_input$nfleets))
            rownames(LF) <- time
            colnames(LF) <- mids
            for(f in 1:dat_input$nfleets){
                LFsub <- length_raw[[f]]
                colnames(LFsub) <- bins_dim
                rownames(LFsub) <- rownames(length_raw[[f]])
                LF[which(rownames(LF) %in% rownames(LFsub)),which(colnames(LF) %in% colnames(LFsub)),f] <- LFsub
            }

            dat_input$LF <- LF
        }

        if(max(mids) < lh$linf){
            add_vals <- seq(from=max(mids)+bw, to=lh$linf*1.3, by=bw)
            LF_new <- array(0, dim=c(nrow(dat_input$LF), (ncol(dat_input$LF)+length(add_vals)), lh$nfleets))
            for(i in 1:lh$nfleets){
                LF_new[,1:ncol(dat_input$LF),i] <- dat_input$LF[,,i]
            }
            colnames(LF_new) <- c(colnames(dat_input$LF), add_vals)
            rownames(LF_new) <- rownames(dat_input$LF)
            dat_input$LF <- LF_new
        }

        dat_input$mids <- as.numeric(colnames(dat_input$LF))
        dat_input$highs <- dat_input$mids + (bw/2)
        dat_input$lows <- dat_input$mids - (bw/2)

        years_i <- seq_along(dat_input$years)
        years_o <- which(dat_input$years %in% as.numeric(rownames(LF)))
        years_oi <- which(years_i %in% as.numeric(rownames(LF)))
        if(length(years_oi)>0) years_name <- years_oi
        if(length(years_o)>0) years_name <- years_o
        rownames(LF) <- years_name
        dat_input$LF <- LF
        # dat_input$years_i <- years_i

        if(is.null(dat_input$Nyears)) dat_input$Nyears <- length(dat_input$years)


    return(dat_input)
}