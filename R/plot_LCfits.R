#' plot length composition data and model fits
#'
#' \code{plot_LCfits} plot length composition data with option for model fits from LIME and LB-SPR
#'
#' @author M.B. Rudd
#' @param LF_df list of length frequency matrices
#' @param Inputs LIME input file; default NULL if only plotting length data
#' @param Report LIME report file; default NULL if only plotting length data
#' @param LBSPR LBSPR results - must have pLF = probability of being harvested in a length bin; default NULL
#' @param ylim ylim for plot; default NULL
#' @param dim dimensions of plot; default NULL if not specified, will approximate best based on number of years to plot
#' @param n if TRUE, will display sample size of length comp data; default FALSE
#' @param true_years optional vector of true years to label years of length comp
#' @importFrom graphics abline axis barplot box legend lines mtext par
#' 
#' @return figure with length composition data and model fits if Report or LBSPR are specified
#' 
#' @export
plot_LCfits <- function(LF_df=NULL, binwidth=1, Inputs=NULL, Report=NULL, LBSPR=NULL, ylim=NULL, dim=NULL, n=FALSE, true_years=NULL){
	# dev.new()

	if(all(is.null(Inputs))){
		lengths <- unique(LF_df$Length)[order(unique(LF_df$Length))]
		bw <- binwidth
		bins <- seq(bw, max(lengths), by=bw)
	}
	if(all(is.null(Inputs))==FALSE){
		LF_array <- Inputs$Data$LF_tlf
		LF_df <- LFreq_df(LF_array)
		bins <- as.numeric(colnames(LF_array))
		bw <- bins[1]
	}
	years <- unique(LF_df$Year)[order(unique(LF_df$Year))]
	nyears <- length(years)
	n_yr <- rep(0, nyears)
	for(i in 1:nyears){
		sub <- LF_df %>% dplyr::filter(Year==years[i])
		n_yr[i] <- nrow(sub)
	}

	nf <- length(unique(LF_df$Fleet))
	LCyrs <- lapply(1:nf, function(x){
		sub <- LF_df %>% dplyr::filter(Fleet==x)
		yrs <- unique(sub$Year)[order(unique(sub$Year))]
		return(yrs)
	})
	all_lc_years <- min(as.numeric(unlist(LCyrs))):max(as.numeric(unlist(LCyrs)))
	if(all(is.null(true_years))) true_years <- all_lc_years


	if(all(is.null(Inputs))) Tyrs <- all_lc_years
	if(all(is.null(Inputs))==FALSE) Tyrs <- 1:Inputs$Data$n_t

	if(all(is.null(Report))==FALSE){
		# pred <- Report$plb
			pred <- lapply(1:nf, function(x){
				return(Report$plb[,,x])
			})
	}
	if(all(is.null(Report))){
		pred <- NULL
	}
	if(all(is.null(LBSPR))==FALSE){
		pred2 <- t(LBSPR@pLCatch)
		rownames(pred2) <- which(rowSums(LFlist[[1]]) > 0)
		colnames(pred2) <- LBSPR@LMids
		LF2_new <- matrix(0, nrow=nrow(pred[[1]]), ncol=ncol(pred[[1]]))
		colnames(LF2_new) <- bins
		rownames(LF2_new) <- Tyrs
		for(i in 1:nrow(LF2_new)){
			if(Tyrs[i] %in% rownames(pred2)){
				for(j in 1:ncol(LF2_new)){
					for(k in 1:length(LBSPR@LMids)){
						if(j==1){
							if(LBSPR@LMids[k] <= bins[j]){	
								LF2_new[i,j] <- pred2[which(rownames(pred2)==Tyrs[i]),k]
							}
						}
						if(j>1){
							if(LBSPR@LMids[k] <= bins[j] & LBSPR@LMids[k] >= bins[j-1]){
								LF2_new[i,j] <- pred2[which(rownames(pred2)==Tyrs[i]),k]
							}
						}
					}
				}
			}
		}
	}
	if(all(is.null(LBSPR))){
		pred2 <- NULL
		LF2_new <- NULL
	}

p <- qplot(data = LF_df,
      x = Length,
      y = ..count../sum(..count..),
      color = Fleet,
      fill = Fleet,
      geom = "histogram",
      binwidth = bw) + 
	facet_wrap(~Year) +
	xlab("Length bin (cm)") + ylab("Proportion")
if(nf==1) p <- p + guides(color=FALSE, fill=FALSE)
p

return(p)
}