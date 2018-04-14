#' plot length composition data and model fits
#'
#' \code{plot_LCfits} plot length composition data with option for model fits from LIME and LB-SPR
#'
#' @author M.B. Rudd
#' @param LFlist list of length frequency matrices
#' @param Inputs LIME input file; default NULL if only plotting length data
#' @param Report LIME report file; default NULL if only plotting length data
#' @param LBSPR LBSPR results - must have pLF = probability of being harvested in a length bin; default NULL
#' @param ylim ylim for plot; default NULL
#' @param ML50 length at 50 percent maturity for comparison; default NULL
#' @param SL50 length at 50 percent selectivity for comparison; default NULL
#' @param dim dimensions of plot; default NULL if not specified, will approximate best based on number of years to plot
#' @param n if TRUE, will display sample size of length comp data; default FALSE
#' @param true_years optional vector of true years to label years of length comp
#' @importFrom graphics abline axis barplot box legend lines mtext par
#' 
#' @return figure with length composition data and model fits if Report or LBSPR are specified
#' 
#' @export
plot_LCfits <- function(LFlist=NULL, Inputs=NULL, Report=NULL, LBSPR=NULL, ylim=NULL, ML50=NULL, SL50=NULL, dim=NULL, n=FALSE, true_years=NULL,xlim=NULL){
	# dev.new()

	if(all(is.null(Inputs))) LFlist=LFlist
	if(all(is.null(Inputs))==FALSE){
		LF_array <- Inputs$Data$LF_tlf
		LFlist <- list()
		for(i in 1:Inputs$Data$n_f){
			LFlist[[i]] <- matrix(LF_array[,,i], nrow=nrow(LF_array), ncol=ncol(LF_array))
			colnames(LFlist[[i]]) <- colnames(LF_array)
			rownames(LFlist[[i]]) <- rownames(LF_array)
		}
		bins <- as.numeric(colnames(LF_array))
		bw <- bins[1]
	}

	nf <- length(LFlist)
	LCyrs <- lapply(1:nf, function(x){
		rownames(LFlist[[x]])
	})
	lbhighs <- as.numeric(colnames(LFlist[[1]]))
	all_lc_years <- min(as.numeric(unlist(LCyrs))):max(as.numeric(unlist(LCyrs)))
	if(all(is.null(true_years))) true_years <- all_lc_years

	if(all(is.null(dim))) dim <- c(ceiling(sqrt(length(all_lc_years))), ceiling(sqrt(length(all_lc_years))))


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
	}
	if(all(is.null(LBSPR))) pred2 <- NULL


	par(mfcol=dim, mar=c(0,0,0,0), omi=c(1,1,1,0.2))

	if(nf>1){
      colfn <- colorRampPalette(c("red","blue"))
      cols <- colfn(nf)
    }
    if(nf==1) cols <- "#228B22"

	if(all(is.null(ylim))) ylim <- c(0, 0.1)
	if(all(is.null(xlim))) xlim <- c(0,length(lbhighs))
	for(i in seq_along(all_lc_years)){
		yr <- i
		for(f in 1:nf){
			if(f==1){
				flcyrs <- seq_along(as.numeric(LCyrs[[f]]))
				barplot(as.numeric(LFlist[[f]][which(flcyrs==yr),]/sum(LFlist[[f]][which(flcyrs==yr),])), xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=ylim, xlim=xlim, col=paste0(cols[1],"50"), border=NA, space=0)
					lines(pred[[f]][which(Tyrs==yr),], col=cols[1], lwd=4)
					lines(pred2[which(rownames(pred2)==yr),], col="#AA00AA", lwd=4)
				box()
			}
			if(f>1 & all_lc_years[yr] %in% LCyrs[[f]]){
				par(new=TRUE)
				flcyrs <- seq_along(as.numeric(LCyrs[[f]]))
				barplot(as.numeric(LFlist[[f]][which(flcyrs==yr),]/sum(LFlist[[f]][which(flcyrs==yr),])), border=NA, space=0, col=paste0(cols[f],"50"), xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=ylim, xlim=xlim)
				if(sum(LFlist[[f]][which(flcyrs==yr),])>0) lines(pred[[f]][which(Tyrs==yr),], col=cols[f], lwd=4)
			}
		}

			xlabs <- pretty(seq_along(lbhighs))
			plot_labs <- rep(NA, length(xlabs))
			if(xlabs[1]!=0) warning("Should start length bin labels at 0")
			plot_labs[1] <- 0
			elabs <- as.numeric(lbhighs[xlabs][which(is.na(lbhighs[xlabs])==FALSE)])
			plot_labs[2:(length(elabs)+1)] <- elabs
			if(is.na(plot_labs[length(plot_labs)])) plot_labs[length(plot_labs)] <- max(elabs) + elabs[1]	

			if(i %% dim[1] == 0 | i==length(all_lc_years)) axis(1, at=xlabs, labels=plot_labs, cex.axis=2)
			if(i %in% 1:dim[1]) axis(2, at=pretty(ylim), las=2, cex.axis=2)
			if(all(is.null(true_years))==FALSE & length(true_years)!=length(true_years)){
				mtext(side=3, true_years[i], font=2, cex=2, line=-2)
				warning("Input years for length composition data do not match number of years in analysis")
			}
			if(all(is.null(true_years))) mtext(side=3, line=-3, true_years[i], font=2, cex=2)
			if(all(is.null(true_years))==FALSE) mtext(side=3, line=-2, true_years[i], font=2)
			if(all(is.null(SL50))==FALSE){
				if(length(SL50)==1) abline(v=SL50, lty=2, lwd=2)
				if(length(SL50)>1) abline(v=SL50[i], lty=2, lwd=2)
			}
			if(all(is.null(ML50))==FALSE) abline(v=ML50, col="orange", lty=2, lwd=2)
			if(n==TRUE) print_letter(xy=c(0.15,0.9), paste0("n = ", sum(obs[i,,])), cex=2)
		}
	mtext(side=1, "Length bin (cm)", outer=TRUE, line=4, cex=1.5)
	mtext(side=2, "Proportion", outer=TRUE, line=5, cex=1.5)

}