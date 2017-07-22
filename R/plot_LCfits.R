#' plot length composition data and model fits
#'
#' \code{plot_LCfits} plot length composition data with option for model fits from LIME and LB-SPR
#'
#' @param Inputs LIME input file, must have "LF"=length frequency and "T_yrs"=total vector of years
#' @param Inputs2 optional LIME input file for length data comparison; default NULL
#' @param Inputs3 optional LIME input file for length data comparison; default NULL
#' @param Inputs4 optional LIME input file for length data comparison; default NULL
#' @param Report LIME report file; default NULL if only plotting length data
#' @param Report2 optional LIME report file for model fit comparison; default NULL
#' @param LBSPR LBSPR results - must have pLF = probability of being harvested in a length bin; default NULL
#' @param true_lc_years true years of length comp data for plotting; default NULL
#' @param ylim ylim for plot; default NULL
#' @param ML50 length at 50% maturity for comparison; default NULL
#' @param SL50 length at 50% selectivity for comparison; default NULL
#' @param dim dimensions of plot; default NULL if not specified, will approximate best based on number of years to plot
#' @param n if TRUE, will display sample size of length comp data; default FALSE
#' 
#' @return figure with length composition data and model fits if Report or LBSPR are specified
#' 
#' @export
plot_LCfits <- function(Inputs, Inputs2=NULL, Inputs3=NULL, Inputs4=NULL, Report=NULL, Report2=NULL, LBSPR=NULL, true_lc_years=NULL, ylim=NULL, ML50=NULL, SL50=NULL, dim=NULL, n=FALSE){
	# dev.new()

	obs <- Inputs$LF
	lbhighs <- colnames(obs)
	lc_years <- rownames(obs)

	if(all(is.null(Inputs2))==FALSE){
		obs2 <- Inputs2$LF
		lbhighs2 <- colnames(obs2)
		lc_years2 <- rownames(obs2)
	}
	if(all(is.null(Inputs3))==FALSE){
		obs3 <- Inputs3$LF
		lbhighs3 <- colnames(obs3)
		lc_years3 <- rownames(obs3)
	}
	if(all(is.null(Inputs4))==FALSE){
		obs4 <- Inputs4$LF
		lbhighs4 <- colnames(obs4)
		lc_years4 <- rownames(obs4)
	}


	if(all(is.null(dim))) dim <- c(ceiling(sqrt(length(lc_years))), ceiling(sqrt(length(lc_years))))


	if(all(is.null(Report))==FALSE){
		pred <- Report$plb
		Tyrs <- Inputs$T_yrs
	}
	if(all(is.null(Report))){
		pred <- NULL
		Tyrs <- NULL
	}
	if(all(is.null(LBSPR))==FALSE) pred2 <- t(LBSPR$pLF)
	if(all(is.null(LBSPR))) pred2 <- NULL

	par(mfcol=dim, mar=c(0,0,0,0), omi=c(1,1,1,0.2))

	# find_max <- 0
	# for(i in 1:length(lc_years)){
	# 	yr <- lc_years[i]
	# 		omax <- max(obs[which(lc_years==yr),]/sum(obs[which(lc_years==yr),]))
	# 		if(all(is.null(Inputs2))==FALSE) omax2 <- max(obs2[which(lc_years2==yr),]/sum(obs2[which(lc_years2==yr),]))
	# 		if(all(is.null(Inputs2))) omax2 <- NULL
	# 		if(all(is.null(Report))==FALSE) pmax <- max(pred[which(Tyrs==yr),])
	# 		if(all(is.null(Report))) pmax <- NULL
	# 		max <- max(c(omax,pmax,omax2))
	# 		if(max > find_max) find_max <- max
	# }
	if(all(is.null(ylim))) ylim <- c(0, 0.1)
	for(i in 1:length(lc_years)){
		yr <- lc_years[i]
		barplot(as.numeric(obs[which(lc_years==yr),]/sum(obs[which(lc_years==yr),])), xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=ylim, col="#22222250", border=NA, space=0)
		lines(pred[which(Tyrs==yr),], col="blue", lwd=4)
		lines(pred2[which(lc_years==yr),], col="red", lwd=4)
		box()
		if(all(is.null(Inputs2))==FALSE){
			par(new=TRUE)
			barplot(as.numeric(obs2[which(lc_years2==yr),]/sum(obs2[which(lc_years2==yr),])), border=NA, space=0, col="#DD000050", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=ylim)
		}
		if(all(is.null(Inputs3))==FALSE){
			par(new=TRUE)
			barplot(as.numeric(obs3[which(lc_years2==yr),]/sum(obs3[which(lc_years2==yr),])), border=NA, space=0, col="#0000DD50", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=ylim)
		}
		if(all(is.null(Inputs4))==FALSE){
			par(new=TRUE)
			barplot(as.numeric(obs4[which(lc_years2==yr),]/sum(obs4[which(lc_years2==yr),])), border=NA, space=0, col="#00DD0050", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=ylim)
		}

		xlabs <- pretty(seq_along(lbhighs))
		plot_labs <- rep(NA, length(xlabs))
		if(xlabs[1]!=0) warning("Should start length bin labels at 0")
		plot_labs[1] <- 0
		elabs <- as.numeric(lbhighs[xlabs][which(is.na(lbhighs[xlabs])==FALSE)])
		plot_labs[2:(length(elabs)+1)] <- elabs
		if(is.na(plot_labs[length(plot_labs)])) plot_labs[length(plot_labs)] <- max(elabs) + elabs[1]

		if(i %% dim[1] == 0 | i==length(lc_years)) axis(1, at=xlabs, labels=plot_labs, cex.axis=2)
		if(i %in% 1:dim[1]) axis(2, at=pretty(ylim), las=2, cex.axis=2)
		if(all(is.null(true_lc_years))==FALSE & length(true_lc_years)!=length(lc_years)){
			mtext(side=3, lc_years[i], font=2, cex=2, line=-2)
			warning("Input years for length composition data do not match number of years in analysis")
		}
		if(all(is.null(true_lc_years))) mtext(side=3, line=-3, lc_years[i], font=2, cex=2)
		if(all(is.null(true_lc_years))==FALSE) mtext(side=3, line=-2, true_lc_years[i], font=2)
		if(all(is.null(SL50))==FALSE){
			if(length(SL50)==1) abline(v=SL50, lty=2, lwd=2)
			if(length(SL50)>1) abline(v=SL50[i], lty=2, lwd=2)
		}
		if(all(is.null(ML50))==FALSE) abline(v=ML50, col="orange", lty=2, lwd=2)
		if(n==TRUE) print.letter(xy=c(0.15,0.9), paste0("n = ", sum(Inputs$LF[i,])), cex=2)
	}
	mtext(side=1, "Length bin (cm)", outer=TRUE, line=4, cex=1.5)
	mtext(side=2, "Proportion", outer=TRUE, line=5, cex=1.5)

}