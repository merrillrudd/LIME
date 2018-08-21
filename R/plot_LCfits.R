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
		LF_prop <- array(NA, dim=dim(LF_array))
		for(a in 1:dim(LF_prop)[3]){
			sub <- LF_array[,,a]
			for(y in 1:nrow(sub)){
				LF_prop[y,,a] <- sub[y,]/sum(sub[y,])
			}			
		}
		rownames(LF_prop) <- rownames(LF_array)
		colnames(LF_prop) <- colnames(LF_array)
		LF_df <- reshape2::melt(LF_prop)
		names(LF_df) <- c("Year", "Length", "Fleet", "Proportion")
			LF_df <- LF_df %>%
						dplyr::group_by(Year, Length, Fleet, Proportion)
			LF_df$Year <- factor(LF_df$Year)
			LF_df$Length <- as.numeric(LF_df$Length)
			LF_df$Fleet <- factor(LF_df$Fleet)
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
				sub <- Report$plb[,,x]
				rownames(sub) <- Tyrs
				colnames(sub) <- Inputs$Data$lbmids
				return(sub)
			})

			pred_df <- reshape2::melt(pred)
			names(pred_df) <- c("Year", "Length", "Proportion", "Fleet")
			pred_df <- pred_df %>%
						dplyr::group_by(Year, Length, Proportion, Fleet)
			pred_df$Year <- factor(pred_df$Year)
			pred_df$Length <- as.numeric(pred_df$Length)
			pred_df$Proportion <- as.numeric(pred_df$Proportion)
			pred_df$Fleet <- factor(pred_df$Fleet)
	}
	if(all(is.null(Report))){
		pred <- NULL
	}
	if(all(is.null(LBSPR))==FALSE){
		pred2 <- t(LBSPR@pLCatch)
		rownames(pred2) <- years
		colnames(pred2) <- LBSPR@LMids

			pred_df_mod2 <- reshape2::melt(pred2)
			names(pred_df_mod2) <- c("Year", "Length", "Proportion")
			pred_df_mod2 <- pred_df_mod2 %>%
						dplyr::group_by(Year, Length, Proportion)
			pred_df_mod2$Year <- factor(pred_df_mod2$Year)
			pred_df_mod2$Length <- as.numeric(pred_df_mod2$Length)
			pred_df_mod2$Proportion <- as.numeric(pred_df_mod2$Proportion)
	}
	if(all(is.null(LBSPR))){
		pred2 <- NULL
	}

if(all(is.null(Inputs))){
p <- ggplot(LF_df) + 
	geom_histogram(aes(x=Length, y=..count../sum(..count..), color=Fleet, fill=Fleet), binwidth=binwidth, alpha=0.6) +
	scale_fill_brewer(palette="Set1") +
	facet_wrap(Year~., ncol=5, dir="v") +
	ylab("Proportion") + xlab("Length bin (cm)")
}

if(all(is.null(Report))==FALSE){

	LF_df2 <- LF_df %>% mutate("Type"="Observed") %>% mutate("Model"="Data")
	pred_df2 <- pred_df %>% mutate("Type"="Predicted") %>% mutate("Model"="LIME")
	df_all <- rbind(LF_df2, pred_df2)

	if(all(is.null(LBSPR))==FALSE){
		pred_df2_mod2 <- pred_df_mod2 %>% mutate("Type"="Predicted") %>% mutate("Model"="LBSPR") %>% mutate("Fleet"=factor(1))
		df_all <- dplyr::bind_rows(df_all, pred_df2_mod2)
	}

	p <- ggplot(df_all) + 
		geom_ribbon(data=df_all %>% filter(Type=="Observed"), aes(x=Length, ymin=0, ymax=Proportion, fill=Fleet), alpha=0.6) +
		geom_line(data=df_all %>% filter(Type=="Predicted"), aes(x=Length, y=Proportion, color=Model), lwd=1.2) +
		scale_fill_brewer(palette="Set1") +
		scale_color_brewer(palette="Set1", direction=-1) +
		facet_wrap(Year~., ncol=5, dir="v")  +
		xlab("Length bin (cm)") + ylab("Proportion")

}
if(nf==1) p <- p + guides(fill=FALSE)
p

return(p)
}