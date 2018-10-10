#' plot length composition data and model fits
#'
#' \code{plot_LCfits} plot length composition data with option for model fits from LIME and LB-SPR
#'
#' @author M.B. Rudd
#' @param LF_df list of length frequency matrices
#' @param binwidth bin width for plotting data with LF_df
#' @param Inputs LIME input file; default NULL if only plotting length data
#' @param Report LIME report file; default NULL if only plotting length data
#' @param LBSPR LBSPR results - must have pLF = probability of being harvested in a length bin; default NULL
#' @param plot_fit default = TRUE if Report file or LBSPR file is included, plot the model fits, otherwise plot data only
#' @param n if TRUE, will display sample size of length comp data; default FALSE
#' @param time_labels default=NULL, otherwise a vector of time series names
#' @importFrom graphics abline axis barplot box legend lines mtext par
#' 
#' @return figure with length composition data and model fits if Report or LBSPR are specified
#' 
#' @export
plot_LCfits <- function(LF_df=NULL, binwidth=1, Inputs=NULL, Report=NULL, LBSPR=NULL, plot_fit=TRUE, n=FALSE, time_labels=NULL){
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
			sub <- matrix(LF_array[,,a], nrow=nrow(LF_array), ncol=ncol(LF_array))
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
	years <- unique(as.numeric(LF_df$Year))[order(unique(as.numeric(LF_df$Year)))] #unique(LF_df$Year)[order(as.numeric(unique(LF_df$Year)))]
	nyears <- length(years)
	# n_yr <- rep(0, nyears)
	# for(i in 1:nyears){
	# 	sub <- LF_df %>% dplyr::filter(Year==years[i])
	# 	n_yr[i] <- nrow(sub)
	# }
	if(all(is.null(time_labels))==FALSE){
		time <- sapply(1:nrow(LF_df), function(x){
			return(time_labels[LF_df$Year[x]])
		})
		LF_df <- data.frame(LF_df)
		LF_df$Month <- sapply(1:nrow(LF_df), function(x) strsplit(time[x],"_")[[1]][1])
		LF_df$Year2 <- sapply(1:nrow(LF_df), function(x) strsplit(time[x],"_")[[1]][2])
	}

	nf <- length(unique(LF_df$Fleet))
	LCyrs <- lapply(1:nf, function(x){
		sub <- LF_df %>% dplyr::filter(Fleet==x)
		yrs <- unique(sub$Year)[order(unique(sub$Year))]
		return(yrs)
	})
	all_lc_years <- min(as.numeric(unlist(LCyrs))):max(as.numeric(unlist(LCyrs)))


	if(all(is.null(Inputs))) Tyrs <- all_lc_years
	if(all(is.null(Inputs))==FALSE) Tyrs <- 1:Inputs$Data$n_t

	if(all(is.null(Report))==FALSE){
		# pred <- Report$plb
			pred <- lapply(1:nf, function(x){
				sub <- matrix(Report$plb[,,x], nrow=length(Tyrs))
				rownames(sub) <- Tyrs
				colnames(sub) <- colnames(Inputs$Data$LF_tlf)
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
		colnames(pred2) <- LBSPR@LMids[1:ncol(pred2)]

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
		geom_histogram(aes(x=Length, y=..density.., color=Fleet, fill=Fleet), position="identity", binwidth=binwidth, alpha=0.6) +
		scale_fill_brewer(palette="Set1") +
		scale_color_brewer(palette="Set1") +
		ylab("Proportion") + xlab("Length bin (cm)")
	if("Month" %in% colnames(LF_df)){
		p <- p + facet_wrap(Year2~factor(Month), dir="v")
	}
	if("Month" %in% colnames(LF_df)==FALSE){
		p <- p + facet_wrap(Year~., dir="v")
	}
	if(nf==1) p <- p + guides(color=FALSE, fill=FALSE)
}

if(all(is.null(Report))==FALSE){

	LF_df2 <- LF_df %>% mutate("Type"="Observed") %>% mutate("Model"="Data")
	pred_df2 <- pred_df %>% mutate("Type"="Predicted") %>% mutate("Model"="LIME")
	df_all <- rbind(LF_df2, pred_df2)

	if(all(is.null(LBSPR))==FALSE){
		pred_df2_mod2 <- pred_df_mod2 %>% mutate("Type"="Predicted") %>% mutate("Model"="LBSPR") %>% mutate("Fleet"=factor(1))
		df_all <- dplyr::bind_rows(df_all, pred_df2_mod2)
	}

  if(length(unique(df_all$Model))>1){
	p <- ggplot(df_all) + 
		geom_ribbon(data=df_all %>% filter(Type=="Observed"), aes(x=Length, ymin=0, ymax=Proportion, fill=Fleet), alpha=0.6) +
		scale_fill_brewer(palette="Set1") +
		xlab("Length bin (cm)") + ylab("Proportion")
	if("Month" %in% colnames(LF_df)){
		p <- p + facet_wrap(Year2~factor(Month), dir="v")
	}
	if("Month" %in% colnames(LF_df)==FALSE){
		p <- p + facet_wrap(Year~., dir="v")
	}
	if(plot_fit==TRUE){
		p <- p + geom_line(data=df_all %>% filter(Type=="Predicted"), aes(x=Length, y=Proportion, color=Model), lwd=1.2) +
				scale_color_brewer(palette="Set1", direction=-1)
	}
  }
  if(length(unique(df_all$Fleet))>1){
  	p <- ggplot(df_all) + 
		geom_ribbon(data=df_all %>% filter(Type=="Observed"), aes(x=Length, ymin=0, ymax=Proportion, fill=Fleet), alpha=0.6) +
		scale_fill_brewer(palette="Set1") +
		facet_wrap(Year~., ncol=5, dir="v")  +
		xlab("Length bin (cm)") + ylab("Proportion")
	if(plot_fit==TRUE){
		p <- p + geom_line(data=df_all %>% filter(Type=="Predicted"), aes(x=Length, y=Proportion, color=Fleet), lwd=1.2) +
				scale_color_brewer(palette="Set1")
	}
  }
if(nf==1) p <- p + guides(fill=FALSE)

}
p

return(p)
}