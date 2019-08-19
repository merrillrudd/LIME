rm(list=ls())

#####################
## Load packages
#####################
## install LIME from github
## installation instructions available at https://github.com/merrillrudd/LIME/wiki/2---Introduction-and-installation
## (will need to download Rtools and some R package dependencies available on CRAN)
devtools::install_github("merrillrudd/LIME")
library(LIME)
library(ggplot2)
library(dplyr)

LIME_dir <- "C:\\merrill\\LIME\\examples"
demo_dir <- file.path(LIME_dir, "demo_pulldata")

##----------------------------------------------------------------
## Step 1: Specify biological inputs and parameter starting values
##----------------------------------------------------------------
## single fleet
lh <- create_lh_list(vbk=0.21, ## von Bertalanffy growth parameters (required)
					 linf=65, 
					 t0=-0.01,
					 lwa=0.0245, ## length-weight parameters (required)
					 lwb=2.79, 
					 S50=20, ## logistic selectivity-at-length parameters (starting values)
					 S95=26, 
					 selex_input="length",
					 selex_type=c("logistic"), 
					 M50=34,	## logistic maturity-at-length parameters (required)
					 M95=39,
					 maturity_input="length",
					 M=0.38,  ## natural mortality (required)
					 SigmaR=0.6, ## recruitment standard deviation (starting value)
					 SigmaF=0.2, ## fishing mortality standard deviation penalty (assumed, fixed)
					 SigmaC=0.1, ## catch observation error (assumed, fixed)
					 SigmaI=0.1, ## index observation error (assumed, fixed)
					 R0=1,  ## equilibrium recruitment (starting value only when catch data included)
					 qcoef=1e-5, ## catchability coefficient (starting value only when index data included)
					 binwidth=1, ## length data bin width
					 nseasons=1, ## number of seasons per year (must have data at this level)
					 nfleets=1) ## number of fleets (when data available for multiple fleets)

## plot life history by age or length
par(mfrow=c(2,2), mar=c(4,4,3,1))
plot(lh$L_a, type="l", lwd=4, xlab="Age", ylab="Length (cm)")
plot(lh$W_a, type="l", lwd=4, xlab="Age", ylab="Weight (g)")
plot(lh$Mat_l, type="l", lwd=4, xlab="Length (cm)", ylab="Proportion mature")
## compare with selectivity
plot(lh$S_fl[1,], type="l", lwd=4, xlab="Length (cm)", ylab="Proportion selected to gear")

## plot probability of being an length given age
plba <- with(lh, age_length(highs, lows, L_a, CVlen))
ramp <- colorRamp(c("purple4", "darkorange"))
col_vec <- rgb( ramp(seq(0, 1, length = nrow(plba))), max = 255)
par(mfrow=c(1,1))
matplot(t(plba[-1,]), type="l", lty=1, lwd=3, col=col_vec, xaxs="i", yaxs="i", ylim=c(0, 0.5), xlab="Length bin (cm)", ylab="Density")
legend("topright", legend=lh$ages[seq(2,length(lh$ages),by=3)], col=col_vec[seq(2,length(lh$ages),by=3)],, lwd=3, title="Age")

##----------------------------------------------------
## Step 2: Setup data input
## ---------------------------------------------------

data_dir <- file.path(demo_dir, "sim_data")

length_data <- read.csv(file.path(data_dir, "LF_20yrs.csv"), row=1)
	## rename column names with length bins
	colnames(length_data) <- seq(from=lh$binwidth/2, by=lh$binwidth, length.out=ncol(length_data))

	## years are row names
	years <- as.numeric(rownames(length_data))

	## make sure length data cast as matrix
	length_matrix <- as.matrix(length_data, nrow=length(years), ncol=ncol(length_data))
	rownames(length_matrix) <- factor(years)

catch_data <- read.csv(file.path(data_dir, "Catch_20yrs.csv"), row=1, header=FALSE)
	## cast as vector with years as named elements
	catch_data <- t(catch_data)
	colnames(catch_data) <- years

	## row is fleet
	rownames(catch_data) <- 1

	## make sure catch data cast as matrix
	catch_matrix <- as.matrix(catch_data, nrow=nrow(catch_data), ncol=ncol(catch_data))

index_data <- read.csv(file.path(data_dir, "Index_20yrs.csv"), row=1, header=FALSE)
	index_data <- t(index_data)
	colnames(index_data) <- years

	## row is fleet
	rownames(index_data) <- 1

	## make sure catch data cast as matrix
	index_matrix <- as.matrix(index_data, nrow=nrow(index_data), ncol=ncol(index_data))

####################
## plot length data
####################
## turn length frequency data to a data-frame
LF_df <- LFreq_df(length_matrix)

## LIME function using argument 'LF_df' to fit length data
plot_LCfits(LF_df=LF_df)

###############################
## plot catch and effort data
###############################

par(mfrow=c(2,1))
plot(years, catch_matrix, type="l", lwd=2, ylim=c(0, max(catch_matrix)*1.2), xlab="Year", ylab="Catch (biomass)")
plot(years, index_matrix, type="l", lwd=2, ylim=c(0, max(index_matrix)*1.2), xlab="Year", ylab="Abundance index")


#######################################
## Create data input list
#######################################
## example with length data only
data_LF <- list("years"=years, "LF"=LF_df)

## create model inputs with life history information and data
## outputs length data as array
inputs_LC <- create_inputs(lh=lh, input_data=data_LF)

#######################################
## Other data type input options
#######################################
data_all <- list("years"=years, "LF"=LF_df, "I_ft"=index_matrix, "C_ft"=catch_matrix)
inputs_all <- create_inputs(lh=lh, input_data=data_all)

##----------------------------------------------------
## Step 3: Run models
## ---------------------------------------------------
#######################################
## length-data only
#######################################

## LBSPR
library(LBSPR)
LB_pars <- new("LB_pars")
LB_pars@MK <- inputs_all$M/inputs_all$vbk
LB_pars@Linf <- inputs_all$linf
LB_pars@L50 <- inputs_all$ML50
LB_pars@L95 <- inputs_all$ML95
LB_pars@Walpha <- inputs_all$lwa
LB_pars@Wbeta <- inputs_all$lwb
LB_pars@R0 <- inputs_all$R0
LB_pars@Steepness <- ifelse(inputs_all$h==1, 0.99, inputs_all$h)
LB_pars@BinWidth <- inputs_all$binwidth

LB_lengths <- new("LB_lengths")
LB_lengths@LMids <- inputs_all$mids
LB_lengths@LData <- t(matrix(inputs_all$LF, ncol=length(inputs_all$mids)))
LB_lengths@Years <- as.numeric(rownames(inputs_all$LF))
LB_lengths@NYears <- ncol(LB_lengths@LData)

lbspr <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths)

## LIME
# inputs_LC$SigmaR <- 0.0001
lc_only <- run_LIME(modpath=NULL, 
				input=inputs_LC,
				data_avail="LC",
				est_rdev_t=c(rep(0,10),rep(1,10)))

## check TMB inputs
Inputs <- lc_only$Inputs

## Report file
Report <- lc_only$Report

## Standard error report
Sdreport <- lc_only$Sdreport

## check convergence
hessian <- Sdreport$pdHess
gradient <- lc_only$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE

##----------------------------------------------------
## Step 4: Plot results
## ---------------------------------------------------
## plot length composition data
plot_LCfits(Inputs=Inputs, 
			Report=Report,
			LBSPR=lbspr)		

## plot model output
plot_output(Inputs=Inputs, 
			Report=Report,
			Sdreport=Sdreport, 
			lh=lh,
			LBSPR=lbspr,
			plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
			set_ylim=list("Fish" =c(0,1), "SPR" = c(0,1), "SB"=c(0,2)))	



#######################################
## Data-rich test
#######################################
## include abundance index, catch, and length composition data with LIME
rich <- run_LIME(modpath=NULL, 
				input=inputs_all,
				data_avail="Index_Catch_LC",
				C_type=2) 

## check TMB inputs
Inputs <- rich$Inputs

## Report file
Report <- rich$Report

## Standard error report
Sdreport <- rich$Sdreport

## check convergence
hessian <- Sdreport$pdHess
gradient <- rich$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE

## plot length composition data and fits
plot_LCfits(Inputs=Inputs, 
			Report=Report,
			LBSPR=lbspr)	

## plot model output
plot_output(Inputs=Inputs, 
			Report=Report,
			Sdreport=Sdreport, 
			lh=lh,
			LBSPR=lbspr,
			plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
			set_ylim=list("SPR" = c(0,1)))


#######################################
## Catch + length data
#######################################
## include catch and length data with LIME
catch_lc <- run_LIME(modpath=NULL, 
				input=inputs_all,
				data_avail="Catch_LC",
				C_type=2)

## check TMB inputs
Inputs <- catch_lc$Inputs

## Report file
Report <- catch_lc$Report

## Standard error report
Sdreport <- catch_lc$Sdreport

## check convergence
hessian <- Sdreport$pdHess
gradient <- catch_lc$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE

## plot length composition data
plot_LCfits(Inputs=Inputs, 
			Report=Report)		

## plot model output
plot_output(Inputs=Inputs, 
			Report=Report,
			Sdreport=Sdreport, 
			lh=lh,
			plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
			set_ylim=list("SPR" = c(0,1), "SB"=c(0,2)))		

#######################################
## Index + length data
#######################################
index_lc <- run_LIME(modpath=NULL, 
				input=inputs_all,
				data_avail="Index_LC")

## check TMB inputs
Inputs <- index_lc$Inputs

## Report file
Report <- index_lc$Report

## Standard error report
Sdreport <- index_lc$Sdreport

## check convergence
hessian <- Sdreport$pdHess
gradient <- index_lc$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE


## plot length composition data
plot_LCfits(Inputs=Inputs, 
			Report=Report)		

## plot model output
plot_output(Inputs=Inputs, 
			Report=Report,
			Sdreport=Sdreport, 
			lh=lh,
			plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
			set_ylim=list("SPR" = c(0,1), "SB"=c(0,2)))		


##----------------------------------------------------
## Shorter time series and data gaps
## ---------------------------------------------------

## still modeling over 20 years but with data gaps
## all data types must still be 20 years long
## 		but with 0s for years with no length composition data
## 		-1 for years with no catch or abundance data
## remove first 15 years of length data
length_5 <- length_matrix
length_5[1:15,] <- 0

## remove first 10 years of index data
index_10 <- index_matrix
index_10[1:10] <- -1

## remove every other year of catch data
catch_alt <- catch_matrix
catch_alt[seq(2,length(catch_matrix),by=2)] <- -1

## plot length data
lf5_df <- LFreq_df(length_5)
plot_LCfits(lf5_df)

## plot catch and effort
par(mfrow=c(2,1))
cplot <- catch_alt
cplot[which(cplot < 0)] <- NA
plot(years, cplot, type="o", pch=19, lwd=2, ylim=c(0, max(catch_matrix)*1.2), xlab="Year", ylab="Catch (biomass)")

iplot <- index_10
iplot[which(iplot < 0)] <- NA
plot(years, iplot, type="o", pch=19, lwd=2, ylim=c(0, max(index_matrix)*1.2), xlab="Year", ylab="Abundance index")

## setup data list
data_gaps <- list("years"=years, "LF"=length_5, "C_ft"=catch_alt, "I_ft"=index_10)

## data inputs
inputs_gaps <- create_inputs(lh=lh, input_data=data_gaps)


#######################################
## Data-rich test
#######################################
rich_gaps <- run_LIME(modpath=NULL, 
				input=inputs_gaps,
				data_avail="Index_Catch_LC",
				C_type=2) 

## check that all parameters are identifiable
check <- TMBhelper::Check_Identifiable(rich_gaps$obj)

## check TMB inputs
Inputs <- rich_gaps$Inputs

## Report file
Report <- rich_gaps$Report

## Standard error report
Sdreport <- rich_gaps$Sdreport

## check convergence
hessian <- Sdreport$pdHess
gradient <- rich_gaps$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE

### compare with LBSPR using last 5 years of length data
## LBSPR
LB_pars <- new("LB_pars")
LB_pars@MK <- inputs_gaps$M/inputs_gaps$vbk
LB_pars@Linf <- inputs_gaps$linf
LB_pars@L50 <- inputs_gaps$ML50
LB_pars@L95 <- inputs_gaps$ML95
LB_pars@Walpha <- inputs_gaps$lwa
LB_pars@Wbeta <- inputs_gaps$lwb
LB_pars@R0 <- inputs_gaps$R0
LB_pars@Steepness <- ifelse(inputs_gaps$h==1, 0.99, inputs_gaps$h)
LB_pars@BinWidth <- inputs_gaps$binwidth

LB_lengths <- new("LB_lengths")
LB_lengths@LMids <- inputs_gaps$mids

LF5_lbspr <- matrix(inputs_gaps$LF[-which(rowSums(inputs_gaps$LF)==0),,1], ncol=length(inputs_gaps$mids))
rownames(LF5_lbspr) <- which(rowSums(inputs_gaps$LF[,,1])>0)
LB_lengths@LData <- t(LF5_lbspr)
LB_lengths@Years <- as.numeric(rownames(LF5_lbspr))
LB_lengths@NYears <- nrow(LF5_lbspr)

lbspr2 <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths)

## plot model output
## points show years with length data
plot_output(Inputs=Inputs, 
			Report=Report,
			Sdreport=Sdreport, 
			lh=lh,
			LBSPR=lbspr2,
			plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
			set_ylim=list("SPR" = c(0,1), "SB"=c(0,2)))		
