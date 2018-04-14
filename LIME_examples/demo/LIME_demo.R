rm(list=ls())

## Packages
devtools::install_github("merrillrudd/LIME", dependencies=TRUE, ref="multifleet")
library(LIME)

devtools::install_github("kaskr/TMB_contrib_R/TMBhelper", dependencies=TRUE)
library(TMBhelper)

library(ggplot2)
library(dplyr)

##----------------------------------------------------------------
## Step 1: Read in length data
##----------------------------------------------------------------

setwd("C:\\merrill\\LIME\\LIME_examples\\demo")
data <- read.csv("Ex3_LBSPR_dat.csv", header=TRUE)

## identify length bins
bins <- data[,1]

## setup length frequency matrix
lf <- matrix(data[,2], nrow=1, ncol=nrow(data))
rownames(lf) <- "2016"
colnames(lf) <- bins

## plot length composition
plot_LCfits(LFlist=list("LF"=lf), ylim=c(0,0.2))

##----------------------------------------------------------------
## Step 2: Specify biological inputs and parameter starting values
##----------------------------------------------------------------
## single fleet
lh <- create_lh_list(vbk=0.117, 
					 linf=38, 
					 t0=-0.01,
					 lwa=3.4e-5, 
					 lwb=2.87, 
					 M50=32,
					 M95=34,
					 maturity_input="length",
					 M=0.119, 
					 h=0.65,
					 S50=20, ## starting value
					 S95=27, ## starting value
					 selex_input="length",
					 selex_type=c("logistic"),
					 CVlen=0.1,
					 SigmaR=0.5,
					 SigmaF=0.1,
					 binwidth=2,
					 nfleets=1)

ggplot(lh$df %>% dplyr::filter(By=="Age")) + 
geom_line(aes(x=X, y=Value, color=Fleet), lwd=2) +
facet_wrap(~Variable, scale="free_y") +
xlab("Age")

##----------------------------------------------------------------
## Step 3: Run LIME
##----------------------------------------------------------------

## Single year only
data_list <- list("years"=2016, "LF"=lf)

input_data <- create_inputs(lh=lh, input_data=data_list)

res <- run_LIME(modpath=NULL, input=input_data, data_avail="LC")

## check TMB inputs
Inputs <- res$Inputs

## Report file
Report <- res$Report

## Standard error report
Sdreport <- res$Sdreport

## check convergence
hessian <- Sdreport$pdHess
gradient <- res$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE


## SPR
spr <- Report$SPR
spr

## standard error
sd_spr <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),2]

## lower confidence limit
max(0,spr - 1.96 * sd_spr)

## upper confidence limit
spr + 1.96 * sd_spr

## plot length composition data and fits
plot_LCfits(Inputs=Inputs, 
			Report=Report,
			ylim=c(0,0.3),
			true_years=input_data$years)		

## plot model output
plot_output(Inputs=Inputs, 
			Report=Report,
			Sdreport=Sdreport, 
			lh=lh,
			plot=c("Selex"))
abline(v=lh$linf/bw, lwd=2, lty=2, col="red")
abline(v=sum((Inputs$Data$LF_tlf[,,1] * Inputs$Data$lbmids))/sum(Inputs$Data$LF_tlf[,,1]) / bw, lwd=2, lty=2, col="blue")
legend("topleft", legend=c("Selectivity", "Mean length in catch", "Linf"), col=c("#00AA00", "blue", "red"), lty=c(1,2,2), lwd=2)

##----------------------------------------------------------------
## Step 4: Compare with LB-SPR
##----------------------------------------------------------------
library(LBSPR)

LB_pars <- new("LB_pars")
LB_pars@MK <- 1
LB_pars@Linf <- lh$linf
LB_pars@L50 <- lh$ML50 
LB_pars@L95 <- lh$ML95
LB_pars@Walpha <- lh$lwa
LB_pars@Wbeta <- lh$lwb
LB_pars@BinWidth <- lh$binwidth
LB_pars@Steepness <- lh$h
LB_pars@R0 <- 1
LB_pars@L_units<-"cm"
LB_lengths <- new("LB_lengths")
LB_lengths@LMids <- lh$mids
LB_lengths@LData <- t(data_list$LF)
LB_lengths@Years <- as.numeric(rownames(data_list$LF))
LB_lengths@NYears <- as.numeric(length(rownames(data_list$LF)))


lbspr_res <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths, yrs=NA, Control=list(modtype="GTG"))

spr2 <- lbspr_res@SPR
spr2

sd_spr2 <- sqrt(lbspr_res@Vars[1,"SPR"])

## lower confidence limit
sapply(1:length(spr2), function(x) max(0,spr2[x] - 1.96 * sd_spr2))

## upper confidence limit
spr2 + 1.96 * sd_spr2

## plot length composition data and fits
plot_LCfits(Inputs=Inputs, 
			Report=Report,
			ylim=c(0,0.3),
			LBSPR=lbspr_res,
			true_years=input_data$years)		

## plot model output
plot_output(Inputs=Inputs, 
			Report=Report,
			Sdreport=Sdreport, 
			LBSPR=lbspr_res,
			lh=lh,
			plot=c("Selex"))
abline(v=lh$linf/bw, lwd=2, lty=2, col="red")
abline(v=sum((Inputs$Data$LF_tlf[,,1] * Inputs$Data$lbmids))/sum(Inputs$Data$LF_tlf[,,1]) / bw, lwd=2, lty=2, col="blue")
legend("topleft", legend=c("LIME Selectivity", "LBSPR Selectivity", "Mean length in catch", "Linf"), col=c("#00AA00", "#AA00AA", "blue", "red"), lty=c(1,1,2,2), lwd=2)

##----------------------------------------------------------------
## Explore multifleet
##----------------------------------------------------------------

	##----------------------------------------------------------------
	## Step 1: Read in length data
	##----------------------------------------------------------------

data2 <- read.csv("Data_multifleet.csv", header=TRUE)

## identify length bins
bins <- data2[,1]

## setup length frequency matrix
LF <- list()
LF[[1]] <- matrix(data2[,2], nrow=1, ncol=length(bins))
LF[[2]] <- matrix(data2[,3], nrow=1, ncol=length(bins))
rownames(LF[[1]]) <- rownames(LF[[2]]) <- "2016"
colnames(LF[[1]]) <- colnames(LF[[2]]) <- bins

## plot length composition
plot_LCfits(LFlist=LF, ylim=c(0,0.25))

##----------------------------------------------------------------
## Step 2: Specify biological inputs and parameter starting values
##----------------------------------------------------------------
## single fleet
lh2 <- create_lh_list(vbk=0.117, 
					 linf=38, 
					 t0=-0.01,
					 lwa=3.4e-5, 
					 lwb=2.87, 
					 M50=32,
					 M95=34,
					 maturity_input="length",
					 M=0.119, 
					 h=0.65,
					 S50=c(20,25),
					 S95=c(27,30), 
					 selex_input="length",
					 selex_type=c("logistic","logistic"),
					 CVlen=0.1,
					 SigmaR=0.5,
					 SigmaF=0.1,
					 binwidth=2,
					 nfleets=2)

ggplot(lh2$df %>% dplyr::filter(By=="Age")) + 
geom_line(aes(x=X, y=Value, color=Fleet), lwd=2) +
facet_wrap(~Variable, scale="free_y") +
xlab("Age")

	##----------------------------------------------------------------
	## Step 3: Run LIME
	##----------------------------------------------------------------

data_list <- list("years"=2016, "LF"=LF)

input_data <- create_inputs(lh=lh2, input_data=data_list)

res <- run_LIME(modpath=NULL, input=input_data, data_avail="LC", est_totalF=TRUE, LFdist=0)

## check TMB inputs
Inputs <- res$Inputs

## Report file
Report <- res$Report

## Standard error report
Sdreport <- res$Sdreport

## check convergence
hessian <- Sdreport$pdHess
gradient <- res$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE

res2 <- get_converged(res)

## check TMB inputs
Inputs <- res2$Inputs

## Report file
Report <- res2$Report

## Standard error report
Sdreport <- res2$Sdreport

## check convergence
hessian <- Sdreport$pdHess
gradient <- res2$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE


## SPR
Report$SPR

## plot length composition data and fits
plot_LCfits(Inputs=Inputs, 
			Report=Report,
			ylim=c(0,0.3),
			true_years=input_data$years)		

## plot model output
plot_output(Inputs=Inputs, 
			Report=Report,
			Sdreport=Sdreport, 
			lh=lh2,
			LBSPR=lbspr_res,
			plot=c("Selex"))
abline(v=lh$linf/bw, lwd=2, lty=2, col="red")
abline(v=sum((Inputs$Data$LF_tlf[,,1] * Inputs$Data$lbmids))/sum(Inputs$Data$LF_tlf[,,1]) / bw, lwd=2, lty=2, col="blue")
legend("topleft", legend=c("LIME Fleet 1", "LIME Fleet 2", "LBSPR Selectivity", "Mean length in catch", "Linf"), col=c("red", "blue", "#AA00AA", "blue", "red"), lty=c(1,1,1,2,2), lwd=2)
