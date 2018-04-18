rm(list=ls())

## Packages
devtools::install_github("merrillrudd/LIME", dependencies=TRUE, ref="multifleet")
library(LIME)

devtools::install_github("kaskr/TMB_contrib_R/TMBhelper", dependencies=TRUE)
library(TMBhelper)

library(ggplot2)
library(dplyr)

##----------------------------------------------------------------
## *** LIME, single year, single fleet
##----------------------------------------------------------------

	##----------------------------------------------------------------
	## Step 1: Read in length data
	##----------------------------------------------------------------

setwd("C:\\merrill\\LIME\\LIME_examples\\demo")
data <- read.csv("Length_singleyear.csv", header=TRUE)

## identify length bins
bins <- data[,1]

## setup length frequency matrix
lf <- matrix(data[,2], nrow=1, ncol=nrow(data))
rownames(lf) <- "2016"
colnames(lf) <- bins


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
					 S50=26, ## starting value
					 S95=34, ## starting value
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

LFlist <- list()
LFlist[[1]] <- matrix(input_data$LF[,,1], nrow=length(input_data$years))
colnames(LFlist[[1]]) <- input_data$highs
rownames(LFlist[[1]]) <- input_data$years

plot_LCfits(LFlist=LFlist, ylim=c(0,0.2), true_years=input_data$years)

res <- run_LIME(modpath=NULL, input=input_data, data_avail="LC", Rdet=TRUE)

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

	##----------------------------------------------------------------
	## Step 4: Examine output
	##----------------------------------------------------------------

## SPR
spr <- Report$SPR_t

## standard error
sd_spr <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),2]

## lower confidence limit
lcl_spr <- max(0,spr - 1.96 * sd_spr)

## upper confidence limit
ucl_spr <- spr + 1.96 * sd_spr

save <- data.frame("model"="LIME", "run"="singleyear", "spr"=spr, "lcl_spr"=lcl_spr, "ucl_spr"=ucl_spr)

F50 <- with(lh, uniroot(f=calc_ref, lower=0, upper=3, ages=ages, Mat_a=Mat_a, W_a=W_a, M=M, ref=0.5, type="SPR"))$root
Report$F_y/F50
(1-spr)/(1-F50)

Report$D_t


## plot length composition data and fits
plot_LCfits(Inputs=Inputs, 
			Report=Report,
			ylim=c(0,0.2),
			true_years=input_data$years)
abline(v=lh$linf, col="red", lwd=2, lty=2)		

## plot model output
plot_output(Inputs=Inputs, 
			Report=Report,
			Sdreport=Sdreport, 
			lh=lh,
			plot="Selex")
abline(v=lh$linf, lwd=2, lty=2, col="red")
abline(v=Report$ML_ft[length(Report$ML_ft)], lwd=2, lty=2, col="blue")
legend("topleft", legend=c("Selectivity", "Mean length in catch", "Linf"), col=c("#00AA00", "blue", "red"), lty=c(1,2,2), lwd=2)



##----------------------------------------------------------------
## *** LB-SPR, single year, single fleet
##----------------------------------------------------------------
library(LBSPR)

	##----------------------------------------------------------------
	## Step 1: Read in length data
	##----------------------------------------------------------------
LB_lengths <- new("LB_lengths")
LB_lengths@LMids <- data[,1]
LB_lengths@LData <- matrix(data[,2], ncol=1)
LB_lengths@Years <- 2016
LB_lengths@NYears <- 1

	##----------------------------------------------------------------
	## Step 2: Specify biological inputs and parameter starting values
	##----------------------------------------------------------------
LB_pars <- new("LB_pars")
LB_pars@MK <- 1
LB_pars@Linf <- lh$linf
LB_pars@L50 <- lh$ML50
LB_pars@L95 <- lh$ML95
LB_pars@CVLinf <- 0.1
LB_pars@FecB <- 3
LB_pars@Mpow <- 0
LB_pars@Walpha <- lh$lwa
LB_pars@Wbeta <- lh$lwb
LB_pars@BinWidth <- lh$binwidth


	##----------------------------------------------------------------
	## Step 3: Run LBSPR
	##----------------------------------------------------------------
lbspr_res <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths, Control=list(modtype=c("GTG")))


	##----------------------------------------------------------------
	## Step 4: Examine output and compare
	##----------------------------------------------------------------
spr <- lbspr_res@SPR

sd_spr <- sqrt(lbspr_res@Vars[1,"SPR"])

## lower confidence limit
lcl_spr <- max(0,spr - 1.96 * sd_spr)

## upper confidence limit
ucl_spr <- spr + 1.96 * sd_spr


save <- rbind.data.frame(save, data.frame("model"="LBSPR", "run"="singleyear", "spr"=spr, "lcl_spr"=lcl_spr, "ucl_spr"=ucl_spr))


ggplot(save) +
geom_point(aes(x=run, y=spr, color=model), cex=5) +
geom_linerange(aes(x=run, ymin=lcl_spr, ymax=ucl_spr, color=model)) +
ylim(c(0,max(ucl_spr))) +
xlab("Run") + ylab("SPR")

(1-spr2)/(1-F50)


## plot length composition data and fits
plot_LCfits(Inputs=Inputs, 
			Report=Report,
			ylim=c(0,0.3),
			LBSPR=lbspr_res,
			true_years=input_data$years)		
abline(v=lh$linf, lwd=2, lty=2, col="red")


## plot model output
plot_output(Inputs=Inputs, 
			Report=Report,
			Sdreport=Sdreport, 
			LBSPR=lbspr_res,
			lh=lh,
			true_years=data_list$years,
			plot=c("Selex"))
abline(v=lh$linf, lwd=2, lty=2, col="red")
abline(v=Report$ML_ft[length(Report$ML_ft)], lwd=2, lty=2, col="blue")
legend("topleft", legend=c("LIME Selectivity", "LBSPR Selectivity", "Mean length in catch", "Linf"), col=c("#00AA00", "#AA00AA", "blue", "red"), lty=c(1,1,2,2), lwd=2)


##---------------------------------------------------------------
## *** LIME, Multiple years, single-fleet
##----------------------------------------------------------------

	##----------------------------------------------------------------
	## Step 1: Weight length data
	##----------------------------------------------------------------
catch <- read.csv("Catch_multifleet.csv", header=TRUE)

data_mymf <- read.csv("Length_multiyear_multifleet.csv", header=TRUE)

years <- unique(data_mymf$Year)

cldata <- full_join(data_mymf, catch)

## observed samples by bin by year
nsamps <- sum(cldata %>% select(grep("X", colnames(cldata))))

## multiply observed number of fish in bins by catch, keep fleet information
lweight1 <- ((cldata %>% select(grep("X", colnames(cldata)))) * cldata$Catch) %>% mutate(Fleet=cldata$Fleet)

## add weighted lengths from each fleet
lweight2 <- ((lweight1 %>% filter(Fleet==1)) + (lweight1 %>% filter(Fleet==2))) %>% select(-Fleet)

## calculate proportions weighted by catch
lweight3 <- matrix(t(sapply(1:length(years), function(x) as.numeric(lweight2[x,]/sum(lweight2[x,])))), nrow=length(years), ncol=length(bins))

LF_new <- lweight3 * nsamps
rownames(LF_new) <- years
colnames(LF_new) <- bins

plot_LCfits(LFlist=list("LF"=LF_new), ylim=c(0,0.3), dim=c(5,3))

	##----------------------------------------------------------------
	## Step 3: Run LIME
	##----------------------------------------------------------------

data_list <- list("years"=years, "LF"=LF_new)

## use lh - single fleet
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

	##----------------------------------------------------------------
	## Step 4: Examine results
	##----------------------------------------------------------------

spr <- Report$SPR_t[length(Report$SPR_t)]

## standard error
sd_spr <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),2][length(Report$SPR_t)]

## lower confidence limit
lcl_spr <- max(0,spr - 1.96 * sd_spr)

## upper confidence limit
ucl_spr <- min(1,spr + 1.96 * sd_spr)


save <- rbind.data.frame(save, data.frame("model"="LIME", "run"="multiyear_singlefleet", "spr"=spr[length(spr)], "lcl_spr"=lcl_spr[length(spr)], "ucl_spr"=ucl_spr[length(spr)]))

ggplot(save) +
geom_point(aes(x=run, y=spr, color=model), cex=5) +
geom_linerange(aes(x=run, ymin=lcl_spr, ymax=ucl_spr, color=model)) +
scale_y_continuous(limits=c(0,1)) +
xlab("Run") + ylab("SPR")

plot_LCfits(Inputs=Inputs, Report=Report, ylim=c(0,0.25), true_years=years)

## plot model output
plot_output(Inputs=Inputs, 
			Report=Report,
			Sdreport=Sdreport, 
			lh=lh,
			true_years=years)
abline(v=lh$linf, lwd=2, lty=2, col="red")
abline(v=Report$ML_ft[length(Report$ML_ft)], lwd=2, lty=2, col="blue")

##---------------------------------------------------------------
## *** LBSPR, Multiple years, single-fleet
##----------------------------------------------------------------
LB_lengths <- new("LB_lengths")
LB_lengths@LMids <- bins
LB_lengths@LData <- t(LF_new)
LB_lengths@Years <- years
LB_lengths@NYears <- length(years)


lbspr_res <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths, Control=list(modtype=c("GTG")))

spr <- lbspr_res@SPR[length(lbspr_res@SPR)]

sd_spr <- sqrt(lbspr_res@Vars[1,"SPR"])

## lower confidence limit
lcl_spr <- max(0,spr - 1.96 * sd_spr)

## upper confidence limit
ucl_spr <- spr + 1.96 * sd_spr

save <- rbind.data.frame(save, data.frame("model"="LBSPR", "run"="multiyear_singlefleet", "spr"=spr, "lcl_spr"=lcl_spr, "ucl_spr"=ucl_spr))

ggplot(save) +
geom_point(aes(x=run, y=spr, color=model), cex=5) +
geom_linerange(aes(x=run, ymin=lcl_spr, ymax=ucl_spr, color=model)) +
scale_y_continuous(limits=c(0,NA)) +
xlab("Run") + ylab("SPR")

plot_LCfits(Inputs=Inputs, Report=Report, LBSPR=lbspr_res, ylim=c(0,0.25))

## plot model output
plot_output(Inputs=Inputs, 
			Report=Report,
			Sdreport=Sdreport, 
			lh=lh,
			LBSPR=lbspr_res,
			true_years=years)
abline(v=lh$linf, lwd=2, lty=2, col="red")
abline(v=Report$ML_ft[length(Report$ML_ft)], lwd=2, lty=2, col="blue")
legend("topleft", legend=c("LIME Fleet 1", "LIME Fleet 2", "LBSPR Selectivity", "Mean length in catch", "Linf"), col=c("red", "blue", "#AA00AA", "blue", "red"), lty=c(1,1,1,2,2), lwd=2)


##----------------------------------------------------------------
## *** LIME, single year, multiple fleets
##----------------------------------------------------------------

	##----------------------------------------------------------------
	## Step 1: Read in length data
	##----------------------------------------------------------------

data_1ymf <- read.csv("Length_multifleet.csv", header=TRUE)

## identify length bins
bins <- data_1ymf[,1]

## setup length frequency matrix
LF <- list()
LF[[1]] <- matrix(data_1ymf[,2], nrow=1, ncol=length(bins))
LF[[2]] <- matrix(data_1ymf[,3], nrow=1, ncol=length(bins))
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

res2 <- run_LIME(modpath=NULL, input=input_data, data_avail="LC", est_totalF=TRUE, LFdist=0, Rdet=TRUE)

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

	##----------------------------------------------------------------
	## Step 4: Examine results
	##----------------------------------------------------------------

## not converged, but let's look at results anyway
spr <- Report$SPR_t

## standard error
sd_spr <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),2]

## lower confidence limit
lcl_spr <- max(0,spr - 1.96 * sd_spr)

## upper confidence limit
ucl_spr <- spr + 1.96 * sd_spr


save <- rbind.data.frame(save, data.frame("model"="LIME", "run"="singleyear_multifleet", "spr"=spr, "lcl_spr"=lcl_spr, "ucl_spr"=ucl_spr))

ggplot(save) +
geom_point(aes(x=run, y=spr, color=model), cex=5) +
geom_linerange(aes(x=run, ymin=lcl_spr, ymax=ucl_spr, color=model)) +
scale_y_continuous(limits=c(0,1)) +
xlab("Run") + ylab("SPR")

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
			true_years=data_list$years,
			plot=c("Selex"))
abline(v=lh$linf, lwd=2, lty=2, col="red")
abline(v=Report$ML_ft[length(Report$ML_ft)], lwd=2, lty=2, col="blue")
legend("topleft", legend=c("LIME Fleet 1", "LIME Fleet 2", "LBSPR Selectivity", "Mean length in catch", "Linf"), col=c("red", "blue", "#AA00AA", "blue", "red"), lty=c(1,1,1,2,2), lwd=2)


##----------------------------------------------------------------
## *** LIME, Multiple years, multiple fleets
##----------------------------------------------------------------

	##----------------------------------------------------------------
	## Step 1: Read in length data
	##----------------------------------------------------------------

LF1 <- as.matrix(data_mymf %>% filter(Fleet==1) %>% select(-c(Fleet, Year)))
LF2 <- as.matrix(data_mymf %>% filter(Fleet==2) %>% select(-c(Fleet, Year)))

LFmymf <- list()
LFmymf[[1]] <- LF1
LFmymf[[2]] <- LF2
for(i in 1:2){
	rownames(LFmymf[[i]]) <- years
	colnames(LFmymf[[i]]) <- bins
}


## plot length composition
plot_LCfits(LFlist=LFmymf, ylim=c(0,0.25), dim=c(5,3))



	##----------------------------------------------------------------
	## Step 3: Run LIME
	##----------------------------------------------------------------
data_list <- list("years"=years, "LF"=LFmymf)

## use lh2 - multifleet inputs
input_data <- create_inputs(lh=lh2, input_data=data_list)

## need to decrease the degree to which F can change between years -- I found this one worked but the estiamtes change if this is decreased to 0.20
input_data$SigmaF <- 0.025

res <- run_LIME(modpath=NULL, input=input_data, data_avail="LC", LFdist=0, est_totalF=TRUE)

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

	##----------------------------------------------------------------
	## Step 4: Examine results
	##----------------------------------------------------------------

spr <- Report$SPR_t[length(Report$SPR_t)]

## standard error
sd_spr <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),2][length(Report$SPR_t)]

## lower confidence limit
lcl_spr <- max(0,spr - 1.96 * sd_spr)

## upper confidence limit
ucl_spr <- min(1,spr + 1.96 * sd_spr)


save <- rbind.data.frame(save, data.frame("model"="LIME", "run"="multiyear_multifleet", "spr"=spr[length(spr)], "lcl_spr"=lcl_spr[length(spr)], "ucl_spr"=ucl_spr[length(spr)]))

ggplot(save) +
geom_point(aes(x=run, y=spr, color=model), cex=5) +
geom_linerange(aes(x=run, ymin=lcl_spr, ymax=ucl_spr, color=model)) +
scale_y_continuous(limits=c(0,NA)) +
xlab("Run") + ylab("SPR")

plot_LCfits(Inputs=Inputs, Report=Report, ylim=c(0,0.25), true_years=years)

## plot model output
plot_output(Inputs=Inputs, 
			Report=Report,
			Sdreport=Sdreport, 
			lh=lh2,
			true_years=years)
abline(v=lh$linf, lwd=2, lty=2, col="red")
abline(v=Report$ML_ft[length(Report$ML_ft)], lwd=2, lty=2, col="blue")




