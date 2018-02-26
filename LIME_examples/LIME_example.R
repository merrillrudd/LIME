rm(list=ls())

## Packages

devtools::install_github("merrillrudd/LIME", dependencies=TRUE, ref="master")
library(LIME)

devtools::install_github("kaskr/TMB_contrib_R/TMBhelper", dependencies=TRUE)
library(TMBhelper)

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

##----------------------------------------------------------------
## Step 1: Specify biological inputs and parameter starting values
##----------------------------------------------------------------

lh <- create_lh_list(vbk=0.21, 
					 linf=65, 
					 t0=-0.01,
					 lwa=0.0245, 
					 lwb=2.79, 
					 S50=c(20), 
					 S95=c(26), 
					 selex_input="length",
					 selex_type=c("logistic"),
					 M50=34,
					 M95=NULL,
					 maturity_input="length",
					 M=0.27, 
					 binwidth=2,
					 CVlen=0.1,
					 SigmaR=0.2,
					 SigmaF=0.2,
					 SigmaC=0.2,
					 SigmaI=0.2,
					 R0=1,
					 qcoef=1e-5,
					 start_ages=0,
					 rho=0.43,
					 nseasons=1)

par(mfrow=c(2,2))
plot(lh$L_a, type="l", lwd=3, col="forestgreen")
plot(lh$W_a, type="l", lwd=3, col="forestgreen")
plot(lh$Mat_a, type="l", lwd=3, col="forestgreen")
plot(lh$S_l, type="l", lwd=3, col="forestgreen")

##----------------------------------------------------
## Step 2: Setup data input
## ---------------------------------------------------
## Demonstrate data generation option
## specify model path to save true population/generated data
true <- generate_data(modpath=NULL,
					  itervec=1, 
					  Fdynamics="Endogenous",
					  Rdynamics="Constant",
					  lh=lh,
					  Nyears=20,
					  Nyears_comp=20,
					  comp_sample=200,
					  init_depl=0.7,
					  seed=28)


## plot simulated data
par(mfrow=c(2,2))
plot(true$F_t, type="l", lwd=3, col="steelblue", ylim=c(0, max(true$F_t)*1.2))
plot(true$R_t, type="l", lwd=3, col="steelblue", ylim=c(0, max(true$R_t)*1.2))
plot(true$SPR_t, type="l", lwd=3, col="steelblue", ylim=c(0, 1))
plot(true$D_t, type="l", lwd=3, col="steelblue", ylim=c(0, max(true$D_t)*1.2))


#######################################
## Length comp data input options
#######################################
LF_matrix <- true$LF

## example with length data only
data_LF <- list("years"=1:true$Nyears, "LF"=LF_matrix)

## plot generated length data
plot_LCfits(LFlist=list("LF"=LF_matrix))


##----------------------------------------------------
## Step 3: Run Model
## ---------------------------------------------------

##--------------------------
## build inputs and object
##--------------------------

res <- run_LIME(modpath=NULL,
				lh=lh,  
				input=data_LF,
				data_avail="LC",
				LFdist=1,
				C_opt=0,
				est_sigma="log_sigma_R",
				fix_param=FALSE,
				fix_param_t=FALSE,
				randomR=TRUE,
				newtonsteps=3,
				F_up=10,
				S50_up=lh$linf,
				derive_quants=FALSE,
				itervec=NULL,
				rewrite=TRUE,
				simulation=FALSE)


## run LIME - may take a few minutes
## looking for outer mgc to minimize and ustep moving towards 1 for well-behaved model
## specify model path to save results
## length comp only
# start <- Sys.time()
# res <- run_LIME(modpath=NULL,
# 				lh=lh,
# 				input_data=data_LF,
# 				est_sigma="log_sigma_R", 
# 				data_avail="LC")
# end <- Sys.time() - start

## check convergence
check <- res$df

## check for other issues
issues <- res$issues

## check TMB inputs
Inputs <- res$Inputs

## Report file
Report <- res$Report

## Standard error report
Sdreport <- res$Sdreport

##----------------------------------------------------
## Step 4: Plot fits
## ---------------------------------------------------
## plot length composition data
plot_LCfits(LFlist=list("LF"=LF_matrix), 
			Inputs=Inputs, 
			Report=Report)

## plot model output
plot_output(Inputs=Inputs, 
			Report=Report,
			Sdreport=Sdreport, 
			lh=lh,
			True=true, 
			plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
			set_ylim=list("Fish" = c(0,1), "SPR" = c(0,1)))




##----------------------------------------------------------------------------
## Step 5: Run with other data types - catch and abundance index examples
## ---------------------------------------------------------------------------
Cw_t <- true$Cw_t ## (optional) catch data, with elements of vector named with the year observed, e.g. 1:20 or 1998:2017 (can be any years within years to model)
I_t <- true$I_t ## (optional) abundance index, with elements of vector named with the year observed, e.g. 1:20 or 1998:2017 (can be any years within years to model)

data_LF_Catch <- list("years"=years, "LF"=LF, "C_t"=Cw_t) ## length comp + catch
data_LF_Index <- list("years"=years, "LF"=LF, "I_t"=I_t) ## length comp + index
data_rich <- list("years"=years, "LF"=LF, "I_t"=I_t, "C_t"=Cw_t) ## length comp + index


## length comp + index + catch
## add option to specify catch in biomass (2) or numbers (1)
res <- run_LIME(modpath=NULL,
				lh=lh,
				input_data=data_rich,
				est_sigma="log_sigma_R", 
				data_avail="Index_Catch_LC",
				C_opt=2) 

## check convergence
check <- res$df
## check for other issues
issues <- res$issues
## check TMB inputs
Inputs <- res$Inputs
## Report file
Report <- res$Report
## Standard error report
Sdreport <- res$Sdreport

## plot length composition data
plot_LCfits(LFlist=list("LF"=LF), 
			Inputs=Inputs, 
			Report=Report,
			dim=c(5,2))


## plot model output
plot_output(Inputs=Inputs, 
			Report=Report,
			Sdreport=Sdreport, 
			lh=lh,
			True=true, 
			plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
			set_ylim=list("Fish" = c(0,1), "SPR" = c(0,1)))



## length comp + catch
## add option to specify catch in biomass (2) or numbers (1)
res <- run_LIME(modpath=NULL,
				lh=lh,
				input_data=data_LF_Catch,
				est_sigma="log_sigma_R", 
				data_avail="Catch_LC",
				C_opt=2)

## check convergence
check <- res$df
## check for other issues
issues <- res$issues
## check TMB inputs
Inputs <- res$Inputs
## Report file
Report <- res$Report
## Standard error report
Sdreport <- res$Sdreport

## plot length composition data
plot_LCfits(LFlist=list("LF"=LF), 
			Inputs=Inputs, 
			Report=Report,
			dim=c(5,2))


## plot model output
plot_output(Inputs=Inputs, 
			Report=Report,
			Sdreport=Sdreport, 
			lh=lh,
			True=true, 
			plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
			set_ylim=list("Fish" = c(0,1), "SPR" = c(0,1)))



## length comp + index
res <- run_LIME(modpath=NULL,
				lh=lh,
				input_data=data_LF_Index,
				est_sigma="log_sigma_R", 
				data_avail="Index_LC")

## check convergence
check <- res$df
## check for other issues
issues <- res$issues
## check TMB inputs
Inputs <- res$Inputs
## Report file
Report <- res$Report
## Standard error report
Sdreport <- res$Sdreport

## plot length composition data
plot_LCfits(LFlist=list("LF"=LF), 
			Inputs=Inputs, 
			Report=Report,
			dim=c(5,2))


## plot model output
plot_output(Inputs=Inputs, 
			Report=Report,
			Sdreport=Sdreport, 
			lh=lh,
			True=true, 
			plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
			set_ylim=list("Fish" = c(0,1), "SPR" = c(0,1)))
