rm(list=ls())

## Packages

devtools::install_github("merrillrudd/LIME", dependencies=TRUE, build_vignettes=TRUE)
library(LIME)

devtools::install_github("kaskr/TMB_contrib_R/TMBhelper", dependencies=TRUE)
library(TMBhelper)

##----------------------------------------------------------------
## Step 1: Specify biological inputs and parameter starting values
##----------------------------------------------------------------
lh <- create_lh_list(vbk=0.21, 
					 linf=65, 
					 t0=-0.01,
					 lwa=0.0245, 
					 lwb=2.79, 
					 S50=20, 
					 S95=26, 
					 selex_input="length",
					 M50=34,
					 M95=NULL,
					 maturity_input="length",
					 M=0.27, 
					 binwidth=1,
					 CVlen=0.1,
					 SigmaR=0.737,
					 SigmaF=0.2,
					 SigmaC=0.2,
					 SigmaI=0.2,
					 R0=1,
					 qcoef=1e-5,
					 start_ages=0,
					 rho=0.43,
					 nseasons=1)

par(mfrow=c(2,2))
plot(lh$L_a, type="l", lwd=3, col="forestgreen", xlab="Age", ylab="Length")
plot(lh$W_a, type="l", lwd=3, col="forestgreen", xlab="Age", ylab="Weight")
plot(lh$Mat_l, type="l", lwd=3, col="forestgreen", xlab="Length", ylab="Proportion mature")
plot(lh$S_l, type="l", lwd=3, col="forestgreen", xlab="Length", ylab="Proportion vulnerable to gear")

##----------------------------------------------------
## Step 2: Setup data input
## ---------------------------------------------------
## Demonstrate data generation option
## specify model path to save true population/generated data
true <- generate_data(modpath=NULL,
					  itervec=1, 
					  Fdynamics="Ramp",
					  Rdynamics="AR",
					  lh=lh,
					  Nyears=20,
					  Nyears_comp=10,
					  comp_sample=200,
					  init_depl=0.8,
					  seed=143)

## Data input components
years <- true$years ## total years to model, can be 1:20 or 1998:2017
LF <- true$LF ## length composition data, years along rows and length bin along columns. year names should match 'year' quantity, e.g. 11:20 or 2008:2017 (can be any years within years to model)

## input data list
data_LF <- list("years"=years, "LF"=LF) ## length comp only

## plot length composition data
plot_LCfits(Inputs=data_LF) ## "Inputs" argument just must be a list with "LF" as one of the components, e.g. plot_LCfits(Inputs=list("LF"=true$LF))

##----------------------------------------------------
## Step 3: Run Model
## ---------------------------------------------------

## run LIME - may take a few minutes
## looking for outer mgc to minimize and ustep moving towards 1 for well-behaved model
## specify model path to save results
## length comp only
start <- Sys.time()
res <- run_LIME(modpath=NULL,
				lh=lh,
				input_data=data_LF,
				est_sigma="log_sigma_R", 
				data_avail="LC")
end <- Sys.time() - start

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
plot_LCfits(Inputs=Inputs$Data, 
			Report=Report,
			true_lc_years=2008:2017, 
			ylim=c(0,0.15), 
			dim=c(5,2))

## plot model output
plot_output(all_years=1:20,
			lc_years=11:20, 
			Inputs=Inputs, 
			Report=Report, 
			Sdreport=Sdreport, 
			lh=lh, 
			true_years=1998:2017, 
			True=true, 
			plot=c("Fish","Rec","SPR","Selex"),
			set_ylim=list("Fish"=c(0,2)))




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
plot_LCfits(Inputs=Inputs$Data, 
			Report=Report,
			true_lc_years=2008:2017, 
			ylim=c(0,0.15), 
			dim=c(5,2))

## plot model output
plot_output(all_years=1:20,
			lc_years=11:20, 
			Inputs=Inputs, 
			Report=Report, 
			Sdreport=Sdreport, 
			lh=lh, 
			true_years=1998:2017, 
			True=true, 
			plot=c("Fish","Rec","SPR","Selex"),
			set_ylim=list("Fish"=c(0,2)))


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
plot_LCfits(Inputs=Inputs$Data, 
			Report=Report,
			true_lc_years=2008:2017, 
			ylim=NULL, 
			ML50=lh$ML50, 
			SL50=Report$S50,
			dim=c(5,2))

## plot model output
plot_output(all_years=1:20,
			lc_years=11:20, 
			Inputs=Inputs, 
			Report=Report, 
			Sdreport=Sdreport, 
			lh=lh, 
			true_years=1998:2017, 
			True=true, 
			plot=c("Fish","Rec","SPR","Selex"),
			set_ylim=list("Fish"=c(0,2)))


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
plot_LCfits(Inputs=Inputs$Data, 
			Report=Report,
			true_lc_years=2008:2017, 
			ylim=NULL, 
			ML50=lh$ML50, 
			SL50=Report$S50,
			dim=c(5,2))

## plot model output
plot_output(all_years=1:20,
			lc_years=11:20, 
			Inputs=Inputs, 
			Report=Report, 
			Sdreport=Sdreport, 
			lh=lh, 
			true_years=1998:2017, 
			True=true, 
			plot=c("Fish","Rec","SPR","Selex"),
			set_ylim=list("Fish"=c(0,2)))
