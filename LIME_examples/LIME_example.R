rm(list=ls())

## Packages

devtools::install_github("merrillrudd/LIME", dependencies=TRUE, ref="multifleet")
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
# lh <- create_lh_list(vbk=0.21, 
# 					 linf=65, 
# 					 t0=-0.01,
# 					 lwa=0.0245, 
# 					 lwb=2.79, 
# 					 S50=c(20,30), 
# 					 S95=c(26,36), 
# 					 selex_input="length",
# 					 selex_type=c("logistic","logistic"),
# 					 M50=34,
# 					 M95=NULL,
# 					 maturity_input="length",
# 					 M=0.27, 
# 					 binwidth=1,
# 					 CVlen=0.1,
# 					 SigmaR=0.737,
# 					 SigmaF=0.2,
# 					 SigmaC=0.2,
# 					 SigmaI=0.2,
# 					 R0=1,
# 					 qcoef=1e-5,
# 					 start_ages=0,
# 					 rho=0.43,
# 					 nseasons=1,
# 					 nfleets=2)

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
					 binwidth=1,
					 CVlen=0.1,
					 SigmaR=0.2,
					 SigmaF=0.2,
					 SigmaC=0.2,
					 SigmaI=0.2,
					 R0=1,
					 qcoef=1e-5,
					 start_ages=0,
					 rho=0.43,
					 nseasons=1,
					 nfleets=1)

## by age
p <- ggplot(lh$df %>% 
			filter(Variable %in% c("Length","Weight","Selectivity","Maturity")) %>% 
			filter(By == "Age")) + 
	geom_line(aes(x=X, y=Value, color=Fleet), lwd=2) +
	facet_wrap(~Variable, scale="free_y") +
	xlab("Age") +
	theme_lsd()

## by length
p <- ggplot(lh$df %>% 
			filter(Variable %in% c("Selectivity","Maturity")) %>% 
			filter(By == "Length")) + 
	geom_line(aes(x=X, y=Value, color=Fleet), lwd=2) +
	facet_grid(Variable~., scale="free_y") +
	xlab("Length") +
	theme_lsd()


##----------------------------------------------------
## Step 2: Setup data input
## ---------------------------------------------------
## Demonstrate data generation option
## specify model path to save true population/generated data
# true <- generate_data(modpath=NULL,
# 					  itervec=1, 
# 					  Fdynamics=c("Constant","Endogenous"),
# 					  Rdynamics="Constant",
# 					  lh=lh,
# 					  Nyears=20,
# 					  Nyears_comp=c(20,10),
# 					  comp_sample=200,
# 					  init_depl=0.7,
# 					  seed=44,
# 					  fleet_percentage=c(0.5,0.5))
true <- generate_data(modpath=NULL,
					  itervec=1, 
					  Fdynamics=c("Constant"),
					  Rdynamics="Constant",
					  lh=lh,
					  Nyears=20,
					  Nyears_comp=c(20),
					  comp_sample=200,
					  init_depl=0.7,
					  seed=28,
					  fleet_percentage=1)


## plot simulated data
dfsim <- true$dfsim
p <- ggplot(dfsim %>% filter(Variable %in% c("LengthComp", "TotalBiomass", "SpawningBiomass") == FALSE)) +
	geom_line(aes(x=X, y=Value, colour=Fleet), lwd=2) +
	facet_wrap(~Variable, scale='free_y') +
	expand_limits(y=0) +
	xlab("Year") +
	theme_lsd()

#######################################
## Length comp data input options
#######################################
## Option 1: Length comp array
LF_array <- true$LF ## array with rows = years, columns = upper length bins, 3rd dimension = fleets

## Option 2: Length comp list
LF_list <- lapply(1:lh$nfleets, function(x) true$LF[,,x]) ##list with 1 element per fleet, and each element is a matrix with rows = years, columns = upper length bins

	## plot length composition data using LF_list
	plot_LCfits(LFlist=LF_list) ## "Inputs" argument just must be a list with "LF" as one of the components, e.g. plot_LCfits(Inputs=list("LF"=true$LF))

## Option 3: Data frame
LF_df <- true$dfsim %>% filter(Variable == "LengthComp") ## long-form data frame where "X" = year, "Value"=length measurement, and "Fleet"=discrete variables representing a fleet. 

## example with length data only
data_LF <- list("years"=1:true$Nyears, "LF"=LF_array)

##if using multinomial distribution, must specify annual effective sample size by fleet
data_LF_neff <- list("years"=1:true$Nyears, "LF"=LF_array, "neff_ft"=true$obs_per_year)

## create model inputs with life history information and data
## outputs length data as array
inputs_LC <- create_inputs(lh=lh, input_data=data_LF)


##----------------------------------------------------
## Step 3: Run Model
## ---------------------------------------------------
##-------------------------
## template file
# ##--------------------------
src_dir <- file.path("C:\\merrill\\LIME\\src")
setwd(src_dir)
compile("LIME.cpp")

dyn.load( dynlib("LIME") )


##--------------------------
## build inputs and object
##--------------------------

res <- run_LIME(modpath=NULL, 
				input=inputs_LC,
				data_avail="LC",
				Fpen=1,
				SigRpen=1,
				SigRprior=c(0.737,0.3),
				LFdist=1,
				C_type=0,
				est_more=FALSE,
				fix_more=FALSE,
				f_startval_ft=NULL,
				rdev_startval_t=NULL,
				est_selex_f=TRUE,
				randomR=TRUE,
				newtonsteps=FALSE,
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
plot_LCfits(LFlist=LF_list, 
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
