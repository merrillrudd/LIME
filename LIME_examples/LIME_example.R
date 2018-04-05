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

###************************************
## Section 1: Single fleet
###************************************
##----------------------------------------------------------------
## Step 1: Specify biological inputs and parameter starting values
##----------------------------------------------------------------
## single fleet
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
					 SigmaR=0.737,
					 SigmaF=0.2,
					 SigmaC=0.1,
					 SigmaI=0.1,
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
	xlab("Age")

## by length
p <- ggplot(lh$df %>% 
			filter(Variable %in% c("Selectivity","Maturity")) %>% 
			filter(By == "Length")) + 
	geom_line(aes(x=X, y=Value, color=Fleet), lwd=2) +
	facet_grid(Variable~., scale="free_y") +
	xlab("Length")


##----------------------------------------------------
## Step 2: Setup data input
## ---------------------------------------------------
#######################################
## Simulation feature
#######################################
true <- generate_data(modpath=NULL,
					  itervec=1, 
					  Fdynamics=c("Endogenous"),
					  Rdynamics="Constant",
					  lh=lh,
					  Nyears=20,
					  Nyears_comp=c(20),
					  comp_sample=200,
					  init_depl=0.7,
					  seed=132,
					  fleet_proportions=1)


## plot simulated data
par(mfrow=c(3,2))
plot(true$SPR_t, type="l", ylim=c(0,1), lwd=2, xlab="Time", ylab="SPR")
plot(true$R_t, type="l", ylim=c(0,3), lwd=2, xlab="Time", ylab="Recruitment")
plot(x=1,y=1,type="n", ylim=c(0,1), xlim=c(1,length(true$SPR_t)), xlab="Time", ylab="Fishing mortality")
lty <- ifelse(lh$nfleets==1,1,2)
for(f in 1:lh$nfleets){
	lines(true$F_ft[f,], lwd=2, lty=lty)
}
plot(true$D_t, type="l", ylim=c(0,2), lwd=2, xlab="Time", ylab="Relative spawning biomass")
plot(x=1, y=1, type="n", ylim=c(0,max(true$Cw_ft)), xlim=c(1,length(true$SPR_t)), xlab="Time", ylab="Catch (biomass)")
for(f in 1:lh$nfleets){
	lines(true$Cw_ft[f,], lwd=2, lty=lty)
}
plot(x=1, y=1, type="n", ylim=c(0,max(true$I_ft)), xlim=c(1,length(true$SPR_t)), xlab="Time", ylab="Abundance index")
for(f in 1:lh$nfleets){
	lines(true$I_ft[f,], lwd=2, lty=lty)
}

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
# LF_df <- true$dfsim %>% filter(Variable == "LengthComp") ## long-form data frame where "X" = year, "Value"=length measurement, and "Fleet"=discrete variables representing a fleet. 

## example with length data only
data_LF <- list("years"=1:true$Nyears, "LF"=LF_array)

##if using multinomial distribution, must specify annual effective sample size by fleet
data_LF_neff <- list("years"=1:true$Nyears, "LF"=LF_array, "neff_ft"=true$obs_per_year)

## create model inputs with life history information and data
## outputs length data as array
inputs_LC <- create_inputs(lh=lh, input_data=data_LF)

#######################################
## Other data type input options
#######################################
data_all <- list("years"=1:true$Nyears, "LF"=LF_array, "I_ft"=true$I_ft, "C_ft"=true$Cw_ft, "neff_ft"=true$obs_per_year)
inputs_all <- create_inputs(lh=lh, input_data=data_all)

##----------------------------------------------------
## Step 3: Run model
## ---------------------------------------------------
#######################################
## Data-rich test
#######################################
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



##----------------------------------------------------
## Step 4: Plot results
## ---------------------------------------------------
## plot length composition data and fits
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
			set_ylim=list("SPR" = c(0,1)))

#######################################
## Length-data only
#######################################

lc_only <- run_LIME(modpath=NULL, 
				input=inputs_LC,
				data_avail="LC")

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
			set_ylim=list("Fish" =c(0,2), "SPR" = c(0,1), "SB"=c(0,2)))		

#######################################
## Catch + length data
#######################################
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
			set_ylim=list("SPR" = c(0,1), "SB"=c(0,2)))		

###************************************
## Section 2: Multiple fleets
###************************************
##----------------------------------------------------------------
## Step 1: Specify biological inputs and parameter starting values
##----------------------------------------------------------------
## specify nfleets=<number of fleets>, and add starting values for selectivity, selectivity type for nfleets
lh <- create_lh_list(vbk=0.21, 
					 linf=65, 
					 t0=-0.01,
					 lwa=0.0245, 
					 lwb=2.79, 
					 S50=c(20,30), 
					 S95=c(26,36), 
					 selex_input="length",
					 selex_type=c("logistic","logistic"),
					 M50=34,
					 M95=NULL,
					 maturity_input="length",
					 M=0.27, 
					 h=0.7,
					 binwidth=1,
					 CVlen=0.1,
					 SigmaR=0.737,
					 SigmaF=0.2,
					 SigmaC=0.1,
					 SigmaI=0.1,
					 R0=1,
					 qcoef=1e-5,
					 start_ages=0,
					 rho=0.43,
					 nseasons=1,
					 nfleets=2)

##----------------------------------------------------
## Step 2: Setup data input
## ---------------------------------------------------
#######################################
## Simulation feature
#######################################
## specify Fdynamics, number of years of composition data, composition samples annually from each fleet, and fleet_proportions (e.g. 60% of catch from fleet 1, 40% of catch from fleet 2, although catch may  not be observed)
true <- generate_data(modpath=NULL,
					  itervec=1, 
					  Fdynamics=c("Constant","Endogenous"),
					  Rdynamics="Constant",
					  lh=lh,
					  Nyears=20,
					  Nyears_comp=c(20,10),
					  comp_sample=200,
					  init_depl=0.7,
					  seed=132,
					  fleet_proportions=c(0.6,0.4))

## plot simulated data
par(mfrow=c(3,2))
plot(true$SPR_t, type="l", ylim=c(0,1), lwd=2, xlab="Time", ylab="SPR")
plot(true$R_t, type="l", ylim=c(0,3), lwd=2, xlab="Time", ylab="Recruitment")
plot(x=1,y=1,type="n", ylim=c(0,1), xlim=c(1,length(true$SPR_t)), xlab="Time", ylab="Fishing mortality")
lty <- ifelse(lh$nfleets==1,1,2)
for(f in 1:lh$nfleets){
	lines(true$F_ft[f,], lwd=2, lty=lty)
}
plot(true$D_t, type="l", ylim=c(0,2), lwd=2, xlab="Time", ylab="Relative spawning biomass")
plot(x=1, y=1, type="n", ylim=c(0,max(true$Cw_ft)), xlim=c(1,length(true$SPR_t)), xlab="Time", ylab="Catch (biomass)")
for(f in 1:lh$nfleets){
	lines(true$Cw_ft[f,], lwd=2, lty=lty)
}
plot(x=1, y=1, type="n", ylim=c(0,max(true$I_ft)), xlim=c(1,length(true$SPR_t)), xlab="Time", ylab="Abundance index")
for(f in 1:lh$nfleets){
	lines(true$I_ft[f,], lwd=2, lty=lty)
}

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
# LF_df <- true$dfsim %>% filter(Variable == "LengthComp") ## long-form data frame where "X" = year, "Value"=length measurement, and "Fleet"=discrete variables representing a fleet. 

## example with length data only
data_LF <- list("years"=1:true$Nyears, "LF"=LF_array)

##if using multinomial distribution, must specify annual effective sample size by fleet
data_LF_neff <- list("years"=1:true$Nyears, "LF"=LF_array, "neff_ft"=true$obs_per_year)

## create model inputs with life history information and data
## outputs length data as array
inputs_LC <- create_inputs(lh=lh, input_data=data_LF)

#######################################
## Other data type input options
#######################################
data_all <- list("years"=1:true$Nyears, "LF"=LF_array, "I_ft"=true$I_ft, "C_ft"=true$Cw_ft, "neff_ft"=true$obs_per_year)
inputs_all <- create_inputs(lh=lh, input_data=data_all)

##----------------------------------------------------
## Step 3: Run model
## ---------------------------------------------------
#######################################
## Data-rich test
#######################################
## dirichlet-multinomial (LFdist=1) currently not working with multiple fleets
## but with more data types than length composition, can estimate F by fleet
rich_mf <- run_LIME(modpath=NULL, 
				input=inputs_all,
				data_avail="Index_Catch_LC",
				C_type=2,
				LFdist=0)

## check TMB inputs
Inputs <- rich_mf$Inputs

## Report file
Report <- rich_mf$Report

## Standard error report
Sdreport <- rich_mf$Sdreport

## check convergence
hessian <- Sdreport$pdHess
gradient <- rich_mf$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE

##----------------------------------------------------
## Step 4: Plot results
## ---------------------------------------------------
## plot length composition data and fits
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
			set_ylim=list("SPR" = c(0,1)))

#######################################
## Length-data only
#######################################
## dirichlet-multinomial (LFdist=1) currently not working with multiple fleets
## with length composition only, model does not converge in estimating F for each fleet
## must estimate total F and specify fleet proportions
lc_only_mf <- run_LIME(modpath=NULL, 
				input=inputs_LC,
				data_avail="LC",
				LFdist=0,
				est_totalF=TRUE,
				prop_f = c(0.6,0.4))

## check TMB inputs
Inputs <- lc_only_mf$Inputs

## Report file
Report <- lc_only_mf$Report

## Standard error report
Sdreport <- lc_only_mf$Sdreport

## check convergence
hessian <- Sdreport$pdHess
gradient <- lc_only_mf$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE


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
			set_ylim=list("Fish" =c(0,2), "SPR" = c(0,1), "SB"=c(0,2)))		

#######################################
## Catch + length data
#######################################
catch_lc_mf <- run_LIME(modpath=NULL, 
				input=inputs_all,
				data_avail="Catch_LC",
				C_type=2,
				LFdist=0)

## check TMB inputs
Inputs <- catch_lc_mf$Inputs

## Report file
Report <- catch_lc_mf$Report

## Standard error report
Sdreport <- catch_lc_mf$Sdreport

## check convergence
hessian <- Sdreport$pdHess
gradient <- catch_lc_mf$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE


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
			set_ylim=list("SPR" = c(0,1), "SB"=c(0,2)))		

#######################################
## Index + length data
#######################################
index_lc_mf <- run_LIME(modpath=NULL, 
				input=inputs_all,
				data_avail="Index_LC",
				LFdist=0,
				est_totalF=TRUE,
				prop_f=c(0.6,0.4))


## check TMB inputs
Inputs <- index_lc_mf$Inputs

## Report file
Report <- index_lc_mf$Report

## Standard error report
Sdreport <- index_lc_mf$Sdreport

## check convergence
hessian <- Sdreport$pdHess
gradient <- index_lc_mf$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE


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
			set_ylim=list("SPR" = c(0,1), "SB"=c(0,2)))		




##----------------------------------------------------
## Check multiple seasons
## ---------------------------------------------------
## life history list
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
					 h=0.7,
					 binwidth=1,
					 CVlen=0.1,
					 SigmaR=0.737,
					 SigmaF=0.2,
					 SigmaC=0.1,
					 SigmaI=0.1,
					 R0=1,
					 qcoef=1e-5,
					 start_ages=0,
					 rho=0.43,
					 nseasons=4,
					 nfleets=1)

## generate data
true <- generate_data(modpath=NULL,
					  itervec=1, 
					  Fdynamics=c("Constant"),
					  Rdynamics="Constant",
					  lh=lh,
					  Nyears=20,
					  Nyears_comp=20,
					  comp_sample=200,
					  init_depl=0.7,
					  seed=132,
					  fleet_proportions=1,
					  pool=TRUE)


## plot simulated data
par(mfrow=c(3,2))
plot(true$SPR_t, type="l", ylim=c(0,1), lwd=2, xlab="Time", ylab="SPR")
plot(true$R_t, type="l", ylim=c(0,3), lwd=2, xlab="Time", ylab="Recruitment")
plot(x=1,y=1,type="n", ylim=c(0,1), xlim=c(1,length(true$SPR_t)), xlab="Time", ylab="Fishing mortality")
lty <- ifelse(lh$nfleets==1,1,2)
for(f in 1:lh$nfleets){
	lines(true$F_ft[f,], lwd=2, lty=lty)
}
plot(true$D_t, type="l", ylim=c(0,2), lwd=2, xlab="Time", ylab="Relative spawning biomass")
plot(x=1, y=1, type="n", ylim=c(0,max(true$Cw_ft)), xlim=c(1,length(true$SPR_t)), xlab="Time", ylab="Catch (biomass)")
for(f in 1:lh$nfleets){
	lines(true$Cw_ft[f,], lwd=2, lty=lty)
}
plot(x=1, y=1, type="n", ylim=c(0,max(true$I_ft)), xlim=c(1,length(true$SPR_t)), xlab="Time", ylab="Abundance index")
for(f in 1:lh$nfleets){
	lines(true$I_ft[f,], lwd=2, lty=lty)
}

# #######################################
# ## Length comp data input options
# #######################################
# ## Option 1: Length comp array
# LF_array <- true$LF ## array with rows = years, columns = upper length bins, 3rd dimension = fleets

# ## Option 2: Length comp list
# LF_list <- lapply(1:lh$nfleets, function(x) true$LF[,,x]) ##list with 1 element per fleet, and each element is a matrix with rows = years, columns = upper length bins

# 	## plot length composition data using LF_list
# 	plot_LCfits(LFlist=LF_list, ylim=c(0,0.15)) ## "Inputs" argument just must be a list with "LF" as one of the components, e.g. plot_LCfits(Inputs=list("LF"=true$LF))

# ## Option 3: Data frame
# # LF_df <- true$dfsim %>% filter(Variable == "LengthComp") ## long-form data frame where "X" = year, "Value"=length measurement, and "Fleet"=discrete variables representing a fleet. 

# ## example with length data only
# data_LF <- list("years"=1:true$Nyears, "LF"=LF_array)

# ##if using multinomial distribution, must specify annual effective sample size by fleet
# data_LF_neff <- list("years"=1:true$Nyears, "LF"=LF_array, "neff_ft"=true$obs_per_year)

# ## create model inputs with life history information and data
# ## outputs length data as array
# inputs_LC <- create_inputs(lh=lh, input_data=data_LF)

# #######################################
# ## Other data type input options
# #######################################
# data_all <- list("years"=1:true$Nyears, "LF"=LF_array, "I_ft"=true$I_ft, "C_ft"=true$Cw_ft, "neff_ft"=true$obs_per_year)
# inputs_all <- create_inputs(lh=lh, input_data=data_all)

# ## dirichlet-multinomial not currently working for nfleets > 1
# rich <- run_LIME(modpath=NULL, 
# 				input=inputs_all,
# 				data_avail="Index_Catch_LC",
# 				Fpen=1,
# 				SigRpen=1,
# 				SigRprior=c(0.737,0.3),
# 				LFdist=0,
# 				C_type=2,
# 				est_more=FALSE,
# 				fix_more=FALSE,
# 				mirror=NULL,
# 				f_startval_ft=NULL,
# 				rdev_startval_t=NULL,
# 				est_selex_f=TRUE,
# 				randomR=TRUE,
# 				newtonsteps=3,
# 				F_up=10,
# 				S50_up=lh$linf,
# 				derive_quants=FALSE,
# 				itervec=NULL,
# 				rewrite=TRUE,
# 				simulation=FALSE,
# 				est_totalF=FALSE,
# 				prop_f=1)

# ## check TMB inputs
# Inputs <- rich$Inputs

# ## Report file
# Report <- rich$Report

# ## Standard error report
# Sdreport <- rich$Sdreport


# ## check convergence
# hessian <- Sdreport$pdHess
# gradient <- rich$opt$max_gradient <= 0.001
# check <- hessian == TRUE & gradient == TRUE

# ## plot length composition data
# plot_LCfits(LFlist=LF_list, 
# 			Inputs=Inputs, 
# 			Report=Report,
# 			ylim=c(0,0.2))

# ## plot model output
# plot_output(Inputs=Inputs, 
# 			Report=Report,
# 			Sdreport=Sdreport, 
# 			lh=lh,
# 			True=true, 
# 			plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
# 			set_ylim=list("SPR" = c(0,1)))


