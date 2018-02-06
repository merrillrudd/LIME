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
					 nseasons=1,
					 nfleets=2)

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
true <- generate_data(modpath=NULL,
					  itervec=1, 
					  Fdynamics=c("Constant","Endogenous"),
					  Rdynamics="Constant",
					  lh=lh,
					  Nyears=20,
					  Nyears_comp=c(20,10),
					  comp_sample=200,
					  init_depl=0.7,
					  seed=143,
					  fleet_percentage=c(0.5,0.5))

## Data input components
years <- true$years ## total years to model, can be 1:20 or 1998:2017
LF <- true$LF ## length composition data, years along rows and length bin along columns. year names should match 'year' quantity, e.g. 11:20 or 2008:2017 (can be any years within years to model)

## input data list
data_LF <- list("LF"=LF) ## length comp only

## plot length composition data
plot_LCfits(LFlist=lapply(1:dim(LF)[3], function(x) LF[,,x])) ## "Inputs" argument just must be a list with "LF" as one of the components, e.g. plot_LCfits(Inputs=list("LF"=true$LF))

## plot simulated data
dfsim <- true$dfsim
p <- ggplot(dfsim %>% filter(Variable %in% c("LengthComp", "TotalBiomass", "SpawningBiomass") == FALSE)) +
	geom_line(aes(x=X, y=Value, colour=Fleet), lwd=2) +
	facet_wrap(~Variable, scale='free_y') +
	expand_limits(y=0) +
	xlab("Year") +
	theme_lsd()

##----------------------------------------------------
## Step 3: Run Model
## ---------------------------------------------------
##-------------------------
## template file
##--------------------------
src_dir <- file.path("C:\\merrill\\LIME\\src")
setwd(src_dir)
compile("LIME.cpp")

dyn.load( dynlib("LIME") )


##--------------------------
## build inputs and object
##--------------------------

Data <- list("n_t"=dim(LF)[1],
			 "n_lb"=dim(LF)[2],
			 "n_fl"=dim(LF)[3],
			 "n_a"=length(lh$ages),
			 "LF_tlf"=LF,
			 "n_lc_ft"=true$obs_per_year,
			 "I_ft"=as.matrix(0),
			 "C_ft"=as.matrix(0),
			 "C_opt"=0,
			 "ML_ft"=as.matrix(0),
			 "ages"=lh$ages,
			 "L_a"=lh$L_a,
			 "W_a"=lh$W_a,
			 "M"=lh$M,
			 "h"=lh$h,
			 "Mat_a"=lh$Mat_a,
			 "lbhighs"=seq(lh$binwidth, by=lh$binwidth, length=dim(LF)[2]),
			 "lbmids"=seq(lh$binwidth/2, by=lh$binwidth, length=dim(LF)[2]),
			 "Fpen"=1,
			 "SigRpen"=1,
			 "SigRprior"=c(lh$SigmaR,0.3),
			 "selex_type_f"=c(1,1),
			 "LFdist"=1,
			 "S_yrs"=c(sapply(1:dim(LF)[1], function(x) rep(x, lh$nseasons))),
			 "n_s"=lh$nseasons,
			 "n_y"=ceiling(dim(LF)[1]/lh$nseasons))

Params <- list("log_F_ft"=matrix(log(1), nrow=Data$n_fl, ncol=Data$n_t),
				"log_q_f"=rep(log(lh$qcoef), Data$n_fl),
				"beta"=log(lh$R0),
				"log_sigma_R"=log(lh$SigmaR),
				"log_S50_f"=log(lh$SL50),
				"log_Sdelta_f"=log(lh$SL95 - lh$SL50),
				"log_sigma_F"=log(lh$SigmaF),
				"log_sigma_C"=log(lh$SigmaC),
				"log_sigma_I"=log(lh$SigmaI),
				"log_CV_L"=log(lh$CVlen),
				"log_theta"=log(rep(1, Data$n_fl)),
				"Nu_input"=rep(0, Data$n_t))

Map <- list()
    # Map[["logsigma"]] <- NA
    # Map[["logsigma"]] <- factor(Map[["logsigma"]])

Inputs <- list("Data"=Data, "Parameters"=Params, "Map"=Map, "Random"="Nu_input")
Obj <- MakeADFun( data=Inputs$Data, parameters=Inputs$Parameters, map=Inputs$Map, Random=Inputs$Random, DLL="LIME")


Upr <- rep(Inf, length(Obj$par))
Upr[which(names(Obj$par)=="log_S50_f")] <- log(lh$linf)


Opt <- TMBhelper::Optimize( obj=Obj, loopnum=3, upper=Upr )
# Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, upper=Upr, lower=Lwr)
Opt[["final_gradient"]] = Obj$gr( Opt$opt$par ) 
df <- Opt[["final_gradient"]]
colnames(df) <- names(Obj$par)
which(abs(df) > 0.001)
 
Report <- Obj$report()
Sdreport <- sdreport(Obj)

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
plot_LCfits(LFlist=lapply(1:dim(LF)[3], function(x) LF[,,x]), 
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
