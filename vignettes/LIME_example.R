rm(list=ls())

## Packages

devtools::install_github("merrillrudd/LIME", build.vignettes=TRUE, dependencies=TRUE)
library(LIME)

devtools::install_github("kaskr/TMB_contrib_R/TMBhelper", dep=TRUE)
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

##----------------------------------------------------
## Step 2: Setup data input
## ---------------------------------------------------

## Demonstrate data generation option
true <- generate_data(modpath=NULL,
					  data_avail="Index_Catch_LC",
					  itervec=1, 
					  Fdynamics="Ramp",
					  Rdynamics="AR",
					  lh=lh,
					  Nyears=20,
					  Nyears_comp=10,
					  comp_sample=200,
					  init_depl=0.4)

## Data input components
years <- true$years ## total years to model, can be 1:20 or 1998:2017
LF <- true$LF ## length composition data, years along rows and length bin along columns. year names should match 'year' quantity, e.g. 11:20 or 2008:2017 (can be any years within years to model)
C_t <- true$C_t ## (optional) catch data, with elements of vector named with the year observed, e.g. 1:20 or 1998:2017 (can be any years within years to model)
I_t <- true$I_t ## (optional) abundance index, with elements of vector named with the year observed, e.g. 1:20 or 1998:2017 (can be any years within years to model)

## input data list
data_LF <- list("years"=years, "LF"=LF) ## length comp only
data_LF_Catch <- list("years"=years, "LF"=LF, "C_t"=C_t) ## length comp + catch
data_LF_Index <- list("years"=years, "LF"=LF, "I_t"=I_t) ## length comp + index

## plot length composition data
plot_LCfits(Inputs=data_LF) ## "Inputs" argument just must be a list with "LF" as one of the components, e.g. plot_LCfits(Inputs=list("LF"=true$LF))

##----------------------------------------------------
## Step 3: Run Model
## ---------------------------------------------------

## run LIME - may take a few minutes
## looking for outer mgc to minimize and ustep moving towards 1 for well-behaved model

## length comp only
res <- run_LIME(modpath=NULL,
				lh=lh,
				input_data=data_LF,
				est_sigma="log_sigma_R", 
				data_avail="LC")

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
			dim=c(5,2), 
			n=FALSE)

## plot model output
plot_output(all_years=1:20,
			lc_years=11:20, 
			Inputs=Inputs, 
			Report=Report, 
			Sdreport=Sdreport, 
			lh=lh, 
			true_years=1998:2017, 
			True=true, 
			plot=c("Fish","Rec","SPR","ML","SB","Selex"))






## length comp + catch
res <- run_LIME(modpath=NULL,
				lh=lh,
				input_data=data_LF_Catch,
				est_sigma="log_sigma_R", 
				data_avail="Catch_LC")

## length comp + index
res <- run_LIME(modpath=NULL,
				lh=lh,
				input_data=data_LF_Index,
				est_sigma="log_sigma_R", 
				data_avail="Index_LC")

#######################################################################
## ---------------- model comparison ------------------
#######################################################################

aic <- calc_AIC(modpath_vec=alldirs)

## move forward with model 1 -- including abundance index and fixing growth CV

#######################################################################
## ---------------- Figures ---------------------------
#######################################################################

## figure directory setup
fig_dir <- file.path(res_dir, "figs")
dir.create(fig_dir, showWarnings=FALSE)

## results of chosen model
dd <- 1 
choose_dir <- alldirs[dd]
	Inputs <- readRDS(file.path(choose_dir, "Inputs2.rds"))
	Report <- readRDS(file.path(choose_dir, "Report.rds"))
	Sdreport <- readRDS(file.path(choose_dir, "Sdreport.rds"))
	Quants <- readRDS(file.path(choose_dir, "Derived_quants.rds"))
	flag <- ifelse(file.exists(file.path(choose_dir, "NAs_final_gradient.txt"))|file.exists(file.path(choose_dir, "high_final_gradient.txt")), TRUE, FALSE)

	## save=FALSE displays plot
	## save=TRUE saves into figure directory
	LIME_fits(Inputs, Report, Sdreport, data, save=TRUE)

### length composition model fits
png(file.path(fig_dir, "LC_fits.png"), width=8, height=7, res=200, units="in")
par(mfrow=c(3,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
LC_yrs <- Inputs$Data$LC_yrs
T_yrs <- Inputs$Data$T_yrs
for(t in 1:length(LC_yrs)){
	plot(data$LFprop[t,], pch=19, lwd=2, ylim=c(0, max(data$LFprop)), xaxt="n", yaxt="n")
	lines(Report$plb[LC_yrs[t],], lwd=4, col="red")
	if(t %in% c(7,8,9)) axis(1, cex.axis=1.2)
	if(t %in% c(1,4,7)) axis(2, cex.axis=1.2)
	print.letter(data$years[t]+10, c(0.1,0.925), cex=1.5)
}
mtext("Year", side=1, line=3.5, cex=1.5, outer=TRUE)
mtext("Proportion", side=2, line=3.5, cex=1.5, outer=TRUE)
dev.off()


#######################################################################
## ---------------- Sensitivity ------------------------
#######################################################################

sens_dir <- file.path(dir, "sensitivities")
dir.create(sens_dir, showWarnings=TRUE)

## natural mortality directory
M_sens_dir <- file.path(sens_dir, "M")
dir.create(M_sens_dir, showWarnings=FALSE)

## vector of possible M values
M_vec <- seq(0.1,0.8, by=0.05)

## loop over possible M values
for(mm in 1:length(M_vec)){

	dir_new <- file.path(M_sens_dir, M_vec[mm])
	dir.create(dir_new, showWarnings=FALSE)

	## can use this function to adjust some parameters
	## make own list of life history with the same values in output list
	lh_new <- choose_lh_list(species="CRSNAP", selex="asymptotic", param_adjust=c("CVlen","M"), val=c(0.1, M_vec[mm]))

	data_avail <- "Index_LC"
	est_sigma <- c("log_sigma_R")
	ignore <- runModel(modpath=dir_new, itervec=NULL, est_sigma=est_sigma, data_avail=data_avail, lh_list=lh_new, rewrite=TRUE, start_f=0, simulation=FALSE, input_data=data)
}

## likelihood profile
M_prof <- rep(NA, length(M_vec))
for(cc in 1:length(M_vec)){
	rep <- readRDS(file.path(M_sens_dir, M_vec[cc], "Report.rds"))
	M_prof[cc] <- rep$jnll
}
plot(x=M_vec, y=M_prof, pch=19, ylim=c(0, max(M_prof)*1.05), xlab="Value of M", ylab="NLL")
points(x=cr_lh$M, y=Report$jnll, col="red", pch=19)
abline(h=min(M_prof), lty=2)