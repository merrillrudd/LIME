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
					  Fdynamics="Constant",
					  Rdynamics="AR",
					  lh=lh,
					  Nyears=20,
					  comp_sample=200,
					  init_depl=0.4)




#######################################################################
## ---------------- Model settings and directories ------------------
#######################################################################

## data availability scenarios (can also set up for "Catch_LC" when both catch and length comp are available -- just need to make sure the data type is in the name -- with Index, Catch, and LC being the options)
# avail_set <- c("Index_LC", "LC")
avail_set <- "LC"

## estimate variances -- always estimating Recruitment variation (log_sigma_R), but could add on other variance parameters (match variance names exactly ** update manual) -- in this case could estimate the CV for the growth curve 
## these tags are used in the directory names and what they mean can be specified when model is run
estsigma_set <- c("RecVar", "RecGrowthVars")

## setup combinations of models to run
modcombos <- as.matrix(expand.grid("Data_avail"=avail_set, "Est_variance"=estsigma_set))

## setup results directory
res_dir <- file.path(dir, "results")
dir.create(res_dir, showWarnings=FALSE)

## transform model combinations into directory names
alldirs <- model_paths(modcombos=modcombos, res_dir=res_dir)

#######################################################################
## ---------------- Assessment model ------------------
#######################################################################

start_run <- Sys.time()

## loop over possible models
for(dd in 1:length(alldirs)){

	## get available data types from model path name, used for formatting TMB input
	data_avail <- ifelse(grepl("Index", alldirs[dd]), avail_set[which(grepl("Index", avail_set))], avail_set[which(grepl("Index", avail_set)==FALSE)])

	## get variance parameters to estimate from model name, used for formatting TMB input
	id_sigma <- estsigma_set[which(sapply(1:length(estsigma_set), function(x) grepl(estsigma_set[x], alldirs[dd])))]
	if(id_sigma=="RecVar") est_sigma <- "log_sigma_R"
	if(id_sigma=="RecGrowthVars") est_sigma <- c("log_sigma_R", "log_CV_L")
	if(id_sigma=="RecGrowthIndexVars") est_sigma <- c("log_sigma_R", "log_CV_L", "log_sigma_I")

	## run assessment model and save final gradients, parameter names, and estimates to check convergence to directory
	ignore <- runModel(modpath=alldirs[dd], itervec=NULL, est_sigma=est_sigma, data_avail=data_avail, lh_list=lh, rewrite=TRUE, start_f=0, simulation=FALSE, input_data=data_new)
}

end_run <- Sys.time() - start_run

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