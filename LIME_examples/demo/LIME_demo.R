rm(list=ls())

## Packages
devtools::install_github("merrillrudd/LIME", dependencies=TRUE, ref="multifleet")
library(LIME)

devtools::install_github("kaskr/TMB_contrib_R/TMBhelper", dependencies=TRUE)
library(TMBhelper)

##----------------------------------------------------------------
## Step 1: Read in length data
##----------------------------------------------------------------
## single fleet
lh <- create_lh_list(vbk=0.117, 
					 linf=37.4, 
					 t0=-0.01,
					 lwa=3.4e-5, 
					 lwb=2.87, 
					 M50=26,
					 maturity_input="length",
					 M=0.119, 
					 S50=c(20), 
					 S95=c(26), 
					 selex_input="length",
					 selex_type=c("logistic"),
					 CVlen=0.1,
					 SigmaR=0.5,
					 SigmaF=0.1,
					 binwidth=1,
					 nfleets=1)