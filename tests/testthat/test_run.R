# context("Running LIME")

# test_that("final gradient is below 0.001", {

# 	lh <- create_lh_list(vbk=0.21, 
# 						 linf=65, 
# 						 t0=-0.01,
# 						 lwa=0.0245, 
# 						 lwb=2.79, 
# 						 S50=20, 
# 						 S95=26, 
# 						 selex_input="length",
# 						 M50=34,
# 						 M95=NULL,
# 						 maturity_input="length",
# 						 M=0.27, 
# 						 binwidth=1,
# 						 CVlen=0.1,
# 						 SigmaR=0.737,
# 						 SigmaF=0.2,
# 						 SigmaC=0.2,
# 						 SigmaI=0.2,
# 						 R0=1,
# 						 qcoef=1e-5,
# 						 start_ages=0,
# 						 rho=0.43,
# 						 nseasons=1)	

# 	## Demonstrate data generation option
# 	set.seed(143)
# 	true <- generate_data(modpath=NULL,
# 						  data_avail="Index_Catch_LC",
# 						  itervec=1, 
# 						  Fdynamics="Ramp",
# 						  Rdynamics="AR",
# 						  lh=lh,
# 						  Nyears=20,
# 						  Nyears_comp=10,
# 						  comp_sample=200,
# 						  init_depl=0.4)	

# 	## Data input components
# 	years <- true$years ## total years to model, can be 1:20 or 1998:2017
# 	LF <- true$LF ## length composition data, years along rows and length bin along columns. year names should match 'year' quantity, e.g. 11:20 or 2008:2017 (can be any years within years to model)
# 	C_t <- true$C_t ## (optional) catch data, with elements of vector named with the year observed, e.g. 1:20 or 1998:2017 (can be any years within years to model)
# 	I_t <- true$I_t ## (optional) abundance index, with elements of vector named with the year observed, e.g. 1:20 or 1998:2017 (can be any years within years to model)	

# 	## input data list
# 	data_LF <- list("years"=years, "LF"=LF) ## length comp only

# 	res <- run_LIME(modpath=NULL,
# 				lh=lh,
# 				input_data=data_LF,
# 				est_sigma="log_sigma_R", 
# 				data_avail="LC")

# 	check <- res$df
# 	expect_false(all(is.na(Report)))
# 	expect_false(all(is.na(Sdreport)))
# 	expect_false(all(is.na(res$opt)))
# 	expect_true(max(abs(check[,1]))<0.001)
# })

