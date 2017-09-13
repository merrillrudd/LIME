# context("Life history")

# test_that("Monthly life history setup is working",{

# 	lh1 <- create_lh_list(vbk=0.21, 
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

# 	lh2 <- create_lh_list(vbk=0.21, 
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
# 						 nseasons=12)

# 	index <- which(lh2$ages %in% lh1$ages)
# 	L_index <- lh2$L_a[index]

# 	expect_true(length(lh2$ages)==length(lh1$ages)*12)
# 	expect_identical(lh1$L_a, L_index)

# })