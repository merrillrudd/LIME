#' Calculate biological reference points
#'
#' \code{calc_derived_quants} calculates derived quanties for status or productivity
#'
#' @param Obj, the fitted TMB object
#' @param lh, list of life history information

#' @return List, a tagged list of potentially useful benchmarks
#' @export
calc_derived_quants = function( Obj, lh ){
  # Extract elements
  Data = Obj$env$data
  ParHat = Obj$env$parList()
  Report = Obj$report()

  if(max(Data$ages)<10){
    ## update life history with estimated values, but more refined time step for life history type
    lh_new <- with(lh, create_lh_list(vbk=vbk, linf=linf, lwa=lwa, lwb=lwb, S50=Report$S50, S95=Report$S95, selex_input="length", M50=ML50, M95=ML95, maturity_input="length", selex_type="logistic", binwidth=binwidth, t0=t0, CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=Report$sigma_R, SigmaF=SigmaF, R0=exp(Report$beta), h=Report$h, qcoef=Report$qcoef, M=M, AgeMax=AgeMax, start_ages=ages[1], nseasons=12))
  }
  if(max(Data$ages)>=10){
      lh_new <- with(lh, create_lh_list(vbk=vbk, linf=linf, lwa=lwa, lwb=lwb, S50=Report$S50, S95=Report$S95, selex_input="length", M50=ML50, M95=ML95, maturity_input="length", selex_type="logistic", dome_sd=dome_sd, binwidth=binwidth, t0=t0, CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=Report$sigma_R, SigmaF=SigmaF, R0=exp(Report$beta), h=Report$h, qcoef=Report$qcoef, M=M, AgeMax=AgeMax, start_ages=ages[1], nseasons=1))
  }

  SPR <- with(lh, calc_ref(ages=ages, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, F=Report$F_t[length(Report$F_t)], ref=FALSE))
  F30 <- tryCatch(with(lh_new, uniroot(calc_ref, lower=0, upper=50, ages=ages, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, ref=0.3)$root), error=function(e) NA) * lh_new$nseasons
  F40 <- tryCatch(with(lh_new, uniroot(calc_ref, lower=0, upper=50, ages=ages, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, ref=0.4)$root), error=function(e) NA) * lh_new$nseasons
  FF30 <- FF40 <- NULL
  if(is.na(F30)==FALSE) FF30 <- Report$F_y[length(Report$F_y)]/F30
  if(is.na(F40)==FALSE) FF40 <- Report$F_y[length(Report$F_y)]/F40

  # Total biomass
  TB_t = as.vector( Report$W_a %*% t(Report$N_ta) )
  Cw_t <- as.vector(Report$W_a %*% t(Report$Cn_ta))

  # MSY calculations
  fmsy <- optimize(calc_msy, ages=lh_new$ages, S_a=lh_new$S_a, M=lh_new$M, R0=exp(Report$beta), W_a=lh_new$W_a, Mat_a=lh_new$Mat_a, lower=0, upper=10, maximum=TRUE)$objective
  msy <- calc_msy(F=fmsy, ages=lh_new$ages, S_a=lh_new$S_a, M=lh_new$M, R0=exp(Report$beta), W_a=lh_new$W_a, Mat_a=lh_new$Mat_a)

  # Yield_Fn = function( Fmean, Return_type="Yield" ){
  #   # Modify data
  #   Data_new = Data
  #   Data_new[["n_t"]] = 1000
  #   Data_new[["n_c"]] = 1000
  #   Data_new[["n_i"]] = 1000
  #   Data_new[["n_lc"]] = 1000
  #   Data_new[["T_yrs"]] = 1:1000
  #   Data_new[["C_yrs"]] = 1:1000
  #   Data_new[["I_yrs"]] = 1:1000
  #   Data_new[["LC_yrs"]] = 1:1000
  #   Data_new[["obs_per_yr"]] = rep(Data[["obs_per_yr"]][1], 1000)
  #   Data_new[["RecDev_biasadj"]] = rep(0, Data_new[["n_t"]])
  #   Data_new[["C_t"]] = rep(1, Data_new[["n_t"]])
  #     names(Data_new[["C_t"]]) <- Data_new[["C_yrs"]]
  #   Data_new[["LF"]] = matrix(0, nrow=Data_new[["n_t"]], ncol=Data_new[["n_lb"]])
  #     rownames(Data_new[["LF"]]) <- Data_new[["LC_yrs"]]
  #   Data_new[["I_t"]] = cbind(1, rep(1,Data_new[["n_t"]]))
  #     names(Data_new[["I_t"]]) <- Data_new[["I_yrs"]]
  #   # Modify parameters
  #   ParHat_new = ParHat
  #   ParHat_new[["log_F_t_input"]] = rep( log(Fmean+1e-10), Data_new[["n_t"]])
  #   ParHat_new[["Nu_input"]] = rep(0, Data_new[["n_t"]])
  #   Obj_new = MakeADFun(data=Data_new[1:length(Data_new)], parameters=ParHat_new[1:length(ParHat_new)], inner.control=list(maxit=1e3), DLL="LIME" )
  #   # Extract and return stuff
  #   Report_new = Obj_new$report()
  #   if(Return_type=="Yield") Return = rev(Report_new$Cw_t_hat)[1]
  #   if(Return_type=="Report") Return = Report_new
  #   return(Return)
  # }

  # # Calculate Fmsy
  Fmsy = fmsy
  TBmsy = calc_equil_abund(ages=lh_new$ages, S_a=lh_new$S_a, M=lh_new$M, F=fmsy/lh_new$nseasons, R0=exp(Report$beta)) * lh_new$W_a
  SBmsy = calc_equil_abund(ages=lh_new$ages, S_a=lh_new$S_a, M=lh_new$M, F=fmsy/lh_new$nseasons, R0=exp(Report$beta)) * lh_new$W_a * lh_new$Mat_a
  MSY = msy
  TB0 = calc_equil_abund(ages=lh_new$ages, S_a=lh_new$S_a, M=lh_new$M, F=0, R0=exp(Report$beta)) * lh_new$W_a
  SB0 = calc_equil_abund(ages=lh_new$ages, S_a=lh_new$S_a, M=lh_new$M, F=0, R0=exp(Report$beta)) * lh_new$W_a * lh_new$Mat_a

  # Return
  Return <- list("SPR"=SPR, "F30"=F30, "F40"=F40, "FF30"=FF30, "FF40"=FF40, "Fmsy"=Fmsy, "FFmsy"=Report$F_t[length(Report$F_t)]/Fmsy, "SB0"=SB0, "TB0"=TB0, "TB_t"=TB_t, "SB_t"=Report$SB_t, "MSY"=MSY, "TBmsy"=TBmsy, "SBmsy"=SBmsy, "TBBmsy"=TB_t[length(TB_t)]/TBmsy, "SBBmsy"=Report$SB_t[length(Report$SB_t)]/SBmsy)
  return( Return )
}

