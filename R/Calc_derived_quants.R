#' Calculate biological reference points
#'
#' \code{calc_derived_quants} calculates derived quanties for status or productivity
#'
#' @param Obj, the fitted TMB object

#' @return List, a tagged list of potentially useful benchmarks
#' @export
calc_derived_quants = function( Obj ){
  # Extract elements
  Data = Obj$env$data
  ParHat = Obj$env$parList()
  Report = Obj$report()

  SPR <- with(Report, calc_ref(Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, F=F_t[length(F_t)], ref=FALSE))
  F30 <- tryCatch(with(Report, uniroot(calc_ref, lower=0, upper=50, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, ref=0.3)$root), error=function(e) NA)
  F40 <- tryCatch(with(Report, uniroot(calc_ref, lower=0, upper=50, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, ref=0.4)$root), error=function(e) NA)
  FF30 <- FF40 <- NULL
  if(is.na(F30)==FALSE) FF30 <- Report$F_t[length(Report$F_t)]/F30
  if(is.na(F40)==FALSE) FF40 <- Report$F_t[length(Report$F_t)]/F40

  # Total biomass
  TB_t = as.vector( Report$W_a %*% t(Report$N_ta) )
  Cw_t <- as.vector(Report$W_a %*% t(Report$Cn_ta));

  # MSY calculations
  Yield_Fn = function( Fmean, Return_type="Yield" ){
    # Modify data
    Data_new = Data
    Data_new[["n_t"]] = 1000
    Data_new[["n_c"]] = 1000
    Data_new[["n_i"]] = 1000
    Data_new[["n_lc"]] = 1000
    Data_new[["T_yrs"]] = 1:1000
    Data_new[["C_yrs"]] = 1:1000
    Data_new[["I_yrs"]] = 1:1000
    Data_new[["LC_yrs"]] = 1:1000
    Data_new[["obs_per_yr"]] = rep(Data[["obs_per_yr"]][1], 1000)
    Data_new[["RecDev_biasadj"]] = rep(0, Data_new[["n_t"]])
    Data_new[["C_t"]] = rep(1, Data_new[["n_t"]])
      names(Data_new[["C_t"]]) <- Data_new[["C_yrs"]]
    Data_new[["LF"]] = matrix(0, nrow=Data_new[["n_t"]], ncol=Data_new[["n_lb"]])
      rownames(Data_new[["LF"]]) <- Data_new[["LC_yrs"]]
    Data_new[["I_t"]] = cbind(1, rep(1,Data_new[["n_t"]]))
      names(Data_new[["I_t"]]) <- Data_new[["I_yrs"]]
    # Modify parameters
    ParHat_new = ParHat
    ParHat_new[["log_F_t_input"]] = rep( log(Fmean+1e-10), Data_new[["n_t"]])
    ParHat_new[["Nu_input"]] = rep(0, Data_new[["n_t"]])
    Obj_new = MakeADFun(data=Data_new[1:length(Data_new)], parameters=ParHat_new[1:length(ParHat_new)], inner.control=list(maxit=1e3), DLL="LIME" )
    # Extract and return stuff
    Report_new = Obj_new$report()
    if(Return_type=="Yield") Return = rev(Report_new$Cw_t_hat)[1]
    if(Return_type=="Report") Return = Report_new
    return(Return)
  }

  # Calculate Fmsy
  Fmsy = optimize( f=Yield_Fn, interval=c(0,5), maximum=TRUE)$maximum
  Report_msy = Yield_Fn( Fmean=Fmsy, Return_type="Report" )
  Report_0 = Yield_Fn( Fmean=0, Return_type="Report" )
  TBmsy = rev(as.vector(Report$W_a %*% t(Report_msy$N_ta)))[1]
  SBmsy = rev(Report_msy$SB_t)[1]
  MSY = rev(Report_msy$Cw_t_hat)[1]
  TB0 = rev(as.vector(Report$W_a %*% t(Report_0$N_ta)))[1]
  SB0 = rev(Report_0$SB_t)[1]

  # Return
  Return <- list("SPR"=SPR, "F30"=F30, "F40"=F40, "FF30"=FF30, "FF40"=FF40, "Fmsy"=Fmsy, "FFmsy"=Report$F_t[length(Report$F_t)]/Fmsy, "SB0"=SB0, "TB0"=TB0, "TB_t"=TB_t, "SB_t"=Report$SB_t, "MSY"=MSY, "TBmsy"=TBmsy, "SBmsy"=SBmsy, "TBBmsy"=TB_t[length(TB_t)]/TBmsy, "SBBmsy"=Report$SB_t[length(Report$SB_t)]/SBmsy)
  return( Return )
}

