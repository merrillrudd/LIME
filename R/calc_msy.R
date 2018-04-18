#' Calculate MSY reference point
#'
#' \code{calc_msy} calculates yield per recruit, when maximized across levels of F produces Fmsy
#'
#' @author M.B. Rudd
#' @param F fishing mortality rate
#' @param ages ages in the population
#' @param M natural mortality
#' @param R0 equilibrium recruits
#' @param W_a weight-at-age

#' @return yield per recruit, in biomass
#' @export
calc_msy <- function(F, ages, M, R0, W_a){

	Nage <- calc_equil_abund(ages=ages, M=M, F=F, R0=R0)
	YPR <- sum(Nage * W_a * (1 - exp(-M - F)) * (F) / (M + F))

	return(YPR)
}