#' Calculate abundance at equilibrium
#'
#' \code{calc_equil_abund} calculates expected abundance given constant fishing mortality
#'
#' @param ages vector of ages
#' @param S_a selectivity-at-age
#' @param M natural mortality
#' @param F fishing mortality
#' @param R0 recruitment

#' @return vector of expected abundance at age
#' @export

calc_equil_abund <- function(ages, S_a, M, F, R0){
	N_a <- rep(NA, length(ages))
	N_a[1] <- R0
	for(i in 2:length(ages)){
		if(i<length(ages)) N_a[i] <- N_a[i-1]*exp(-M-S_a[i-1]*F)
		if(i==length(ages)) N_a[i] <- (N_a[i-1]*exp(-M-S_a[i-1]*F))/(1-exp(-M-S_a[i-1]*F))
	}
	return(N_a)
}