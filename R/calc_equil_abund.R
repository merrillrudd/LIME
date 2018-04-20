#' Calculate abundance at equilibrium
#'
#' \code{calc_equil_abund} calculates expected abundance given constant fishing mortality
#'
#' @author M.B. Rudd
#' @param ages vector of ages
#' @param M natural mortality
#' @param F fishing mortality
#' @param S_fa selectivity by fleet by age
#' @param R0 recruitment

#' @return vector of expected abundance at age
#' @export

calc_equil_abund <- function(ages, M, F, S_fa, R0){
	N_a <- rep(NA, length(ages))
	N_a[1] <- R0

	nfleets <- nrow(S_fa)
	Fmat <- matrix(NA, nrow=nfleets, ncol=length(ages))
	for(i in 1:nfleets){
		Fmat[i,] <- S_fa[i,]*F
	}
	Ftotal <- colSums(Fmat)

	for(i in 2:length(ages)){
		if(i<length(ages)) N_a[i] <- N_a[i-1]*exp(-M-Ftotal[i])
		if(i==length(ages)) N_a[i] <- (N_a[i-1]*exp(-M-Ftotal[i]))/(1-exp(-M-Ftotal[i]))
	}
	return(N_a)
}