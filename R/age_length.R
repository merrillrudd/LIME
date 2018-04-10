#' Uncertainty in age given length
#'
#' \code{age_length} Calculates probability of being in a length bin given age
#' @author M.B. Rudd

#' @param highs upper length bins
#' @param lows lower length bins
#' @param L_a Predicted length-at-age
#' @param CVlen CV of predicted length at age
#' @importFrom stats pnorm
#' 
#' @return matrix - probability of being in a length bin given age
#' @export

age_length <- function(highs, lows, L_a, CVlen){
	lbprobs <- function(mnl, sdl) return(pnorm(highs, mnl, sdl) - pnorm(lows, mnl, sdl))
	vlprobs <- Vectorize(lbprobs, vectorize.args = c("mnl", "sdl"))
	plba <- t(vlprobs(L_a, L_a * CVlen))
	plba <- plba/rowSums(plba)
	return(plba)
}