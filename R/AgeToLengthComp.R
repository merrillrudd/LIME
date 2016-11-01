#' Age-converted-to-length structure
#'
#' \code{AgeToLengthComp} Converts vulnerable numbers at age to length composition of catch
#'
#' @param L_a vector, growth curve: predicted length at each age
#' @param CVlen coefficient of variation for growth curve
#' @param highs vector of upper end of length bins
#' @param lows vector of lower end of length bins
#' @param tyears number of years of data
#' @param N_at matrix of numbers in the population at each age over time
#' @param S_a vector of selectivity at age
#' @param comp_sample vector of number of individuals sampled each year (set as 1 for proportions)
#'
#' @return data frame - number of individuals in each length bin in each year
#' @export
AgeToLengthComp <- function(L_a, CVlen, highs, lows, tyears, N_at, S_a, comp_sample){

	################################################
	## Probability being in a length bin given age
	################################################
    lbprobs <- function(mnl,sdl) return(pnorm(highs,mnl,sdl)-pnorm(lows,mnl,sdl))
    vlprobs <- Vectorize(lbprobs,vectorize.args=c("mnl","sdl"))
    plba <- t(vlprobs(L_a,L_a*CVlen))
    plba <- plba/rowSums(plba)

    ################################################
	## Probability being in harvested at an age
	################################################
    page <- matrix(ncol=dim(plba)[1], nrow=tyears)
    for (y in 1:tyears) page[y,] <- N_at[,y] * S_a
    page <- page/rowSums(page)    

    ################################################
	## Probability of sampling a given length bin
	################################################
    plb <- matrix(ncol=length(highs), nrow=tyears)
    for (y in 1:tyears) plb[y,] <- page[y,] %*% plba
    plb <- plb/rowSums(plb)    

    #######################
	## Length frequencies 
	#######################
    LF <- array(0, dim=dim(plb))
      rownames(LF) <- 1:tyears
    for(y in 1:tyears){
    	LF[y,] <- rmultinom(n=1, size=comp_sample[y], prob=plb[y,])
    }

    Outs <- NULL
    Outs$plba <- plba
    Outs$plb <- plb
    Outs$page <- page
    Outs$LF <- LF

    return(Outs)
}