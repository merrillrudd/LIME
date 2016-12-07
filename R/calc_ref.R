#' Per recruit calculation and spawning potential ratio
#'
#' \code{calc_ref} Calculates derived spawning potential ratio: lifetime total egg production in fished:unfished states
#'
#' @param Mat_a from report file / true file for simulation
#' @param W_a from report file / true file for simulation
#' @param M from report file / true file for simulation
#' @param S_a from report file / true file for simulation
#' @param F typically terminal estimated/true F, can be any year
#' @param ref FALSE outputs SPR, ref= a value between 0 and 1 can be used with uniroot to find the F at which SPR=ref

#' @return List, a tagged list of potentially useful benchmarks
#' @details Use this function with uniroot to find the value of F that results in SPR matching the specified reference value (e.g. 0.30 to find F30)
#' @export
calc_ref <- function(Mat_a, W_a, M, S_a, F, ref=FALSE){

    ## calculate spawning biomass per recruit in fished and unfished condition
    ## a function of specified level of fishing mortality and ability to estimate selectivity parameters
        Na0 <- Naf <- rep(NA, length(W_a))
        Na0[1] <- Naf[1] <- 1
        for(a in 2:length(W_a)){
            if(a<length(W_a)){
                Na0[a] <- Na0[a-1]*exp(-M)
                Naf[a] <- Naf[a-1]*exp(-M-S_a[a-1]*F)
            }
            if(a==length(W_a)){
                Na0[a] <- (Na0[a-1]*exp(-M))/(1-exp(-M))
                Naf[a] <- (Naf[a-1]*exp(-M-F*S_a[a]))/(1-exp(-M-F*S_a[a]))
            }
        }

        ## ignore recruits
        SB0 <- sum(Na0*Mat_a*W_a)
        SBf <- sum(Naf*Mat_a*W_a)

        ## automatically returns SPR
        SPR <- SBf/SB0
        if(ref==FALSE) return(SPR)
            
        ## can use uniroot function on call to calc_ref to calculate the fishing mortality rate that results in a specified SPR, then compare current fishing mortality with this reference point
        if(ref!=FALSE){
            diff <- ref - SPR
            return(diff)
        }
}