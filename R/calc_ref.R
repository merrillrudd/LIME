#' Per recruit calculation and spawning potential ratio
#'
#' \code{calc_ref} Calculates derived spawning potential ratio: lifetime total egg production in fished:unfished states
#'
#' @author M.B. Rudd
#' @param ages vector of ages
#' @param Mat_a from report file / true file for simulation
#' @param W_a from report file / true file for simulation
#' @param M from report file / true file for simulation
#' @param F proposed F to calculate reference point
#' @param S_fa selectivity by fleet by age
#' @param type "SPR" or "Biomass"
#' @param ref FALSE outputs SPR, ref= a value between 0 and 1 can be used with uniroot to find the F at which SPR=ref

#' @return List, a tagged list of potentially useful benchmarks
#' @details Use this function with uniroot to find the value of F that results in SPR matching the specified reference value (e.g. 0.30 to find F30)
#' @export
calc_ref <- function(ages, Mat_a, W_a, M, F, S_fa, ref=FALSE, type="SBPR"){

    ## calculate spawning biomass per recruit in fished and unfished condition
    ## a function of specified level of fishing mortality and ability to estimate selectivity parameters
    Na0 <- calc_equil_abund(ages=ages, M=M, F=0, S_fa=S_fa, R0=1)
    Naf <- calc_equil_abund(ages=ages, M=M, F=F, S_fa=S_fa, R0=1)

        ## ignore recruits
    if(type=="SBPR"){
        Nofish <- sum(Na0*Mat_a*W_a)
        Fish <- sum(Naf*Mat_a*W_a)
    }
    if(type=="SPR"){
        Nofish <- sum(Na0*Mat_a*W_a*Mat_a*W_a)
        Fish <- sum(Naf*Mat_a*W_a*Mat_a*W_a)
    }
    if(type=="Biomass" | type=="biomass"){
        Nofish <- sum(Na0*W_a)
        Fish <- sum(Naf*W_a)
    }

        ## automatically returns SPR
        ratio <- Fish/Nofish
        if(ref==FALSE) return(ratio)
            
        ## can use uniroot function on call to calc_ref to calculate the fishing mortality rate that results in a specified ratio, then compare current fishing mortality with this reference point
        if(ref!=FALSE){
            diff <- ref - ratio
            return(diff)
        }
}