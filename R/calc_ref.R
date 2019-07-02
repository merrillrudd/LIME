#' Per recruit calculation and spawning potential ratio
#'
#' \code{calc_ref} Calculates derived spawning potential ratio: lifetime total egg production in fished:unfished states
#'
#' @author M.B. Rudd
#' @param ages vector of ages
#' @param Mat_a from report file / true file for simulation
#' @param W_a from report file / true file for simulation
#' @param M from report file / true file for simulation
#' @param F proposed F to calculate reference point, if S_fa has more than 1 fleet, requires F as a vector, one for each fleet
#' @param S_fa selectivity by fleet by age
#' @param ref FALSE outputs SPR, ref= a value between 0 and 1 can be used with uniroot to find the F at which SPR=ref
#' @param fleet_prop fleet proportions

#' @return List, a tagged list of potentially useful benchmarks
#' @details Use this function with uniroot to find the value of F that results in SPR matching the specified reference value (e.g. 0.30 to find F30)
#' @export
calc_ref <- function(ages, Mat_a, W_a, M, F, S_fa, ref=FALSE, fleet_prop=1){

    ## calculate spawning biomass per recruit in fished and unfished condition
    ## a function of specified level of fishing mortality and ability to estimate selectivity parameters

    F_inp <- sapply(1:nrow(S_fa), function(x) F*fleet_prop[x])
    # if(ref!=FALSE) F_inp <- F
    Na0 <- calc_equil_abund(ages=ages, M=M, F=rep(0,length(fleet_prop)), S_fa=S_fa, R0=1)
    Naf <- calc_equil_abund(ages=ages, M=M, F=F_inp, S_fa=S_fa, R0=1)

        ## ignore recruits
    # really sbpr
    # if(type=="SPR"){ 
        Nofish <- sum(Na0*Mat_a*W_a)
        Fish <- sum(Naf*Mat_a*W_a)
    # }
    # really spr
    # if(type=="SPR2"){ 
        # Nofish <- sum(Na0*Mat_a*W_a*Mat_a*W_a)
        # Fish <- sum(Naf*Mat_a*W_a*Mat_a*W_a)
    # }
    # if(type=="Biomass" | type=="biomass"){
    #     Nofish <- sum(Na0*W_a)
    #     Fish <- sum(Naf*W_a)
    # }

        ## automatically returns SPR
        ratio <- Fish/Nofish
        if(ref==FALSE) return(ratio)
            
        ## can use uniroot function on call to calc_ref to calculate the fishing mortality rate that results in a specified ratio, then compare current fishing mortality with this reference point
        if(ref!=FALSE){
            diff <- ref - ratio
            return(diff)
        }
}