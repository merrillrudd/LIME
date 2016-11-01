#' Settings for simulation testing of data availability scenarios
#'
#' \code{data_avail_settings} Specifies 4 distinct life history types for simulation testing
#'
#' @param avail_set names of data availability scenarios, must match names in function, see details
#' @param ESS Effective sample size, default 1000
#' @param simulation default=TRUE, not set up for simulation=FALSE 

#' @return List, a tagged list of model settings
#' @details avail_set: "Rich_LC"= 20 years length comp, index, and catch, high ESS; "Moderate_LC"=20 years length comp, index, and catch, lower ESS;  "Sample_LC"=20 years length comp, index, and catch, but only sampled every few years, and catch is reported at a rate of 20 percent; "Index_LC1"=20 years of abundance index, no catch data, 1 year of Length comp (can also specify "Index_LC10" for 10 years of length comp); "Catch_LC10"=20 years of catch data, 1 year of length comp (can also specify "Catch_LC10" for 10 years of length comp); "LC1" is 1 year of length comp only; can also specify LC2, 5, 10, and 20.
#' @export
data_avail_settings <- function(avail_set, ESS, simulation=TRUE){
        
    ### simulation
    if(simulation==TRUE){
        settings <- list()
        Nyears <- 20
        if("Rich_LC" %in% avail_set){
            settings$Rich_LC <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=Nyears, alt_yrs=FALSE, sample=FALSE)
        }
        if("Moderate_LC" %in% avail_set){
            settings$Moderate_LC <- list(Nyears=Nyears, comp_sample=50, obs_per_yr=rep(ESS,Nyears), Nyears_comp=Nyears, alt_yrs=FALSE, sample=FALSE)
        }
        if("Sample_LC" %in% avail_set){
            settings$Sample_LC <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=Nyears, alt_yrs=TRUE, sample=0.2)
        }
        if("Index_LC1" %in% avail_set){
            settings$Index_LC1 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=1, alt_yrs=FALSE, sample=FALSE)
        }
        if("Catch_LC1" %in% avail_set){
            settings$Catch_LC1 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=1, alt_yrs=FALSE, sample=FALSE)
        }
        if("Index_LC10" %in% avail_set){
            settings$Index_LC10 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=10, alt_yrs=FALSE, sample=FALSE)
        }
        if("Catch_LC10" %in% avail_set){
            settings$Catch_LC10 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=10, alt_yrs=FALSE, sample=FALSE)
        }

        ### must have at least 1 year of composition data specified to test method with 0 years length comp - just don't include it in data supplied to model
        if("Index_LC0" %in% avail_set){
            settings$Index_LC0 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=1, alt_yrs=FALSE, sample=FALSE)
        }
        if("Catch_LC0" %in% avail_set){
            settings$Catch_LC0 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=1, alt_yrs=FALSE, sample=FALSE)
        }

        if("LC1" %in% avail_set){
            settings$LC1 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=1, alt_yrs=FALSE, sample=FALSE)
        }
        if("LC2" %in% avail_set){
            settings$LC2 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=2, alt_yrs=FALSE, sample=FALSE)
        }
        if("LC5" %in% avail_set){
            settings$LC5 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=5, alt_yrs=FALSE, sample=FALSE)
        }
        if("LC10" %in% avail_set){
            settings$LC10 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=10, alt_yrs=FALSE, sample=FALSE)
        }
        if("LC20" %in% avail_set){
            settings$LC20 <- list(Nyears=Nyears, comp_sample=1000, obs_per_yr=rep(ESS,Nyears), Nyears_comp=20, alt_yrs=FALSE, sample=FALSE)
        }
    }

    return(settings)

}
