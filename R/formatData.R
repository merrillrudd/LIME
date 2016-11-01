#' Generic function to format required LIME data
#'
#' \code{formatData} Arguments specify the data and assumed values required for LIME, to get everything in one list
#'
#' @param lfreq matrix of length frequency data (number of fish in each length bin), with years along rows and bins along columns
#' @param lfreq_years vector of years for length frequency data
#' @param lbins length bins for length frequency data
#' @param obs_per_year vector wth sample size of length measurements each year
#' @param index vector of abundance index annually, if available. default=NULL
#' @param index_years vector of years for abundance index, if available. default=NULL
#' @param catch vector of catch annually, if available. default=NULL
#' @param catch_years vector of years for catch data, if available. default=NULL
#' @param meanlen vector of mean length annually, if available instead of length composition data. default=NULL
#' @param meanlen_years vector of years for mean length data, if available instead of length composition data. default=NULL
#' @param model_years number of years to be modeled - must be at least the continuous years for which data is available. can be longer. must be specified

#' @return List, a tagged list of data for input into model
#' @details required output: I_t: index with each element named by year 1-x, C_t: catch with each element named 1-x, LF: length frequency with years 1-x labeled on the rows and length bin labeled on the columns, LFprop: proportions in each length bin, same dimensions as LF, years: actual years of data, years_i: index years 1-x, lbins: length bins, ML_t: mean length with each element named year 1-x, Nyears: number of years, Nyears_comp: number of years o f length composition data, obs_per_yr: effective sample size of length composition annually
#' @export
formatData <- function(lfreq, lfreq_years, lbins,  obs_per_year, index=NULL, index_years=NULL, catch=NULL, catch_years=NULL, meanlen=NULL, meanlen_years=NULL, model_years){

    LF <- lfreq
    rownames(LF) <- which(model_years %in% lfreq_years)
    colnames(LF) <- lbins

    LFprop <- LF/rowSums(LF)
    rownames(LFprop) <- which(model_years %in% lfreq_years)
    colnames(LFprop) <- lbins

    names(obs_per_year) <- which(model_years %in% lfreq_years)

    I_t <- index
    if(is.null(index)==FALSE){
        names(I_t) <- which(model_years %in% index_years)
    }

    C_t <- catch
    if(is.null(catch)==FALSE){
        names(C_t) <- which(model_years %in% catch_years)
    }

    ML_t <- meanlen
    if(is.null(meanlen)==FALSE){
        names(ML_t) <- which(model_years %in% meanlen_years)        
    }

    years_o <- unique(c(lfreq_years, index_years, catch_years, meanlen_years))[order(unique(c(lfreq_years, index_years, catch_years, meanlen_years)))]
    years_i <- which(model_years %in% years_o)

    DataList <- NULL
    DataList$I_t <- I_t
    DataList$C_t <- C_t
    DataList$LF <- LF
    DataList$LFprop <- LFprop
    DataList$years_o <- years_o
    DataList$years_i <- years_i
    DataList$years_t <- model_years
    DataList$lbins <- lbins
    DataList$ML_t <- ML_t
    DataList$Nyears <- length(model_years)
    DataList$Nyears_comp <- nrow(LF)
    DataList$obs_per_year <- obs_per_year
    return(DataList)
}
