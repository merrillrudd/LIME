#' Simulate spatial variation in growth
#'
#' \code{spatialgrowth_sim} simulate 1D Spatial variation in asymptotic length across sites
#'
#' @param n_i number of sites
#' @param Scale Gaussian scale parameter, default 2
#' @param Sigma2 variance in asymptotic length across sites, default 1
#' @param SD_spatial Gaussian variation (standard deviation), default 0.1
#' @param linf average asymptotic length across sites
#' @param beta_y trend in asymptotic length over sites, default 0.02

#' @return data frame of asymptotic length at each site
#' @export
spatialgrowth_sim <- function(n_i, Scale=2, Sigma2=1, SD_spatial=0.1, linf, beta_y=0.02){
    # require(RandomFields)

    ## sample locations
    lat_min <- -4
    lat_max <- -1
    y_i <- runif(n=n_i, min=lat_min, max=lat_max)

    ## simulate spatial process
    RMmodel <- RMgauss(var=SD_spatial^2, scale=Scale)
    linf_i1 <- linf * exp(RFsimulate(model=RMmodel, x=rep(0, n_i), y=y_i)@data[,1] - Sigma2/2) * exp( beta_y*(y_i - mean(c(lat_min, lat_max))))
    linf_i <- linf_i1 + (linf-mean(linf_i1))

    df <- data.frame(linf_i=linf_i, y_i=y_i)
    return(df)
}