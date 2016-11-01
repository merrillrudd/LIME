#' Calculate AIC
#'
#' \code{calc_AIC} calculate AIC for model runs
#'
#' @param modpath_vec path to model run

#' @return matrix with AIC, AICc, deltaAIC, and deltaAICc (in that order by column) for each model (down the rows)
#' @export
calc_AIC <- function(modpath_vec){
    aic_mat <- matrix(NA, nrow=length(modpath_vec), ncol=6)
    colnames(aic_mat) <- c("AIC", "AICc", "deltaAIC", "deltaAICc", "relLikeAIC", "relLikeAICc")

    for(i in 1:length(modpath_vec)){
        input <- readRDS(file.path(modpath_vec[i], "Inputs2.rds"))
        report <- readRDS(file.path(modpath_vec[i], "Report.rds"))

        nll <- report$jnll
        params <- length(as.vector(unlist(input$Parameters))) - length(which(grepl("Nu_input", names(unlist(input$Parameters)))))
        sampsize <- length(which(input$Data$I_t>0)) + sum(input$Data$obs_per_yr) + length(which(input$Data$C_t>0))

        AIC <- 2*params + 2*nll
        AICc <- AIC + (2*params*(params + 1))/(sampsize - params - 1)

        aic_mat[i,1] <- AIC
        aic_mat[i,2] <- AICc
    }

    aic_mat[,"deltaAIC"] <- aic_mat[,"AIC"] - min(aic_mat[,"AIC"])
    aic_mat[,"deltaAICc"] <- aic_mat[,"AICc"] - min(aic_mat[,"AICc"])
    aic_mat[,"relLikeAIC"] <- sapply(1:nrow(aic_mat), function(x) exp((min(aic_mat[,"AIC"]) -aic_mat[x,"AIC"])/2))
    aic_mat[,"relLikeAICc"] <- sapply(1:nrow(aic_mat), function(x) exp((min(aic_mat[,"AICc"]) - aic_mat[x,"AICc"])/2))

    return(aic_mat)
}