#' Set up model paths
#'
#' \code{model_paths} sets up directory paths for model combinations
#'
#' @param res_dir directory to store results
#' @param modcombos data frame of all model combinations

#' @return vector of directory paths
#' @export
model_paths <- function(res_dir, modcombos){

    devo_path <- function(combo, res_dir){
        old <- res_dir
        for(i in 1:length(combo)){
            new <- file.path(old, combo[i])
            old <- new
            dir.create(old, showWarnings=FALSE)
        }
        return(old)
    }

    alldirs <- sapply(1:nrow(modcombos), function(x) devo_path(combo=modcombos[x,], res_dir=res_dir))

    return(alldirs)
}