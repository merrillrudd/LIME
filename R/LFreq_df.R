#' transforms length composition data from matrix, array, or list to data frame of observations
#'
#' \code{LFreq_df} transforms length frequency data in matrix, array, or list to data frame of observations
#'
#' @author M.B. Rudd
#' @param LF length composition data with years along the rows, length bins along the columns, with fleets as the 3rd dimension in an array or as separate elements of a list
#' @return data frame with Fleet, Year, and Length observation
#' 
#' @export
LFreq_df <- function(LF){

	if((is.array(LF) | is.list(LF)) == FALSE) stop("This function converts length composition matrix to a data frame. LF is not a matrix, array, or list.")
	if(is.array(LF)){
		years <- as.numeric(rownames(LF))
		bins <- as.numeric(colnames(LF))
		if(is.matrix(LF)) nfleets <- 1
		if(is.array(LF) & length(dim(LF))>2) nfleets <- dim(LF)[3]

		byFleet <- lapply(1:nfleets, function(x){
			if(is.matrix(LF)){
				sub <- LF
			}
			if(is.array(LF) & length(dim(LF))>2){
				sub <- LF[,,x]
			}
			byYear <- lapply(1:length(years), function(y){
				sub2 <- sub[y,]
				lengths <- unlist(sapply(1:length(bins), function(z) rep(bins[z],sub2[z])))
				df <- data.frame("Fleet"=x, "Year"=years[y], "Length"=lengths)
				return(df)
			})
			byYear <- do.call(rbind, byYear)
			return(byYear)
		})
		byFleet <- do.call(rbind, byFleet)
	}
	if(is.list(LF)){
		years <- as.numeric(rownames(LF[[1]]))
		bins <- as.numeric(colnames(LF[[1]]))
		nfleets <- length(LF)

		byFleet <- lapply(1:nfleets, function(x){
			sub <- LF[[x]]
			byYear <- lapply(1:length(years), function(y){
				sub2 <- sub[y,]
				lengths <- unlist(sapply(1:length(bins), function(z) rep(bins[z],sub2[z])))
				df <- data.frame("Fleet"=factor(x), "Year"=factor(years[y]), "Length"=lengths)
				return(df)
			})
			byYear <- do.call(rbind, byYear)
			return(byYear)
		})
		byFleet <- do.call(rbind, byFleet)
	}

	return(byFleet)
}