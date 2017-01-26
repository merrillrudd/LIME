#' Copy simulated data from one data availability scenario to another
#'
#' \code{copy_sim} If two scenarios have the same population dynamics but different data availability going into the estimation model, will copy simulated population from one folder to another with the correct number of years of length composition data
#'
#' @param fromdir single directory from which to copy data
#' @param fromcombos data frame indicating the details of that single directory from which the data will be copied
#' @param todir list of possible directories to copy the data to. the function will filter and copy to the appropriate directories with the same sample size, life history, fishing mortality, and recruitment dynamics
#' @param itervec iterations of simulated data
#' @param rewrite if TRUE will overwrite the files, if FALSE will not overwrite the files if one already exists.

#' @return saves files in new directory
#' @export
copy_sim <- function(fromdir, fromcombos, todir, itervec, rewrite){

	subto <- todir[grepl(fromcombos$ESS,todir) & grepl(fromcombos$LH,todir) & grepl(fromcombos$Fdyn,todir) & grepl(fromcombos$Rdyn,todir)]
	data_avail <- sapply(1:length(subto), function(x) strsplit(strsplit(subto[x],"equil/")[[1]][2],"/")[[1]][1])

	for(sdir in 1:length(subto)){
		for(iter in itervec){

			if(rewrite==FALSE & file.exists(file.path(subto[sdir],iter,"True.rds"))) next
			file <- readRDS(file.path(fromdir, iter, "True.rds"))
			lc <- file$LF
			if(grepl("LBSPR",data_avail[sdir])==FALSE) avail <- as.numeric(strsplit(data_avail[sdir],"LC")[[1]][2])
			if(grepl("LBSPR",data_avail[sdir])) avail <- as.numeric(strsplit(data_avail[sdir],"LBSPR")[[1]][2])
			index <- (nrow(lc)-avail+1):nrow(lc)
			lc_new <- lc[index,]		
			if(length(index)==1){
				lc_new <- t(as.matrix(lc_new))
				rownames(lc_new) <- nrow(lc)
			}

			lc0_new <- file$LF0[index,]		
			if(length(index)==1){
				lc0_new <- t(as.matrix(lc0_new))
				rownames(lc0_new) <- nrow(lc)
			}

			obs_per_year <- rep(0,file$Nyears)
			obs_per_year[index] <- file$obs_per_year[index]		

			file_new <- file
			file_new$LF <- lc_new
			file_new$LF0 <- lc0_new
			file_new$obs_per_year <- obs_per_year
			file_new$DataScenario <- data_avail[sdir]

			saveRDS(file_new, file.path(subto[sdir], iter, "True.rds"))
		}
	}


}