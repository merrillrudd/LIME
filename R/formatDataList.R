#' Formatting real data for specific species for TMB model
#'
#' \code{formatDataList} Reads data files from directory for specific species and gets the data into a form to input to the TMB model
#'
#' @param species species code, "SIGSUT" or "CRSNAP"
#' @param data_dir Directory to look for data files 

#' @return List, a tagged list of data to input into TMB model
#' @details required output: I_t: index with each element named by year 1-x, C_t: catch with each element named 1-x, LF: length frequency with years 1-x labeled on the rows and length bin labeled on the columns, LFprop: proportions in each length bin, same dimensions as LF, years: actual years of data, years_i: index years 1-x, lbins: length bins, ML_t: mean length with each element named year 1-x, Nyears: number of years, Nyears_comp: number of years o f length composition data, obs_per_yr: effective sample size of length composition annually
#' @export
formatDataList <- function(species, data_dir){

    if(species=="SIGSUT"){
        sp <- "Siganus_sutor"
        data <- read.csv(file.path(data_dir, paste0(sp, "_data.csv")), header=TRUE)
        meltDat <- melt(data, id.vars=c("Date", "Year", "Landing", "Management.type", "Group", 
            "Fishgear","X.Fishers", "Catch.category", "Family", "Name", "Catch.category",
            "Commercial.noncomm."))
        meltDat$value <- as.numeric(meltDat$value)
        dat_new <- meltDat[-which(is.na(meltDat$value)),-which(colnames(meltDat)=="variable")]      
    
        years <- unique(dat_new$Year)[order(unique(dat_new$Year))]
        tyears <- min(years):max(years)
        tyears_i <- 1:length(tyears)
        years_i <- tyears_i[which(tyears %in% years)]
        lbins <- 1:(max(dat_new$value)*1.5)     

        lcomp1 <- matrix(0, nrow=length(years), ncol=length(lbins))
        colnames(lcomp1) <- lbins
        rownames(lcomp1) <- years_i   
        lcomp_p <- lcomp1   

        catch <- cpue <- vector(length=length(years))
        cpue_day_out <- list()      

        obs_high <- obs_low <- rep(0, length(years))
        for(y in 1:length(years)){
            sub <- dat_new[which(dat_new$Year==years[y]),]   
            dates <- unique(as.character(sub$Date))  
            obs_high[y] <- length(dates)
            months <- unique(as.numeric(sapply(1:length(dates), function(x) strsplit(as.character(dates[x]), "/")[[1]][2])))
            obs_low[y] <- length(months)

            ## annual catch
            catch[y] <- nrow(sub)      

            ## annual cpue - average of daily/group cpue
            cpue_sub <- rep(NA, length(unique(sub$Date)))
            ngroups <- rep(NA, length(unique(sub$Date)))
            for(d in 1:length(unique(sub$Date))){
                subd <- sub[which(sub$Date==unique(sub$Date)[d]),]    
                ngroups <- length(unique(subd$Group)[which(is.na(unique(subd$Group))==FALSE)])
                if(ngroups==0) cpue_sub[d] <- NA
                if(ngroups!=0) cpue_sub[d] <- nrow(subd)/ngroups
                # rm(ngroups)
                # rm(subd)
            }
            if(all(is.na(cpue_sub))) cpue[y] <- NA
            if(all(is.na(cpue_sub)==FALSE)) cpue[y] <- mean(cpue_sub, na.rm=TRUE)
            cpue_day_out[[y]] <- cpue_sub
            rm(cpue_sub)       

            for(l in 1:length(lbins)){
                lcomp1[y,l] <- length(which(floor(sub$value)==lbins[l]))
            }
        }       

        names(cpue) <- names(catch) <- years_i
        cpue <- cpue[-c(which(is.na(cpue)),which(cpue==0))]     

        meanlen <- rep(NA, nrow(lcomp1))
        names(meanlen) <- rownames(lcomp1)
        for(i in 1:nrow(lcomp1)){
            lcomp_p[i,] <- lcomp1[i,]/sum(lcomp1[i,])
            meanlen[i] <- sum(sapply(1:ncol(lcomp1), function(x) lcomp1[i,x]*as.numeric(colnames(lcomp1)[x])))/sum(lcomp1[i,])
        }  

        DataList <- NULL
        DataList$I_t <- cpue
        DataList$C_t <- catch
        DataList$LF <- lcomp1
        DataList$LFprop <- lcomp_p
        DataList$years <- years
        DataList$years_i <- years_i
        DataList$lbins <- lbins
        DataList$meanlen <- meanlen
        DataList$Nyears <- length(min(years):max(years))
        DataList$Nyears_comp <- nrow(LFprop)
        DataList$obs_per_year <- obs_high
    }

    if(species=="CRSNAP"){

        lg <<- read.csv(file.path(data_dir, "cr_snapper_filtered.csv"))

        ## subset by gears
        lg_bl <- lg[which(lg$Gear=="Bottom Longline"),]
        lg_g <- lg[which(lg$Gear=="Gillnet"),]      

        ## annual mean length
        ml <- mean_length(data=lg, plot=FALSE)      

        ## life history info
        cr_lh <- choose_lh_list(species="CRSNAP", selex="asymptotic")       

        ## length frequency data
        lf <- length_frequency(binwidth=1, linf=cr_lh$linf, lmat=cr_lh$L50, data=lg, plot=FALSE, weight=TRUE)       

        ## catch and effort data
        fishery_data <- catch_effort(data=lg, sep_index=TRUE)
        catch <- fishery_data$catch
        cpue_bl <- fishery_data$cpue_bl
        cpue_g <- fishery_data$cpue_g
        tyears <- c((fishery_data$years[1]-10):(fishery_data$years[1]-1),fishery_data$years)
        tyears_i <- 1:length(tyears)
        years <- fishery_data$years
        years_i <- tyears_i[which(tyears %in% years)]
        obs_high <- obs_low <- rep(0, length(tyears_i))
        obs_high[which(tyears_i %in% years_i)] <- fishery_data$obs_high ## number of days fished per year 
        obs_low[which(tyears_i %in% years_i)] <- fishery_data$obs_low ## number of months fished per year

        ## bottom longline cpue index
        cpue_input <- cpue_bl[which(is.na(cpue_bl)==FALSE)] ## choose bottom longline cpue only
        cpue_yrs <- which(tyears %in% names(cpue_bl)[which(is.na(cpue_bl)==FALSE)])
        names(cpue_input) <- cpue_yrs       

        ## length frequency
        lf_input <- lf[which(rowSums(lf)>0),]
        lf_yrs_i <- tyears_i[which(tyears %in% names(which(rowSums(lf)>0)))]
        rownames(lf_input) <- lf_yrs_i      

        lcomp_p <- t(apply(lf_input, 1, function(x) x/sum(x)))
        rownames(lf_input) <- lf_yrs_i      

        lbins <- 1:ncol(lcomp_p)        

        ## mean length
        meanlen_input <- ml$all_gears
        meanlen_yrs_i <- tyears_i[which(tyears %in% names(meanlen_input))]
        names(meanlen_input) <- meanlen_yrs_i

                DataList <- NULL
                DataList$I_t <- cpue_input
                DataList$C_t <- NULL
                DataList$LF <- lf_input
                DataList$LFprop <- lcomp_p
                DataList$years <- tyears
                DataList$years_i <- tyears_i
                DataList$lbins <- lbins
                DataList$ML_t <- meanlen_input
                DataList$Nyears <- length(tyears_i)
                DataList$Nyears_comp <- nrow(lcomp_p)
                DataList$obs_per_year <- obs_high


    }
    return(DataList)

}
