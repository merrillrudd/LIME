#' Operating model
#'
#' \code{sim_pop} Age-converted-to-length-based operating model specifying true population dynamics
#'
#' @author M.B. Rudd
#' @param lh list of life history information, from create_lh_list
#' @param Nyears number of years to simulate
#' @param pool if nseasons (in life history list) is greater than one, pool the generated data into annual time steps, or leave at the season-level? FALSE will generate shorter time step life history info, mean length
#' @param Fdynamics Specify name of pattern of fishing mortality dynamics, Constant, Endogenous, Ramp, Increasing, or None. Input number to project forward using a specific F.
#' @param Rdynamics Specify name of pattern of recruitment dynamics, Constant, Pulsed, Pulsed_up, or BH
#' @param Nyears_comp number of years of length composition data
#' @param comp_sample sample size of length composition data annually
#' @param init_depl initial depletion on which to calculate F1; default = 0.99
#' @param seed set seed for generating stochastic time series
#' @param sample_type a character vector specifying if the length comps are sampled from the 'catch' (default) or from the population
#' @param mgt_type removals based on F (default) or catch
#' @param fleet_proportions vector specifying the relative size of each fleet in terms of fishing pressure. must have length = nfleets and sum to 1.
#' @importFrom stats rnorm
#' @return named list of attributes of true population/data
#' @export
sim_pop <-
  function(lh,
           Nyears,
           pool,
           Fdynamics,
           Rdynamics,
           Nyears_comp,
           comp_sample,
           init_depl=0.99,
           seed,
           sample_type = 'catch',
           mgt_type = 'F',
           fleet_proportions) {

    with(lh, {
      ##########################
      ## Initial calcs
      ##########################

      Nyears_real <- Nyears

      Nyears <- Nyears * nseasons
      tyears <- Nyears_real * nseasons
      Nyears_comp <- Nyears_comp * nseasons


      ##########################
      ## Random variables
      ##########################
      set.seed(seed)
    
      ## recruitment deviations
      RecDev <- c(0,rnorm(Nyears_real-1, mean = -(SigmaR ^ 2) / 2, sd = SigmaR))
    
      ## autocorrelated recruitment deviations
      RecDev_AR <- rep(NA, length(RecDev))
      RecDev_AR[1] <- RecDev[1]
      for (t in 2:length(RecDev)) {
        RecDev_AR[t] <-
          RecDev_AR[t - 1] * rho + sqrt(1 - rho ^ 2) * RecDev[t]
      }
      if(nseasons > 1){
        # RecDev2 <- rep(0, tyears)
        # RecDev2[seq(1,length(RecDev2),by=nseasons)] <- RecDev_AR
        RecDev2 <- unlist(lapply(1:Nyears_real, function(x) rep(RecDev_AR[x], nseasons)))
        RecDev_AR <- RecDev2
      }

      ## fishing mortality deviations
      if(length(SigmaF)==1 & nfleets>1) SigmaF <- rep(SigmaF, nfleets)
      FishDev_f <- t(sapply(1:nfleets, function(x){
        rnorm(tyears, mean = -(SigmaF[x] ^ 2) / 2, sd = SigmaF[x])
      }))
    
      ## abundance index observation error
      if(length(SigmaI)==1 & nfleets>1) SigmaI <- rep(SigmaI, nfleets)
      IndexDev_f <- t(sapply(1:nfleets, function(x){
        rnorm(tyears, mean = -(SigmaI[x] ^ 2) / 2, sd = SigmaI[x])
      }))

      ## catch observation error
      if(length(SigmaC)==1 & nfleets>1) SigmaC <- rep(SigmaC, nfleets)
      CatchDev_f <- t(sapply(1:nfleets, function(x){
        rnorm(tyears, mean = -(SigmaC[x] ^ 2) / 2, sd = SigmaC[x])
      }))

      ## catch observation error
      if(length(SigmaC)==1 & nfleets>1) SigmaC <- rep(SigmaC, nfleets)
      CatchDev_f <- sapply(1:nfleets, function(x){
        rnorm(tyears, mean = -(SigmaC[x] ^ 2) / 2, sd = SigmaC[x])
      })



      #########################
      ## Setup recruitment
      #########################
      if (Rdynamics == "Pulsed"){
        choose <- sample(1:2, 1)
        R_t <- c(
          "initial" = rep(R0, floor(tyears / 3)),
          "pulse1" = rep(ifelse(choose==1, (R0 / 2), (R0 * 2)), floor(tyears / 3)),
          "pulse2" = rep(ifelse(choose==1, (R0 *2), (R0 / 2)), tyears - (2 * floor(tyears / 3)))
        ) * exp(RecDev_AR)
        R_t[which(RecDev_AR==0)] <- 0
      }

      if (Rdynamics == "Constant") {
        R_t <- (rep(R0, tyears)) * exp(RecDev_AR)
        R_t[which(RecDev_AR==0)] <- 0
      }

    ## with multi-seasons, currently spreading out recruitment across seasons in one year
    ## to properly calculate initial F if recruitment is only in one season, need to spread out calculation across the year instead of per-recruit in calc_ref
        R_t <- R_t/nseasons


      #####################################
      ## Effort dynamics
      #####################################
      ## reference points
      F40 <-
        tryCatch(
          uniroot(
            calc_ref,
            lower = 0,
            upper = 200,
            ages = ages,
            Mat_a = Mat_a,
            W_a = W_a,
            M = M,
            ref = 0.4,
            type="SPR"
          )$root,
          error = function(e)
            NA
        )
      Finit <-
        tryCatch(
          uniroot(
            calc_ref,
            lower = 0,
            upper = 200,
            ages = ages,
            Mat_a = Mat_a,
            W_a = W_a,
            M = M, 
            ref = init_depl,
            type = "biomass"
          )$root,
          error = function(e)
            NA
        )
      if (is.na(Finit))
        stop("F corresponding to initial depletion does not exist")
      Fmax <-
        tryCatch(
          uniroot(
            calc_ref,
            lower = 0,
            upper = 200,
            ages = ages,
            Mat_a = Mat_a,
            W_a = W_a,
            M = M,
            ref = 0.05,
            type = "biomass"
          )$root,
          error = function(e)
            NA
        )
      if (is.na(Fmax) | Fmax > 3)
        Fmax <- 3

      if(length(Fdynamics)==1 & nfleets>1) Fdynamics <- rep(Fdynamics, nfleets)
      
      ## fishing mortality dynamics
      F_ft <- matrix(NA, nrow=nfleets, ncol=tyears)
      if(any(Fdynamics=="Constant")){
        index <- which(Fdynamics=="Constant")
        for(i in 1:length(index)){
          F_ft[index[i],] <- Finit * exp(FishDev_f[index[i],])
        }
      }
      if(any(Fdynamics=="Oneway")){
        index <- which(Fdynamics=="Oneway")
        for(i in 1:length(index)){
          F_ft[index[i],] <- c(seq(1,by=0.05,length=Nyears)) * Finit * exp(FishDev_f[index[i],])
        }
      }
      if(any(Fdynamics=="Endogenous")){
        index <- which(Fdynamics=="Endogenous")
        for(i in 1:length(index)){
          F_ft[index[i],1] <- Finit
        }
      }
      if (any(Fdynamics == "None")){
        index <- which(Fdynamics=="None")
        for(i in 1:length(index)){
          F_ft[index[i],] <- rep(0, tyears) * exp(FishDev_f[index[i],])
        }
      }

      ## fishing mortality = include selectivity and Finit with effort dynamics and relative weight of fishery to scale each fishery
      F_atf <- array(NA, dim=c(length(ages), tyears, nfleets))
      for(f in 1:nfleets){
        for(t in 1:tyears){
          for(a in 1:length(ages)){
            F_atf[a,t,f] <- F_ft[f,t] * S_fa[f,a]
          }
        }
      }

      F_at <- matrix(NA, nrow=length(ages), ncol=tyears)
      for(t in 1:tyears){
        for(a in 1:length(ages)){
          F_at[a,t] <- sum(F_atf[a,t,])
        }
      }


      ##########################
      ## Data objects
      ##########################

      ## year 1 abundance at age
      N_at <- N_at0 <- matrix(NA, nrow = length(L_a), ncol = tyears)
      for (a in 1:length(L_a)) {
        if (a == 1) {
          N_at[a, 1] <- R_t[1]
          N_at0[a, 1] <- R_t[1]
        }
        if (a > 1 & a < length(L_a)) {
          N_at[a, 1] <- N_at[a - 1, 1] * exp(-M - F_at[a,1])
          N_at0[a, 1] <- N_at0[a - 1, 1] * exp(-M)
        }
        if (a == length(L_a)) {
          N_at[a, 1] <-
            (N_at[a - 1, 1] * exp(-M - F_at[a,1])) / (1 - exp(-M - F_at[a,1]))
          N_at0[a, 1] <- (N_at0[a - 1, 1] * exp(-M)) / (1 - exp(-M))
        }
      }

      ## year 1 biomass quantities
      TB_t <- SB_t <- rep(NA, tyears)
      TB_t[1] <- sum(N_at[, 1] * W_a)
      SB_t[1] <- sum(N_at[, 1] * W_a * Mat_a)

      # if(is.numeric(Fdynamics) & mgt_type=="catch"){
      #   F_t[1] <- max(0.01, getFt(ct=C_t[1], m=M, sa=S_a, wa=W_a, na=N_at[,1]))
      #   F_t[1] <- min(c(Fmax, F_t[1]), na.rm=TRUE)
      # }

      ## year 1 catch
      Cn_atf <- Cw_atf <- array(NA, c(length(L_a), tyears, nfleets))
      for(f in 1:nfleets){
        Cn_atf[,1,f] <- N_at[,1] * (1 - exp(-M - F_atf[,1,f])) * (F_atf[,1,f] / (M + F_atf[a,1,f]))
        Cw_atf[,1,f] <- Cn_atf[,1,f] * W_a
      }

      ## unfished spawning biomass
      SB0 <- sum(R0/nseasons * exp(-M*nseasons * ages) * W_a * Mat_a)      

      ##########################
      ## Projection
      ##########################

      for (t in 2:tyears) {
        ## fishing effort based on spawning biomass
        if (any(Fdynamics == "Endogenous")) {
          index <- which(Fdynamics == "Endogenous")
          for(i in 1:length(index)){
            F_ft[index[i],t] <- (F_ft[index[i],t-1] * (SB_t[t-1] / (Fequil * SB0/2)) ^ Frate) * exp(FishDev_f[index[i],t])
          }
          ## fishing mortality = include selectivity 
          for(i in 1:length(index)){
            for(a in 1:length(ages)){
              F_atf[a,t,index[i]] <- F_ft[index[i],t] * S_fa[index[i],a]
            }            
          }
          for(a in 1:length(ages)){
            F_at[a,t] <- sum(F_atf[a,t,])
          }
        }


        ## age-structured dynamics
        for (a in 1:length(L_a)) {
          if (a == 1) {
            N_at[a, t] <- R_t[t]
            N_at0[a, t] <- R_t[t]
          }
          if (a > 1 & a < length(L_a)) {
            N_at[a, t] <- N_at[a - 1, t - 1] * exp(-M - F_at[a-1,t-1])
            N_at0[a, t] <- N_at0[a - 1, t - 1] * exp(-M)
          }
          if (a == length(L_a)) {
            N_at[a, t] <- (N_at[a - 1, t - 1] * exp(-M - F_at[a-1,t-1])) + (N_at[a, t - 1] * exp(-M - F_at[a,t-1]))
            N_at0[a, t] <- (N_at0[a - 1, t - 1] * exp(-M)) + (N_at0[a, t - 1] * exp(-M))
          }
        }

          ## spawning biomass
          SB_t[t] <- sum((N_at[, t] * W_a * Mat_a))
          TB_t[t] <- sum(N_at[, t] * W_a)

          # if(is.numeric(Fdynamics) & mgt_type=="catch"){
          #   F_t[y] <- max(0.01,getFt(ct=C_t[y], m=M, sa=S_a, wa=W_a, na=N_at[,y]))
          #   # F_t[y] <- min(c(Fmax, F_t[y]), na.rm=TRUE)
          # }

          ## catch
          for(f in 1:nfleets){
            Cn_atf[,t,f] <- N_at[,t] * (1 - exp(-M - F_atf[,t,f])) * (F_atf[,t,f] / (M + F_atf[,t,f]))
            Cw_atf[,t,f] <- Cn_atf[,t,f] * W_a
          }

      }

      F_ft <- t(sapply(1:nfleets, function(x){
        sub <- F_atf[,,x]
        findMax <- sapply(1:tyears, function(y){
          sub2 <- sub[,y]
          return(max(sub2))
        })
        return(findMax)
      }))
      F_t <- colSums(F_ft)


      SPR_t <-
        sapply(1:length(F_t), function(x)
          calc_ref(
            ages = ages,
            Mat_a = Mat_a,
            W_a = W_a,
            M = M,
            F = F_t[x]
          ))
      SPR <- SPR_t[length(SPR_t)]

      Cn_ft <- t(sapply(1:nfleets, function(x) colSums(Cn_atf[,,x])))
      Cw_ft <- t(sapply(1:nfleets, function(x) colSums(Cn_atf[,,x] * W_a)))

          F_fy <- t(sapply(1:nfleets, function(y){
              sapply(1:Nyears_real, function(x) {
                if (nseasons == 1)
                  time_index <- x
                if (nseasons > 1)
                  time_index <- (1:nseasons) + ((x - 1) * nseasons)
                sum(F_ft[y,time_index])
              })
          }))
          F_y <- colSums(F_fy)  

      if(pool==FALSE) N_t <- colSums(N_at[-1,])
      if(pool==TRUE){
          N_t <- colSums(N_at[-1, which(1:tyears %% nseasons == 0)])
          SB_t <- SB_t[which(1:tyears %% nseasons == 0)]
          TB_t <- TB_t[which(1:tyears %% nseasons == 0)]
          SPR_t <- SPR_t[which(1:tyears %% nseasons == 0)]

          Cn_ft <- t(sapply(1:nfleets, function(y){
              sapply(1:Nyears_real, function(x) {
                if (nseasons == 1)
                  time_index <- x
                if (nseasons > 1)
                  time_index <- (1:nseasons) + ((x - 1) * nseasons)
                sum(Cn_ft[y,time_index])
              }) #* exp(CatchDev - (SigmaC^2)/2)
            }))
          Cn_t <- colSums(Cn_ft)
          Cw_ft <- t(sapply(1:nfleets, function(y){
              sapply(1:Nyears_real, function(x) {
                if (nseasons == 1)
                  time_index <- x
                if (nseasons > 1)
                  time_index <- (1:nseasons) + ((x - 1) * nseasons)
                sum(Cw_ft[y,time_index])
              }) #* exp(CatchDev - (SigmaC^2)/2)
          }))
          Cw_t <- colSums(Cw_ft)     

          F_ft <- t(sapply(1:nfleets, function(y){
              sapply(1:Nyears_real, function(x) {
                if (nseasons == 1)
                  time_index <- x
                if (nseasons > 1)
                  time_index <- (1:nseasons) + ((x - 1) * nseasons)
                sum(F_ft[y,time_index])
              })
          }))
          F_t <- colSums(F_ft)  

          R_t <- sapply(1:Nyears_real, function(x) {
            if (nseasons == 1)
              time_index <- x
            if (nseasons > 1)
              time_index <- (1:nseasons) + ((x - 1) * nseasons)
            sum(R_t[time_index])
          })
      }



     

      ## relative spawning biomass (depletion)
      D_t <- SB_t / SB0 

      ## abundance index
      if(length(qcoef)!=nfleets) qcoef <- rep(qcoef, nfleets)
      I_ft <- t(sapply(1:nfleets, function(x) qcoef[x] * TB_t * exp(IndexDev_f[x,])))

      ## age to length comp
      if(length(Nyears_comp)!=nfleets) Nyears_comp <- rep(Nyears_comp, nfleets)

      ## years with observed length comps
      oyears_mat <- matrix(0, nrow=nfleets, ncol=tyears)
      for(f in 1:nfleets){
        oyears <- (tyears-Nyears_comp[f] + 1):tyears
        oyears_mat[f,oyears] <- 1
      }

      obs_per_year <- matrix(0, nrow=nfleets, ncol=tyears)
      for(f in 1:nfleets){
        for(t in 1:tyears){
          if(oyears_mat[f,t]!=0) obs_per_year[f,t] <- (comp_sample/nseasons)
        }
      }

      LFinfo <-lapply(1:nfleets, function(x){
        AgeToLengthComp(
          lh = lh,
          S_a = lh$S_fa[x,],
          tyears = tyears,
          N_at = N_at,
          comp_sample = obs_per_year[x,],
          sample_type = sample_type
        )
      })
      LF0info <- lapply(1:nfleets, function(x){
        AgeToLengthComp(
          lh = lh,
          S_a = lh$S_fa[x,],
          tyears = tyears,
          N_at = N_at0,
          comp_sample = obs_per_year[x,],
          sample_type = sample_type
        )
      })
      plba <- lapply(1:nfleets, function(x) LFinfo[[x]]$plba)
      plb <- lapply(1:nfleets, function(x) LFinfo[[x]]$plb)
      page <- lapply(1:nfleets, function(x) LFinfo[[x]]$page)
      LF <- lapply(1:nfleets, function(x) LFinfo[[x]]$LF)
      LF0 <- lapply(1:nfleets, function(x) LF0info[[x]]$LF)

      LF_tf <- LF0_tf <- NULL
      if (pool == TRUE) {
        for(f in 1:nfleets){
          LF_tf[[f]] <- LF0_tf[[f]] <- matrix(NA, nrow=Nyears_real, ncol=ncol(LF[[1]]))
         for (y in 1:Nyears_real) {
            if (nseasons == 1) {
              LF_tf[[f]][y,] <- LF[[f]][y,]
              LF0_tf[[f]][y,] <- LF0[[f]][y,]
            }
            if (nseasons > 1) {
              time_index <- (1:nseasons) + ((y - 1) * nseasons)
              LF_tf[[f]][y,] <- colSums(LF[[f]][time_index,])
              LF0_tf[[f]][y,] <- colSums(LF0[[f]][time_index,])
            }
          }
          obs_per_year <- t(sapply(1:nfleets, function(y){
            sapply(1:Nyears_real, function(x) {
              if (nseasons == 1)
                time_index <- x
              if (nseasons > 1)
                time_index <- (1:nseasons) + ((x - 1) * nseasons)
              sum(obs_per_year[y,time_index])
            })
          }))
        }
      }
      if (pool == FALSE) {
        LF_tf <- LF
        LF0_tf <- LF0
      }
      for(f in 1:nfleets){
        colnames(LF_tf[[f]]) <- highs
        colnames(LF0_tf[[f]]) <- highs
        if(pool==TRUE){
          rownames(LF_tf[[f]]) <- 1:Nyears_real
          rownames(LF0_tf[[f]]) <- 1:Nyears_real
        }
      }


      ########################################################
      ## Expected mean length in catch
      ########################################################
      ML_ft <- matrix(NA, nrow=nfleets, ncol=tyears)
      for(f in 1:nfleets){
        for (t in 1:tyears) {
          vul_pop <- sum(N_at[, t] * S_fa[f,])
          vul_lengths <- sum(vul_pop * plb[[f]][t,] * mids)
          ML_ft[f,t] <- vul_lengths / vul_pop
        }
      }
      if (pool == TRUE) ML_ft <- t(sapply(1:nfleets, function(x) ML_ft[x,which(1:tyears %% nseasons == 0)]))

      ## generated data
      # I_out <- data.frame("Variable"="Index", "By"="Time", "X"=c(sapply(1:ncol(I_ft), function(x) rep(x, nfleets))), "Value"=c(I_ft), "Fleet"=rep(1:nfleets, ncol(I_ft)))
      # Cn_out <- data.frame("Variable"="Catch_numbers", "By"="Time", "X"=c(sapply(1:ncol(Cn_ft), function(x) rep(x, nfleets))), "Value"=c(Cn_ft), "Fleet"=rep(1:nfleets, ncol(Cn_ft)))
      # Cw_out <- data.frame("Variable"="Catch_biomass", "By"="Time", "X"=c(sapply(1:ncol(Cw_ft), function(x) rep(x, nfleets))), "Value"=c(Cw_ft), "Fleet"=rep(1:nfleets, ncol(Cw_ft)))

      # LFlong <- lapply(1:nfleets, function(z){
      #   LFsub <- LF_tf[[z]]
      #   LFlong <- melt(LFsub)
      #   colnames(LFlong) <- c("X", "LengthBin", "Value")

      #   tyears_vec <- unique(LFlong$X)[order(unique(LFlong$X))]
      #   oyears <- (max(tyears_vec)-Nyears_comp[z] + 1):max(tyears_vec)
      #   LFlonger <- lapply(1:nrow(LFlong), function(y){
      #     if(LFlong$Value[y]>0){
      #       len <- rep(LFlong$LengthBin[y], LFlong$Value[y])
      #       out <- data.frame("By"="Time", "X"=LFlong$X[y], "Variable"="LengthComp","Value"=len)
      #       return(out)
      #     }
      #     if(LFlong$Value[y]==0) return(data.frame("By"="Time", "X"=LFlong$X[y], "Variable"="LengthComp", "Value"=0))
      #   })
      #   LF2 <- do.call(rbind, LFlonger) %>% mutate("Fleet"=z) %>% filter(X %in% oyears)
      #   return(LF2)
      # })
      # LF_out <- do.call(rbind, LFlong) %>% 
      #           filter(Value != 0)

      # LF0long <- lapply(1:nfleets, function(z){
      #   LF0sub <- LF0_tf[[z]]
      #   LF0long <- melt(LF0sub)
      #   colnames(LF0long) <- c("X", "LengthBin", "Value")

      #   tyears_vec <- unique(LF0long$X)[order(unique(LF0long$X))]
      #   oyears <- (max(tyears_vec)-Nyears_comp[z] + 1):max(tyears_vec)

      #   LF0longer <- lapply(1:nrow(LF0long), function(y){
      #     if(LF0long$Value[y]>0){
      #       len <- rep(LF0long$LengthBin[y], LF0long$Value[y])
      #       out <- data.frame("By"="Time", "X"=LF0long$X[y], "Variable"="LengthComp","Value"=len)
      #       return(out)
      #     }
      #     if(LF0long$Value[y]==0) return(data.frame("By"="Time", "X"=LF0long$X[y], "Variable"="LengthComp", "Value"=0))
      #   })
      #   LF02 <- do.call(rbind, LF0longer) %>% mutate("Fleet"=z) %>% filter(X %in% oyears)
      #   return(LF02)
      # })
      # LF0_out <- do.call(rbind, LF0long) %>% 
      #           filter(Value != 0)

      ## population parameters
      # Ff_out <- data.frame("Variable"="F", "By"="Time", "X"=c(sapply(1:ncol(F_ft), function(x) rep(x, nfleets))), "Value"=c(F_ft), "Fleet"=rep(1:nfleets, ncol(F_ft)))
      # ML_out <- data.frame("Variable"="MeanLen", "By"="Time", "X"=c(sapply(1:ncol(ML_ft), function(x) rep(x, nfleets))), "Value"=c(ML_ft), "Fleet"=rep(1:nfleets, ncol(ML_ft)))

      ## not fleet-specific
      # R_out <- data.frame("Variable"="Recruitment", "By"="Time", "X"=1:length(R_t), "Value"=c(R_t), "Fleet"=0)
      # N_out <- data.frame("Variable"="Numbers", "By"="Time", "X"=1:length(N_t), "Value"=c(N_t), "Fleet"=0)
      # SB_out <- data.frame("Variable"="SpawningBiomass", "By"="Time", "X"=1:length(SB_t), "Value"=c(SB_t), "Fleet"=0)
      # D_out <- data.frame("Variable"="RelativeSB", "By"="Time", "X"=1:length(D_t), "Value"=c(D_t), "Fleet"=0)
      # F_out <- data.frame("Variable"="F", "By"="Time", "X"=1:length(F_t), "Value"=c(F_t), "Fleet"=0)
      # Cn_total_out <- data.frame("Variable"="Catch_numbers", "By"="Time", "X"=1:length(Cn_t), "Value"=Cn_t, "Fleet"=0)
      # Cw_total_out <- data.frame("Variable"="Catch_biomass", "By"="Time", "X"=1:length(Cw_t), "Value"=Cw_t, "Fleet"=0)
      # SPR_out <- data.frame("Variable"="SPR", "By"="Time", "X"=1:length(SPR_t), "Value"=c(SPR_t), "Fleet"=0)
      # TB_out <- data.frame("Variable"="TotalBiomass", "By"="Time", "X"=1:length(TB_t), "Value"=c(TB_t), "Fleet"=0)

      # outdf <- rbind(I_out, Cn_out, Cw_out, Cn_total_out, Cw_total_out, LF_out, LF0_out, Ff_out, ML_out, R_out, N_out, SB_out, D_out, F_out, SPR_out, TB_out)
      # outdf$Fleet <- as.factor(outdf$Fleet)


      ## outputs
      # lh$dfsim <- outdf
      lh$plb <- plb
      lh$plba <- plba
      lh$page <- page
      lh$N_at <- N_at
      lh$nlbins <- length(mids)
      if (pool == TRUE) {
        lh$Nyears <- Nyears_real
        lh$years <- 1:Nyears_real
      }
      if (pool == FALSE) {
        lh$Nyears <- Nyears
        lh$years <- 1:Nyears
      }
      lh$obs_per_year <- obs_per_year
      lh$RecDev <- RecDev_AR
      lh$FishDev_f <- FishDev_f
      lh$SB0 <- SB0
      lh$F_ft <- F_ft
      lh$F_fy <- F_fy
      lh$LF_tf <- LF_tf
      lh$LF0_tf <- LF0_tf
      lh$LF <- LF_tf
      lh$ML_ft <- ML_ft
      lh$R_t <- R_t
      lh$N_t <- N_t
      lh$SB_t <- SB_t
      lh$D_t <- D_t
      lh$SPR_t <- SPR_t
      lh$Cn_ft <- Cn_ft
      lh$Cw_ft <- Cw_ft 
      lh$F_t <- F_t
      lh$F_y <- F_y
      lh$F40 <- F40
      lh$Fmax <- Fmax
      lh$I_ft <- I_ft
      lh$fleet_proportions <- fleet_proportions


      return(lh)

    }) ## end with function

  }