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
#' @param fleet_percentage vector specifying the relative size of each fleet in terms of fishing pressure. must have length = nfleets and sum to 1.
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
           fleet_percentage) {
    ## SB_t = spawning biomass over time
    ## F_t = fishing mortality over time
    ## Cn_at = number of individuals that die from fishing mortality
    ## N_at = abundance by number at age over time

    with(lh, {
      ##########################
      ## Initial calcs
      ##########################

      Nyears_real <- Nyears
      Nyears <- Nyears * nseasons
      tyears <- Nyears_real * nseasons


      ##########################
      ## Random variables
      ##########################
      set.seed(seed)
    
      ## recruitment deviations
      RecDev <- rnorm(tyears, mean = -(SigmaR ^ 2) / 2, sd = SigmaR)
    
      ## autocorrelated recruitment deviations
      RecDev_AR <- rep(NA, length(RecDev))
      RecDev_AR[1] <- RecDev[1]
      for (t in 2:length(RecDev)) {
        RecDev_AR[t] <-
          RecDev_AR[t - 1] * rho + sqrt(1 - rho ^ 2) * RecDev[t]
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


      #########################
      ## Setup recruitment
      #########################
      if (Rdynamics == "Pulsed"){
        choose <- sample(1:2, 1)
        R_t <- c(
          "initial" = rep(R0, floor(Nyears / 3)),
          "pulse1" = ifelse(choose==1, rep(R0 / 2, floor(Nyears / 3)), rep(R0 * 2, floor(Nyears / 3))),
          "pulse2" = ifelse(choose==1, rep(R0 *2, Nyears - (2 * floor(Nyears / 3))), rep(R0/2, Nyears - (2 * floor(Nyears / 3))))
        ) * exp(RecDev_AR)
      }

      if (Rdynamics == "Constant") {
        R_t <- (rep(R0, tyears) / nseasons) * exp(RecDev_AR)
      }



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
            ref = 0.4
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
            ref = init_depl
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
            ref = 0.05
          )$root,
          error = function(e)
            NA
        )
      if (is.na(Fmax) | Fmax > 3)
        Fmax <- 3

      if(length(Fdynamics)==1 & nfleets>1) Fdynamics <- rep(Fdynamics, nfleets)
      
      ## effort dynamics
      E_ft <- matrix(NA, nrow=nfleets, ncol=tyears)
      if(any(Fdynamics=="Constant")){
        index <- which(Fdynamics=="Constant")
        for(i in 1:length(index)){
          E_ft[index[i],] <- 1
        }
      }
      if(any(Fdynamics=="Oneway")){
        index <- which(Fdynamics=="Oneway")
        for(i in 1:length(index)){
          E_ft[index[i],] <- c(seq(1,by=0.05,length=Nyears))
        }
      }
      if(any(Fdynamics=="Endogenous")){
        index <- which(Fdynamics=="Endogenous")
        for(i in 1:length(index)){
          E_ft[index[i],1] <- 1
        }
      }
      if (any(Fdynamics == "None")){
        index <- which(Fdynamics=="None")
        for(i in 1:length(index)){
          E_ft[index[i],] <- rep(0, tyears)
        }
      }

      ## include relative catchability by fleet -- not necessarily due to gear efficiency but due to size of fleet in practice
      qE_ft <- t(sapply(1:nfleets, function(x){
        return(E_ft[x,]*fleet_percentage[x])
      }))

      ## fishing mortality = include selectivity and Finit with effort dynamics and relative weight of fishery to scale each fishery
      F_atf <- array(NA, dim=c(length(ages), tyears, nfleets))
      for(f in 1:nfleets){
        for(t in 1:tyears){
          for(a in 1:length(ages)){
            F_atf[a,t,f] <- qE_ft[f,t] * Finit * S_fa[f,a] * exp(FishDev_f[f,t])
          }
        }
      }



      # ## fishing dynamics after burn-in
      # if (any(Fdynamics == "Ramp")){
      #   index <- which(Fdynamics=="Ramp")
      #   for(i in 1:length(index)){
      #     E_ft[index[i],(nburn+1):tyears] <-
      #     c("rampup" = seq(Finit, Fmax, length = floor(Nyears / 2)),
      #       "peak" = rep(Fmax, floor((
      #         Nyears - floor(Nyears / 2)
      #       ) / 2)),
      #       "managed" = rep(Fmax / 2, Nyears - floor(Nyears / 2) - floor((
      #         Nyears - floor(Nyears / 2)
      #       ) / 2))
      #     ) * exp(FishDev_f[index[i],(nburn+1):tyears])
      #   }
      # }

      # if (any(Fdynamics == "Decreasing")){
      #   index <- which(Fdynamics=="Decreasing")
      #   for(i in 1:length(index)){
      #     E_ft[index[i],(nburn+1):tyears] <- seq(Finit, 0, length=Nyears) * exp(FishDev_f[index[i],(nburn+1):tyears])
      #   }
      # }


      # if(is.numeric(Fdynamics) & mgt_type == "F"){ 
      #   E_ft <- t(sapply(1:nfleets, function(x){
      #     rep(Fdynamics[x], tyears) * exp(FishDev_f[x,])
      #   }))
      # }
      # if(is.numeric(Fdynamics) & mgt_type == "catch"){
      #   C_ft <- t(sapply(1:nfleets, function(x){
      #     rep(Fdynamics[x], tyears) * exp(CatchDev_f[x,])
      #   }))
      #   E_ft[,1] <- Finit
      # }




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
          N_at[a, 1] <- N_at[a - 1, 1] * exp(-M - sum(F_atf[a,1,]))
          N_at0[a, 1] <- N_at0[a - 1, 1] * exp(-M)
        }
        if (a == length(L_a)) {
          N_at[a, 1] <-
            (N_at[a - 1, 1] * exp(-M - sum(F_atf[a,1,]))) / (1 - exp(-M - sum(F_atf[a,1,])))
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
      SB0 <- sum(N_at0[, 1] * Mat_a * W_a)

      ## relative biomass (depletion) over time
      D_t <- rep(NA, tyears)
      D_t[1] <- SB_t[1] / SB0

      ##########################
      ## Projection
      ##########################

      for (t in 2:tyears) {
        ## fishing effort based on spawning biomass
        if (any(Fdynamics == "Endogenous")) {
          index <- which(Fdynamics == "Endogenous")
          for(i in 1:length(index)){
            E_ft[index[i],t] <- (E_ft[index[i],t-1] * (SB_t[t-1] / (Fequil * SB0/2)) ^ Frate)
          }
          ## include relative catchability by fleet
          qE_ft <- t(sapply(1:nfleets, function(x){
            return(E_ft[x,]*fleet_percentage[x])
          }))   

          ## fishing mortality = include selectivity and Finit with effort dynamics and relative weight of fishery to scale each fishery
          for(a in 1:length(ages)){
            F_atf[a,t,index[i]] <- qE_ft[index[i],t] * Finit * S_fa[index[i],a] * exp(FishDev_f[index[i],t])
          }
        }


        ## age-structured dynamics
        for (a in 1:length(L_a)) {
          if (a == 1) {
            N_at[a, t] <- R_t[t]
            N_at0[a, t] <- R_t[t]
          }
          if (a > 1 & a < length(L_a)) {
            N_at[a, t] <- N_at[a - 1, t - 1] * exp(-M - sum(F_atf[a-1,t-1,]))
            N_at0[a, t] <- N_at0[a - 1, t - 1] * exp(-M)
          }
          if (a == length(L_a)) {
            N_at[a, t] <- (N_at[a - 1, t - 1] * exp(-M - sum(F_atf[a-1,t-1,]))) + (N_at[a, t - 1] * exp(-M - sum(F_atf[a,t-1,])))
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

          D_t[t] <- SB_t[t] / SB0
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
      N_t <- colSums(N_at[-1, which(1:tyears %% nseasons == 0)])
      SB_t <- SB_t[which(1:tyears %% nseasons == 0)]
      D_t <- D_t[which(1:tyears %% nseasons == 0)]
      TB_t <- TB_t[which(1:tyears %% nseasons == 0)]
      SPR_t <- SPR_t[which(1:tyears %% nseasons == 0)]

      if(length(qcoef)!=nfleets) qcoef <- rep(qcoef, nfleets)
      I_ft <- t(sapply(1:nfleets, function(x) qcoef[x] * TB_t * exp(IndexDev_f[x,])))

      Cn_ft <- t(sapply(1:nfleets, function(y){
          sapply(1:Nyears_real, function(x) {
            if (nseasons == 1)
              time_index <- x
            if (nseasons > 1)
              time_index <- (1:nseasons) + ((x - 1) * nseasons)
            sum(Cn_ft[y,time_index])
          }) #* exp(CatchDev - (SigmaC^2)/2)
        }))
      Cw_ft <- t(sapply(1:nfleets, function(y){
          sapply(1:Nyears_real, function(x) {
            if (nseasons == 1)
              time_index <- x
            if (nseasons > 1)
              time_index <- (1:nseasons) + ((x - 1) * nseasons)
            sum(Cw_ft[y,time_index])
          }) #* exp(CatchDev - (SigmaC^2)/2)
      }))

      F_ft <- t(sapply(1:nfleets, function(y){
          sapply(1:Nyears_real, function(x) {
            if (nseasons == 1)
              time_index <- x
            if (nseasons > 1)
              time_index <- (1:nseasons) + ((x - 1) * nseasons)
            sum(F_ft[y,time_index])
          })
      }))
      R_t <- sapply(1:Nyears_real, function(x) {
        if (nseasons == 1)
          time_index <- x
        if (nseasons > 1)
          time_index <- (1:nseasons) + ((x - 1) * nseasons)
        sum(R_t[time_index])
      })


      ## age to length comp
      if(length(comp_sample)!=nfleets) comp_sample <- rep(comp_sample, nfleets)
      obs_per_year <- t(sapply(1:nfleets, function(x){
          rep(comp_sample[x] / nseasons, tyears)
      }))

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

      if (pool == TRUE) {
        LF_tf <- LF0_tf <- array(NA, dim=c(Nyears_real, ncol(LF[[1]]), nfleets))
        for(f in 1:nfleets){
         for (y in 1:Nyears_real) {
            if (nseasons == 1) {
              LF_tf[y,,f] <- LF[[f]][y,]
              LF0_tf[y,,f] <- LF0[[f]][y,]
            }
            if (nseasons > 1) {
              time_index <- (1:nseasons) + ((y - 1) * nseasons)
              LF_tf[y,,f] <- colSums(LF[[f]][time_index,])
              LF0_tf[y,,f] <- colSums(LF0[[f]][time_index,])
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
      if (pool == FALSE) {
        LF_t <- LF
        LF0_t <- LF0
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
      if (pool == TRUE)
        ML_ft <- t(sapply(1:nfleets, function(x) ML_ft[x,which(1:tyears %% nseasons == 0)]))

      ## outputs
      lh$I_ft <- I_ft
      lh$Cn_ft <- Cn_ft
      lh$Cw_ft <- Cw_ft
      # if(is.numeric(Fdynamics) & mgt_type=="catch") lh$C_t <- C_t
      lh$LF <- LF
      lh$LF0 <- LF0
      lh$R_t <- R_t
      lh$N_t <- N_t
      lh$SB_t <- SB_t
      lh$D_t <- D_t
      lh$F_t <- F_t
      lh$F_ft <- F_ft
      lh$ML_ft <- ML_ft
      lh$plb <- plb
      lh$plba <- plba
      lh$page <- page
      lh$N_at <- N_at
      lh$SPR <- SPR
      lh$SPR_t <- SPR_t
      lh$TB_t <- TB_t
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


      return(lh)

    }) ## end with function

  }
