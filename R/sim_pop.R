#' Operating model
#'
#' \code{sim_pop} Age-converted-to-length-based operating model specifying true population dynamics
#'
#' @author M.B. Rudd
#' @param lh list of life history information, from create_lh_list
#' @param Nyears number of years to simulate
#' @param pool if nseasons (in life history list) is greater than one, pool the generated data into annual time steps, or leave at the season-level? FALSE will generate shorter time step life history info, mean length
#' @param Fdynamics Specify name of pattern of fishing mortality dynamics, Constant, Endogenous, Oneway, or None. Input matrix with dimensions fleets = rows, years = columns
#' @param Rdynamics Specify name of pattern of recruitment dynamics, Constant, Pulsed, Pulsed_up, or BH
#' @param Nyears_comp number of years of length composition data
#' @param comp_sample sample size of length composition data annually
#' @param init_depl default=NULL, can specify a starting depletion, initial F (using init_F), or 2 values that indicate range from which to choose them
#' @param init_F default=NULL, can specify a starting F, an initial depletion (using init_depl), or 2 vaues that indicate range from which to choose them 
#' @param seed set seed for generating stochastic time series
#' @param sample_type a character vector specifying if the length comps are sampled from the 'catch' (default) or from the population
#' @param mgt_type removals based on F (default) or catch
#' @param fleet_proportions vector specifying the relative size of each fleet in terms of fishing pressure. must have length = nfleets and sum to 1.
#' @param nareas number of areas, default = 1, if greater than 1, must be equal to the number of fleets
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
           init_depl = NULL,
           init_F = NULL,
           seed,
           sample_type = 'catch',
           mgt_type = 'F',
           fleet_proportions,
           nareas) {

    with(lh, {
      ##########################
      ## Initial calcs
      ##########################
      if(Nyears == 1) stop("Must simulate at least two years")
      Nyears_real <- Nyears

      Nyears <- Nyears * nseasons
      Nyears_comp <- Nyears_comp * nseasons

      if((nareas > 1) & (nareas != nfleets)) stop("Number of fleets must be equal to the number of areas if nareas > 1")


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
      # if(nseasons > 1){
      #   # RecDev2 <- rep(0, Nyears)
      #   # RecDev2[seq(1,length(RecDev2),by=nseasons)] <- RecDev_AR
      #   RecDev2 <- unlist(lapply(1:Nyears_real, function(x) rep(RecDev_AR[x], nseasons)))
      #   RecDev_AR <- RecDev2
      # }

      ## fishing mortality deviations
      if(length(SigmaF)==1 & nfleets>1) SigmaF <- rep(SigmaF, nfleets)
      FishDev_f <- t(sapply(1:nfleets, function(x){
        c(0, rnorm(Nyears-1, mean = -(SigmaF[x] ^ 2) / 2, sd = SigmaF[x]))
      }))
      
    
      ## abundance index observation error
      if(length(SigmaI)==1 & nfleets>1) SigmaI <- rep(SigmaI, nfleets)
      IndexDev_f <- t(sapply(1:nfleets, function(x){
        rnorm(Nyears, mean = -(SigmaI[x] ^ 2) / 2, sd = SigmaI[x])
      }))

      ## catch observation error
      if(length(SigmaC)==1 & nfleets>1) SigmaC <- rep(SigmaC, nfleets)
      CatchDev_f <- t(sapply(1:nfleets, function(x){
        rnorm(Nyears, mean = -(SigmaC[x] ^ 2) / 2, sd = SigmaC[x])
      }))


      #########################
      ## Setup recruitment
      #########################
      if(length(Rdynamics) > 1){
        R_t <- (rep(R0, Nyears_real)) * exp(Rdynamics)
      }
      if (Rdynamics == "Pulsed"){
        choose <- sample(1:2, 1)
        R_t <- c(
          "initial" = rep(R0, floor(Nyears / 3)),
          "pulse1" = rep(ifelse(choose==1, (R0 / 2), (R0 * 2)), floor(Nyears / 3)),
          "pulse2" = rep(ifelse(choose==1, (R0 *2), (R0 / 2)), Nyears - (2 * floor(Nyears / 3)))
        ) * exp(RecDev_AR)
        # R_t[which(RecDev_AR==0)] <- 0
      }

      if (Rdynamics == "Constant") {
        R_t <- (rep(R0, Nyears_real)) * exp(RecDev_AR)
        # R_t[which(RecDev_AR==0)] <- 0
      }

      if(Rdynamics=="BevertonHolt"){
        R_t <- rep(NA, Nyears_real)
        R_t[1] <- R0 * exp(RecDev_AR[1])
      }

    ## with multi-seasons, currently spreading out recruitment across seasons in one year
    ## to properly calculate initial F if recruitment is only in one season, need to spread out calculation across the year instead of per-recruit in calc_ref
        # R_t <- R_t/nseasons

      ## multi-seasons - one recruitment event per year
      R_t_new <- unlist(lapply(1:Nyears_real, function(x) c(R_t[x], rep(0,nseasons-1))))
      R_t <- R_t_new

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
            S_fa = S_fa, 
            ref = 0.4,
            fleet_prop=fleet_proportions,
            type = "spr")$root,
          error = function(e)
            NA
        )
     if(is.matrix(Fdynamics)==FALSE){
      if(all(is.null(init_F))==FALSE){
        Finit <- init_F
      } else {
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
              S_fa = S_fa, 
              ref = init_depl,
              fleet_prop=fleet_proportions,
              type = "depletion")$root,
            error = function(e)
              NA
          )        
      }
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
            S_fa = S_fa, 
            ref = 0.05,
            fleet_prop=fleet_proportions,
            type = "depletion")$root,
          error = function(e)
            NA
        )
      if (is.na(Fmax) | Fmax > 3)
        Fmax <- 3

      if(length(Fdynamics)==1 & nfleets>1) Fdynamics <- rep(Fdynamics, nfleets)
      
      ## fishing mortality dynamics
      F_ft <- matrix(NA, nrow=nfleets, ncol=Nyears)
      if(any(Fdynamics=="Constant")){
        index <- which(Fdynamics=="Constant")
        for(i in 1:length(index)){
          F_ft[index[i],] <- Finit * fleet_proportions[index[i]] * exp(FishDev_f[index[i],])
        }
      }
      if(any(Fdynamics=="Oneway")){
        index <- which(Fdynamics=="Oneway")
        for(i in 1:length(index)){
          F_ft[index[i],] <- c(seq(1,by=0.05,length=Nyears)) * Finit * fleet_proportions[index[i]] * exp(FishDev_f[index[i],])
        }
      }
      if(any(Fdynamics=="Endogenous")){
        index <- which(Fdynamics=="Endogenous")
        for(i in 1:length(index)){
          F_ft[index[i],1] <- Finit * fleet_proportions[index[i]]
        }
      }
      if (any(Fdynamics == "None")){
        index <- which(Fdynamics=="None")
        for(i in 1:length(index)){
          F_ft[index[i],] <- rep(0, Nyears) * exp(FishDev_f[index[i],])
        }
      }
    }
    if(is.matrix(Fdynamics)){
      F_ft <- Fdynamics
      for(i in 1:nrow(F_ft)){
        F_ft[i,] <- F_ft[i,] * exp(FishDev_f[i,])
      }
    }

      ## fishing mortality = include selectivity and Finit with effort dynamics and relative weight of fishery to scale each fishery
      F_atf <- array(NA, dim=c(length(ages), Nyears, nfleets))
      for(f in 1:nfleets){
        for(t in 1:Nyears){
          for(a in 1:length(ages)){
            F_atf[a,t,f] <- F_ft[f,t] * S_fa[f,a]
          }
        }
      }

      F_at <- matrix(NA, nrow=length(ages), ncol=Nyears)
      for(t in 1:Nyears){
        for(a in 1:length(ages)){
          F_at[a,t] <- sum(F_atf[a,t,])
        }
      }


      ##########################
      ## Data objects
      ##########################

      ## year 1 abundance at age
      N_ats <- N_ats0 <- array(NA, dim=c(length(L_a), Nyears, nareas)) #nrow = length(L_a), ncol = Nyears)
      N_at <- N_at0 <- matrix(0, nrow = length(L_a), ncol = Nyears)

      for (a in 1:length(L_a)) {
        if (a == 1) {
          for(i in 1:nareas){
            N_ats[a,1,i] <- R_t[1]/nareas
            N_ats0[a,1,i] <- R_t[1]/nareas
          }
        }
        if (a > 1 & a < length(L_a)) {
          if(nareas > 1){
            for(i in 1:nareas){
              N_ats[a,1,i] <- N_ats[a - 1, 1, i] * exp(-M - F_atf[a-1,1,i])
              N_ats0[a,1,i] <- N_ats0[a - 1, 1, i] * exp(-M)
            }            
          }
          if(nareas == 1){
            N_ats[a,1,1] <- N_ats[a-1, 1, 1] * exp(-M - F_at[a-1,1])
            N_ats0[a,1,1] <- N_ats[a-1,1,1] * exp(-M)
          }

        }
        if (a == length(L_a)) {
          if(nareas > 1){
            for(i in 1:nareas){
              N_ats[a, 1, i] <- (N_ats[a - 1, 1, i] * exp(-M - F_atf[a-1,1,i])) / (1 - exp(-M - F_atf[a-1,1,i]))
              N_ats0[a, 1, i] <- (N_ats0[a - 1, 1, i] * exp(-M)) / (1 - exp(-M))
            }
          }
          if(nareas == 1){
            N_ats[a, 1,1] <- (N_ats[a - 1, 1,1] * exp(-M - F_at[a-1,1])) / (1 - exp(-M - F_at[a-1,1]))
            N_ats0[a, 1,1] <- (N_ats0[a - 1, 1,1] * exp(-M)) / (1 - exp(-M))        
          }
        }
      }

      for(a in 1:length(L_a)){
        for(i in 1:nareas){
          N_at[a,1] <- N_at[a,1] + N_ats[a,1,i]
          N_at0[a,1] <- N_at0[a,1] + N_ats0[a,1,i]
        }
      }

      ## year 1 biomass quantities
      TB_t <- SB_t <- rep(NA, Nyears)
      TB_ts <- matrix(NA, nrow=Nyears, ncol=nareas)
      TB_t[1] <- sum(N_at[, 1] * W_a)
      SB_t[1] <- sum(N_at[, 1] * W_a * Mat_a) * 0.5
      for(i in 1:nareas){
        TB_ts[,i] <- sum(N_ats[,1,i] * W_a)
      }

      if(is.numeric(Fdynamics) & length(Fdynamics)== 1 & mgt_type=="catch"){
        F_ft[1,] <- max(0.01, getFt(ct=Fdynamics, m=M, sa=S_fa[1,], wa=W_a, na=N_at[,1]))
        # F_ft[1,] <- min(c(Fmax, F_ft[1,1]), na.rm=TRUE)
      }

      ## year 1 catch
      Cn_atf <- Cw_atf <- array(NA, c(length(L_a), Nyears, nfleets))
      for(f in 1:nfleets){
        if(nareas == 1){
          Cn_atf[,1,f] <- N_at[,1] * (1 - exp(-M - F_atf[,1,f])) * (F_atf[,1,f] / (M + F_atf[a,1,f]))
          Cw_atf[,1,f] <- Cn_atf[,1,f] * W_a
        }
        if(nareas > 1){
          Cn_atf[,1,f] <- N_ats[,1,f] * (1 - exp(-M - F_atf[,1,f])) * (F_atf[,1,f] / (M + F_atf[a,1,f]))
          Cw_atf[,1,f] <- Cn_atf[,1,f] * W_a          
        }
      }

      ## unfished spawning biomass
      SB0 <- sum(calc_equil_abund(ages=ages, M=M, F=rep(0,nrow(S_fa)), R0=R0, S_fa=S_fa) * W_a * Mat_a * 0.5)       


      ##########################
      ## Projection
      ##########################

      for (t in 2:Nyears) {
        ## fishing effort based on spawning biomass
        if(Rdynamics=="BevertonHolt"){
          R_t[t] <- ((4 * h * R0 * SB_t[t-1]) / (SB0 * (1-h) + SB_t[t-1] * (5*h-1))) * exp(RecDev_AR[t])
        }
        if (any(Fdynamics == "Endogenous")) {
          index <- which(Fdynamics == "Endogenous")
          for(i in 1:length(index)){
            F_ft[index[i],t] <- (F_ft[index[i],t-1] * (SB_t[t-1] / (Fequil * SB0)) ^ Frate) * exp(FishDev_f[index[i],t])
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
            for(i in 1:nareas){
              N_ats[a, t, i] <- R_t[t]/nareas
              N_ats0[a, t, i] <- R_t[t]/nareas
            }
          }
          if (a > 1 & a < length(L_a)) {
            if(nareas > 1){
              for(i in 1:nareas){
                N_ats[a, t, i] <- N_ats[a - 1, t - 1, i] * exp(-M - F_atf[a-1,t-1,i])
                N_ats0[a, t, i] <- N_ats0[a - 1, t - 1, i] * exp(-M)
              }              
            }
            if(nareas == 1){
                N_ats[a, t, 1] <- N_ats[a - 1, t - 1, 1] * exp(-M - F_at[a-1,t-1])
                N_ats0[a, t, 1] <- N_ats0[a - 1, t - 1, 1] * exp(-M)              
            }
          }
          if (a == length(L_a)) {
            if(nareas > 1){
              for(i in 1:nareas){
                N_ats[a, t, i] <- (N_ats[a - 1, t - 1, i] * exp(-M - F_atf[a-1,t-1,i])) + (N_ats[a, t - 1, i] * exp(-M - F_atf[a,t-1,i]))
                N_ats0[a, t, i] <- (N_ats0[a - 1, t - 1, i] * exp(-M)) + (N_ats0[a, t - 1,i] * exp(-M))
              }
            }
            if(nareas == 1){
              N_ats[a, t,1] <- (N_ats[a - 1, t - 1,1] * exp(-M - F_at[a-1,t-1])) + (N_ats[a, t - 1,1] * exp(-M - F_at[a,t-1]))
              N_ats0[a, t,1] <- (N_ats0[a - 1, t - 1,1] * exp(-M)) + (N_ats0[a, t - 1,1] * exp(-M))
            }
          }
        }

        for(a in 1:length(L_a)){
          for(i in 1:nareas){
            N_at[a,t] <- N_at[a,t] + N_ats[a,t,i]
            N_at0[a,t] <- N_at0[a,t] + N_ats0[a,t,i]
          }
        }

          ## spawning biomass
          SB_t[t] <- sum((N_at[, t] * W_a * Mat_a * 0.5))
          TB_t[t] <- sum(N_at[, t] * W_a)
          for(i in 1:nareas){
            TB_ts[t,i] <- sum(N_ats[,t,i] * W_a)
          }

          if(is.numeric(Fdynamics) & mgt_type=="catch"){
            F_ft[1,y] <- max(0.01,getFt(ct=Fdynamics, m=M, sa=S_fa[1,], wa=W_a, na=N_at[,y]))
            # F_t[y] <- min(c(Fmax, F_t[y]), na.rm=TRUE)
          }

          ## catch
          for(f in 1:nfleets){
            if(nareas == 1){
              Cn_atf[,t,f] <- N_at[,t] * (1 - exp(-M - F_atf[,t,f])) * (F_atf[,t,f] / (M + F_atf[,t,f]))
              Cw_atf[,t,f] <- Cn_atf[,t,f] * W_a
            }
            if(nareas > 1){
              Cn_atf[,t,f] <- N_ats[,t,f] * (1 - exp(-M - F_atf[,t,f])) * (F_atf[,t,f] / (M + F_atf[,t,f]))
              Cw_atf[,t,f] <- Cn_atf[,t,f] * W_a              
            }
          }

      }

      F_ft <- t(sapply(1:nfleets, function(x){
        sub <- F_atf[,,x]
        findMax <- sapply(1:Nyears, function(y){
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
            S_fa = S_fa,
            F = F_t[x],
            fleet_prop=fleet_proportions
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
          N_t <- colSums(N_at[-1, which(1:Nyears %% nseasons == 0)])
          SB_t <- SB_t[which(1:Nyears %% nseasons == 0)]
          TB_t <- TB_t[which(1:Nyears %% nseasons == 0)]
          SPR_t <- SPR_t[which(1:Nyears %% nseasons == 0)]

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
      if(nareas > 1){
        I_ft <- t(sapply(1:nfleets, function(x) qcoef[x] * TB_ts[,x] * exp(IndexDev_f[x,])))
        Effort_ft <- t(sapply(1:nfleets, function(x) Cw_ft[x,]/(qcoef[x]*TB_ts[,x])))
      }
      if(nareas == 1){
        I_ft <- t(sapply(1:nfleets, function(x) qcoef[x] * TB_t * exp(IndexDev_f[x,])))
        Effort_ft <- t(sapply(1:nfleets, function(x) Cw_ft[x,]/(qcoef[x]*TB_t)))
      }

      ## age to length comp
      if(length(Nyears_comp)!=nfleets) Nyears_comp <- rep(Nyears_comp, nfleets)

      ## years with observed length comps
      oyears_mat <- matrix(0, nrow=nfleets, ncol=Nyears)
      for(f in 1:nfleets){
        oyears <- (Nyears-Nyears_comp[f] + 1):Nyears
        oyears_mat[f,oyears] <- 1
      }

      obs_per_year <- matrix(0, nrow=nfleets, ncol=Nyears)
      for(f in 1:nfleets){
        for(t in 1:Nyears){
          if(oyears_mat[f,t]!=0) obs_per_year[f,t] <- (comp_sample[f]/nseasons)
        }
      }

      LFinfo <-lapply(1:nfleets, function(x){
        if(nareas == 1){
          lf <- AgeToLengthComp(
            lh = lh,
            S_a = lh$S_fa[x,],
            tyears = Nyears,
            N_at = N_at,
            comp_sample = obs_per_year[x,],
            sample_type = sample_type)          
        }
        if(nareas > 1){
          lf <- AgeToLengthComp(
            lh = lh,
            S_a = lh$S_fa[x,],
            tyears = Nyears,
            N_at = N_ats[,,x],
            comp_sample = obs_per_year[x,],
            sample_type = sample_type)
        } 
        return(lf)
      })
      LF0info <- lapply(1:nfleets, function(x){
        if(nareas == 1){
         lf <- AgeToLengthComp(
          lh = lh,
          S_a = lh$S_fa[x,],
          tyears = Nyears,
          N_at = N_at0,
          comp_sample = obs_per_year[x,],
          sample_type = sample_type)         
        }
        if(nareas > 1){
         lf <- AgeToLengthComp(
          lh = lh,
          S_a = lh$S_fa[x,],
          tyears = Nyears,
          N_at = N_ats0[,,x],
          comp_sample = obs_per_year[x,],
          sample_type = sample_type)   
        }
        return(lf)
      })
      plba <- lapply(1:nfleets, function(x) LFinfo[[x]]$plba)
      plb <- lapply(1:nfleets, function(x) LFinfo[[x]]$plb)
      page <- lapply(1:nfleets, function(x) LFinfo[[x]]$page)
      LF <- lapply(1:nfleets, function(x) LFinfo[[x]]$LF)
      LF0 <- lapply(1:nfleets, function(x) LF0info[[x]]$LF)

      LF_tf <- LF0_tf <- NULL
      if (pool == TRUE) {
        for(f in 1:nfleets){
          LF_tf[[f]] <- LF0_tf[[f]] <- matrix(NA, nrow=Nyears, ncol=ncol(LF[[1]]))
         for (y in 1:Nyears) {
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
        colnames(LF_tf[[f]]) <- mids
        colnames(LF0_tf[[f]]) <- mids
        if(pool==TRUE){
          rownames(LF_tf[[f]]) <- 1:Nyears_real
          rownames(LF0_tf[[f]]) <- 1:Nyears_real
        }
      }


      ########################################################
      ## Expected mean length in catch
      ########################################################
      ML_ft <- matrix(NA, nrow=nfleets, ncol=Nyears)
      for(f in 1:nfleets){
        for (t in 1:Nyears) {
          vul_pop <- sum(N_at[, t] * S_fa[f,])
          vul_lengths <- sum(vul_pop * plb[[f]][t,] * mids)
          ML_ft[f,t] <- vul_lengths / vul_pop
        }
      }
      if (pool == TRUE) ML_ft <- t(sapply(1:nfleets, function(x) ML_ft[x,which(1:Nyears %% nseasons == 0)]))

      lh$plb <- lapply(1:length(plb), function(x) plb[[x]])
      lh$plba <- plba
      lh$page <- lapply(1:length(page), function(x) page[[x]])
      lh$N_at <- N_at
      lh$N_ats <- N_ats
      lh$nlbins <- length(mids)
      if (pool == TRUE) {
        lh$Nyears <- Nyears_real
        lh$years <- 1:Nyears_real
      }
      if (pool == FALSE) {
        lh$Nyears <- Nyears
        lh$years <- 1:Nyears
      }
      obs_per_year <- matrix(obs_per_year, nrow=nfleets, ncol=Nyears)
      colnames(obs_per_year) <- 1:Nyears
      lh$obs_per_year <- obs_per_year
      lh$RecDev <- RecDev_AR
      lh$FishDev_f <- matrix(FishDev_f, nrow=nfleets, ncol=Nyears_real)
      lh$SB0 <- SB0
      lh$F_ft <- matrix(F_ft, nrow=nfleets, ncol=Nyears)
      lh$F_fy <- matrix(F_fy, nrow=nfleets, ncol=Nyears_real)

      LF_tfout <- lapply(1:nfleets, function(x){
        if(pool==TRUE){
          sub <- matrix(LF_tf[[x]], nrow=Nyears_real, ncol=length(mids))
          colnames(sub) <- mids
          rownames(sub) <- 1:Nyears_real
        }
        if(pool==FALSE){
          sub <- matrix(LF_tf[[x]], nrow=Nyears, ncol=length(mids))
          colnames(sub) <- mids
          rownames(sub) <- 1:Nyears
        }
        return(sub)
      })
      LF0_tfout <- lapply(1:nfleets, function(x){
        if(pool==TRUE){
          sub <- matrix(LF0_tf[[x]], nrow=Nyears_real, ncol=length(mids))
          colnames(sub) <- mids
          rownames(sub) <- 1:Nyears_real
        }
        if(pool==FALSE){
          sub <- matrix(LF0_tf[[x]], nrow=Nyears, ncol=length(mids))
          colnames(sub) <- mids
          rownames(sub) <- 1:Nyears
        }
        return(sub)
      })
      lh$LF_tlf <- LF_tfout
      lh$LF0_tlf <- LF0_tfout
      lh$LF <- LF_tfout
      lh$ML_ft <- matrix(ML_ft, nrow=nfleets, ncol=Nyears)
      lh$R_t <- R_t
      lh$N_t <- N_t
      lh$SB_t <- SB_t
      lh$D_t <- D_t
      lh$SPR_t <- SPR_t
      lh$SPR <- SPR
      lh$Cn_ft <- matrix(Cn_ft, nrow=nfleets, ncol=Nyears)
      lh$Cw_ft <- matrix(Cw_ft, nrow=nfleets,  ncol=Nyears) 
      lh$F_t <- F_t
      lh$F_y <- F_y
      lh$F40 <- F40
      # lh$Fmax <- Fmax
      lh$I_ft <- matrix(I_ft, nrow=nfleets, ncol=Nyears)
      lh$E_ft <- matrix(Effort_ft, nrow=nfleets, ncol=Nyears)
      lh$fleet_proportions <- fleet_proportions


      return(lh)

    }) ## end with function

  }