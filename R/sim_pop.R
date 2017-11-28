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
#' @param init_depl initial depletion; if FALSE, will use F1 from lh list
#' @param nburn number of years of burn-in for operating model
#' @param seed set seed for generating stochastic time series
#' @param mismatch if TRUE, catch and index overlap with length comp only 1 year
#' @param sample_type a character vector specifying if the length comps are sampled from the 'catch' (default) or from the population
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
           init_depl,
           nburn,
           seed,
           mismatch,
           sample_type = 'catch') {
    ## SB_t = spawning biomass over time
    ## F_t = fishing mortality over time
    ## Cn_at = number of individuals that die from fishing mortality
    ## N_at = abundance by number at age over time

    with(lh, {
      ##########################
      ## Initial calcs
      ##########################

      tyears_only <- nburn + Nyears
      Nyears_real <- Nyears
      nburn_real <- nburn
      nburn <- nburn * nseasons
      Nyears <- Nyears * nseasons
      tyears <- tyears_only * nseasons


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
      FishDev <-
        rnorm(tyears, mean = -(SigmaF ^ 2) / 2, sd = SigmaF)

      ## abundance index observation error
      IndexDev <-
        rnorm(tyears, mean = -(SigmaI ^ 2) / 2, sd = SigmaI)

      ## catch observation error
      CatchDev <-
        rnorm(tyears, mean = -(SigmaC ^ 2) / 2, sd = SigmaC)

      ##########################
      ## Data objects
      ##########################
      TB_t <-
        VB_t <- SB_t <- F_t <- R_t <- D_t <- Z_t <- rep(NA, tyears)
      Cn_at <-
        N_at <-
        N_at0 <- matrix(NA, nrow = length(L_a), ncol = tyears)

      #####################################
      ## Fishing and recruitment dynamics
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
            S_a = S_a,
            ref = 0.4
          )$root,
          error = function(e)
            NA
        )
      if (init_depl == FALSE)
        Finit <- F1
      if (init_depl != FALSE)
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
            S_a = S_a,
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
            S_a = S_a,
            ref = 0.2
          )$root,
          error = function(e)
            NA
        )
      if (is.na(Fmax) | Fmax > 3)
        Fmax <- 3

      if (Fdynamics == "Ramp")
        Framp_t <-
        c(
          rep(Finit, nburn),
          "rampup" = seq(Finit, Fmax, length = floor(Nyears / 2)),
          "peak" = rep(Fmax, floor((
            Nyears - floor(Nyears / 2)
          ) / 2)),
          "managed" = rep(Fmax / 2, Nyears - floor(Nyears / 2) - floor((
            Nyears - floor(Nyears / 2)
          ) / 2))
        )
      if (Fdynamics == "Constant")
        Fconstant_t <- rep(Finit, tyears)
      if (Fdynamics == "Increasing")
        Finc_t <-
        c(rep(Finit, nburn), seq(Finit, Fmax, length = Nyears))
      if (Fdynamics == "Decreasing")
        Fdec_t <-
        c(rep(Finit, nburn), seq(Finit, 0, length = Nyears))
      if (Fdynamics == "None")
        F_t <- rep(0, tyears)
      if (Fdynamics == "4010")
        F_t <- rep(NA, tyears)
      if(is.numeric(Fdynamics)) 
        F_t <- rep(Fdynamics, tyears) * exp(FishDev)

      if (Rdynamics == "Pulsed")
        Rpulse_t <- c(
          rep(R0, nburn),
          "initial" = rep(R0, floor(Nyears / 3)),
          "pulse_down" = rep(R0 / 3, floor(Nyears / 3)),
          "pulse_up" = rep(R0, Nyears - (2 * floor(Nyears / 3)))
        )
      if (Rdynamics == "Pulsed_up")
        Rpulse_t <-
        c(
          rep(R0, nburn),
          "initial" = rep(R0, floor(Nyears / 3)),
          "pulse_up" = rep(R0 * 3, floor(Nyears / 3)),
          "pulse_down" = rep(R0, Nyears - (2 * floor(Nyears / 3)))
        )
      if (Rdynamics == "Constant" |
          Rdynamics == "AR")
        Rconstant_t <- rep(R0, tyears)

      if (Fdynamics == "Ramp") {
        F_t <- Framp_t * exp(FishDev)
      }
      if (Fdynamics == "Constant") {
        F_t <- Fconstant_t * exp(FishDev)
      }
      if (Fdynamics == "Increasing") {
        F_t <- Finc_t * exp(FishDev)
      }
      if (Fdynamics == "Decreasing") {
        F_t <- Fdec_t * exp(FishDev)
      }
      if (Fdynamics == "Endogenous") {
        F_t[1] <- Finit
      }
      if (Fdynamics == "4010") {
        F_t[1] <- Finit
      }
      if (Rdynamics == "Constant") {
        R_t <- Rconstant_t / nseasons * exp(RecDev)
      }
      if (Rdynamics == "AR") {
        R_t <- Rconstant_t / nseasons * exp(RecDev_AR)
      }
      if (Rdynamics == "Pulsed") {
        R_t <- Rpulse_t / nseasons * exp(RecDev)
      }
      if (Rdynamics == "Pulsed_up") {
        R_t <- Rpulse_t / nseasons * exp(RecDev)
      }
      if (Rdynamics == "BH") {
        R_t[1] <- R0 / nseasons * exp(RecDev[1])
      }

      ## year 1
      for (a in 1:length(L_a)) {
        if (a == 1) {
          N_at[a, 1] <- R_t[1]
          N_at0[a, 1] <- R_t[1]
        }
        if (a > 1 & a < length(L_a)) {
          N_at[a, 1] <- N_at[a - 1, 1] * exp(-M - F_t[1] * S_a[a - 1])
          N_at0[a, 1] <- N_at0[a - 1, 1] * exp(-M)
        }
        if (a == length(L_a)) {
          N_at[a, 1] <-
            (N_at[a - 1, 1] * exp(-M - F_t[1] * S_a[a])) / (1 - exp(-M - F_t[1] * S_a[a]))
          N_at0[a, 1] <- (N_at0[a - 1, 1] * exp(-M)) / (1 - exp(-M))
        }

      }
      VB_t[1] <- sum(N_at[, 1] * W_a * S_a)
      TB_t[1] <- sum(N_at[, 1] * W_a)
      SB_t[1] <- sum(N_at[, 1] * W_a * Mat_a)
      Cn_at[, 1] <-
        N_at[, 1] * (1 - exp(-M - F_t[1] * S_a)) * (F_t[1] * S_a) / (M + F_t[1] * S_a)

      ##########################
      ## Projection
      ##########################
      SB0 <- sum(N_at0[, 1] * Mat_a * W_a)
      D_t[1] <- SB_t[1] / SB0

      # D_t <- seq(0,1,length=tyears)
      # F_t <- rep(NA, tyears)
      # F_t[1] <- 0
      # for(t in 2:tyears){
      #             if(D_t[t-1] < 0.10) F_t[t] <- 0
      #             if(D_t[t-1] >= 0.40) F_t[t] <- F40 * exp(FishDev[t])
      #             if(D_t[t-1] >= 0.10 & D_t[t-1] < 0.40) F_t[t] <- ((F40/0.3)*D_t[t-1] - ((0.10*F40)/0.30)) * exp(FishDev[t])
      # }

      for (y in 2:tyears) {
        ## fishing effort and recruitment, not dependent on age structure
        if (Fdynamics == "Endogenous") {
          if (y <= nburn)
            F_t[y] <- Finit
          if (y > nburn)
            F_t[y] <-
              F_t[y - 1] * (SB_t[y - 1] / (Fequil * SB0)) ^ Frate * exp(FishDev[y])
        }
        if (Fdynamics == "4010") {
          if (y <= nburn)
            F_t[y] <- Finit
          if (y > nburn) {
            if (D_t[y - 1] < 0.10)
              F_t[y] <- 0
            # if(D_t[y-1] >= 0.40) F_t[y] <- F40 * exp(FishDev[y])
            # if(D_t[y-1] >= 0.10 & D_t[y-1] < 0.40) F_t[y] <- ((F40/0.3)*D_t[y-1] - ((0.10*F40)/0.30)) * exp(FishDev[y])
            if (D_t[y - 1] >= 0.10)
              F_t[y] <-
                ((F40 / 0.3) * D_t[y - 1] - ((0.10 * F40) / 0.30)) * exp(FishDev[y])
          }
        }
        if (Rdynamics == "BH") {
          if (h == 1)
            h_use <- 0.7
          if (h != 1)
            h_use <- h
          R_t[y] <-
            (4 * h_use * R0 * SB_t[y - 1] / (SB0 * (1 - h_use) + SB_t[y - 1] * (5 *
                                                                                  h_use - 1))) / nseasons * exp(RecDev[y])
        }

        ## age-structured dynamics
        for (a in 1:length(L_a)) {
          if (a == 1) {
            N_at[a, y] <- R_t[y]
            N_at0[a, y] <- R_t[y]
          }
          if (a > 1 & a < length(L_a)) {
            N_at[a, y] <- N_at[a - 1, y - 1] * exp(-M - F_t[y - 1] * S_a[a - 1])
            N_at0[a, y] <- N_at0[a - 1, y - 1] * exp(-M)
          }
          if (a == length(L_a)) {
            N_at[a, y] <-
              (N_at[a - 1, y - 1] * exp(-M - F_t[y - 1] * S_a[a - 1])) + (N_at[a, y -
                                                                                 1] * exp(-M - F_t[y - 1] * S_a[a]))
            N_at0[a, y] <-
              (N_at0[a - 1, y - 1] * exp(-M)) + (N_at0[a, y - 1] * exp(-M))
          }

          ## spawning biomass
          SB_t[y] <- sum((N_at[, y] * W_a * Mat_a))
          VB_t[y] <- sum(N_at[, y] * W_a * S_a)
          TB_t[y] <- sum(N_at[, y] * W_a)

          ## catch
          Cn_at[, y] <-
            N_at[, y] * (1 - exp(-M - F_t[y] * S_a)) * (F_t[y] * S_a) / (M + F_t[y] * S_a)

          Z_t[y] <- mean(M + F_t[y] * S_a, na.rm = T)

          D_t <- SB_t / SB0


        }
      }

      ## static SPR
      SPR_t <-
        sapply(1:length(F_t), function(x)
          calc_ref(
            ages = ages,
            Mat_a = Mat_a,
            W_a = W_a,
            M = M,
            S_a = S_a,
            F = F_t[x]
          ))
      SPR <- SPR_t[length(SPR_t)]

      P <- 0.01
      x <-
        seq(from = 0,
            to = 1,
            length.out = length(L_a)) # relative age vector
      EL <-
        (1 - P ^ (x / (M / vbk))) * linf # length at relative age
      rLens <- EL / linf # relative length
      SPR_alt <-
        sum(Mat_a * rowSums(N_at) * rLens ^ 3) / sum(Mat_a * rowSums(N_at0) * rLens ^
                                                       3)

      Cn_t <- colSums(Cn_at)
      Cw_t <- colSums(Cn_at * W_a)
      N_t <- colSums(N_at[-1, which(1:tyears %% nseasons == 0)])
      SB_t <- SB_t[which(1:tyears %% nseasons == 0)]
      D_t <- D_t[which(1:tyears %% nseasons == 0)]
      TB_t <- TB_t[which(1:tyears %% nseasons == 0)]
      VB_t <- VB_t[which(1:tyears %% nseasons == 0)]
      SPR_t <- SPR_t[which(1:tyears %% nseasons == 0)]

      I_t <- qcoef * TB_t #* exp(IndexDev - (SigmaI^2)/2)
      C_t <- sapply(1:tyears_only, function(x) {
        if (nseasons == 1)
          time_index <- x
        if (nseasons > 1)
          time_index <- (1:nseasons) + ((x - 1) * nseasons)
        sum(Cn_t[time_index])
      }) #* exp(CatchDev - (SigmaC^2)/2)
      Cw_t <- sapply(1:tyears_only, function(x) {
        if (nseasons == 1)
          time_index <- x
        if (nseasons > 1)
          time_index <- (1:nseasons) + ((x - 1) * nseasons)
        sum(Cw_t[time_index])
      }) #* exp(CatchDev - (SigmaC^2)/2)
      F_t <- sapply(1:tyears_only, function(x) {
        if (nseasons == 1)
          time_index <- x
        if (nseasons > 1)
          time_index <- (1:nseasons) + ((x - 1) * nseasons)
        sum(F_t[time_index])
      })
      R_t <- sapply(1:tyears_only, function(x) {
        if (nseasons == 1)
          time_index <- x
        if (nseasons > 1)
          time_index <- (1:nseasons) + ((x - 1) * nseasons)
        sum(R_t[time_index])
      })


      ## age to length comp
      obs_per_year <- rep(comp_sample / nseasons, tyears)
      LFinfo <-
        AgeToLengthComp(
          lh = lh,
          tyears = tyears,
          N_at = N_at,
          comp_sample = obs_per_year,
          sample_type = sample_type
        )
      LF0info <-
        AgeToLengthComp(
          lh = lh,
          tyears = tyears,
          N_at = N_at0,
          comp_sample = obs_per_year,
          sample_type = sample_type
        )

      plba <- LFinfo$plba
      plb <- LFinfo$plb
      page <- LFinfo$page
      LF <- LFinfo$LF
      LF0 <- LF0info$LF

      if (pool == TRUE) {
        LF_t <- LF0_t <- matrix(NA, nrow = tyears_only, ncol = ncol(LF))
        for (y in 1:tyears_only) {
          if (nseasons == 1) {
            LF_t[y,] <- LF[y,]
            LF0_t[y,] <- LF0[y,]
          }
          if (nseasons > 1) {
            time_index <- (1:nseasons) + ((y - 1) * nseasons)
            LF_t[y,] <- colSums(LF[time_index,])
            LF0_t[y,] <- colSums(LF0[time_index,])
          }
        }
        obs_per_year <- sapply(1:tyears_only, function(x) {
          if (nseasons == 1)
            time_index <- x
          if (nseasons > 1)
            time_index <- (1:nseasons) + ((x - 1) * nseasons)
          sum(obs_per_year[time_index])
        })
      }
      if (pool == FALSE) {
        LF_t <- LF
        LF0_t <- LF0
      }



      ########################################################
      ## Expected mean length in catch
      ########################################################
      ML_t <- vector(length = tyears)
      for (y in 1:tyears) {
        vul_pop <- sum(N_at[, y] * S_a)
        vul_lengths <- sum(vul_pop * plb[y,] * mids)
        ML_t[y] <- vul_lengths / vul_pop
      }
      if (pool == TRUE)
        ML_t <- ML_t[which(1:tyears %% nseasons == 0)]

      ########################################################
      ## cut out burn-in
      ########################################################

      if (pool == TRUE) {
        LFout <- LF_t[-c(1:nburn_real),]
        rownames(LFout) <- 1:Nyears_real
        LF0out <- LF0_t[-c(1:nburn_real),]
        rownames(LF0out) <- 1:Nyears_real

        ML_tout <- ML_t[-c(1:nburn_real)]

        LFindex <- (Nyears_real - Nyears_comp + 1):Nyears_real
      }
      if (pool == FALSE) {
        LFout <- LF_t[-c(1:nburn),]
        rownames(LFout) <- 1:Nyears
        LF0out <- LF0_t[-c(1:nburn),]
        rownames(LF0out) <- 1:Nyears

        ML_tout <- ML_t[-c(1:nburn)]

        LFindex <- (Nyears - Nyears_comp * nseasons + 1):Nyears
      }

      I_tout <- I_t[-c(1:nburn_real)]
      C_tout <- C_t[-c(1:nburn_real)]
      Cw_tout <- Cw_t[-c(1:nburn_real)]
      names(C_tout) <-
        names(Cw_tout) <- names(I_tout) <- 1:Nyears_real
      R_tout <- R_t[-c(1:nburn_real)]
      N_tout <- N_t[-c(1:nburn_real)]
      SB_tout <- SB_t[-c(1:nburn_real)]
      TB_tout <- TB_t[-c(1:nburn_real)]
      VB_tout <- VB_t[-c(1:nburn_real)]
      D_tout <- D_t[-c(1:nburn_real)]
      F_tout <- F_t[-c(1:nburn_real)]
      SPR_tout <- SPR_t[-c(1:nburn_real)]
      Z_tout <- Z_t[-c(1:nburn_real)]

      LFout <- LFout[LFindex,]
      LF0out <- LF0out[LFindex,]
      if (is.vector(LFout) == FALSE)
        colnames(LFout) <- highs
      if (is.vector(LF0out) == FALSE)
        colnames(LF0out) <- highs
      if (is.vector(LFout)) {
        LFout <- t(as.matrix(LFout))
        rownames(LFout) <- LFindex
      }
      if (is.vector(LF0out)) {
        LF0out <- t(as.matrix(LF0out))
        rownames(LF0out) <- LFindex
      }

      if (mismatch == TRUE)
        myrs <- 1:LFindex[1]
      if (mismatch == FALSE)
        myrs <- 1:max(LFindex)

      ## outputs
      lh$I_t <- I_tout[myrs]
      lh$C_t <- C_tout[myrs]
      lh$Cw_t <- Cw_tout[myrs]
      lh$LF <- LFout
      lh$LF0 <- LF0out
      lh$R_t <- R_tout
      lh$N_t <- N_tout
      lh$SB_t <- SB_tout
      lh$D_t <- D_tout
      lh$F_t <- F_tout
      lh$ML_t <- ML_tout
      lh$plb <- plb
      lh$plba <- plba
      lh$page <- page
      lh$N_at <- N_at
      lh$SPR <- SPR
      lh$SPR_t <- SPR_tout
      lh$SPR_alt <- SPR_alt
      lh$VB_t <- VB_tout
      lh$TB_t <- TB_tout
      lh$nlbins <- length(mids)
      lh$Z_t <- Z_tout
      if (pool == TRUE) {
        lh$Nyears <- Nyears_real
        lh$years <- 1:Nyears_real
      }
      if (pool == FALSE) {
        lh$Nyears <- Nyears
        lh$years <- 1:Nyears
      }
      lh$obs_per_year <- obs_per_year
      if (Rdynamics != "AR")
        lh$RecDev <- RecDev
      if (Rdynamics == "AR")
        lh$RecDev <- RecDev_AR
      lh$FishDev <- FishDev

      return(lh)

    }) ## end with function

  }
