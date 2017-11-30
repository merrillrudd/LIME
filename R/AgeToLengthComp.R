#' Age-converted-to-length structure
#'
#' \code{AgeToLengthComp} Converts vulnerable numbers at age to length composition of catch
#' @author M.B. Rudd

#' @param lh list of life history attributes, output of create_lh_list
#' @param tyears number of years of data
#' @param N_at matrix of numbers in the population at each age over time
#' @param comp_sample vector of number of individuals sampled each year (set as 1 for proportions)
#' @param sample_type a character vector specifying if the length comps are sampled from the 'catch' (default) or from the population
#' @importFrom stats pnorm rmultinom
#' 
#' @return data frame - number of individuals in each length bin in each year
#' @export
AgeToLengthComp <-
  function(lh,
           tyears,
           N_at,
           comp_sample,
           sample_type = 'catch') {
    with(lh, {
      ################################################
      ## Probability being in a length bin given age
      ################################################
      lbprobs <-
        function(mnl, sdl)
          return(pnorm(highs, mnl, sdl) - pnorm(lows, mnl, sdl))
      vlprobs <- Vectorize(lbprobs, vectorize.args = c("mnl", "sdl"))
      plba <- t(vlprobs(L_a, L_a * CVlen))
      plba <- plba / rowSums(plba)

      ################################################
      ## Probability being in harvested at an age
      ################################################
      page <- matrix(ncol = dim(plba)[1], nrow = tyears)
      for (y in 1:tyears)
        # page[y, ] <- N_at[, y] * ifelse(sample_type == 'catch', S_a, 1)
        if(sample_type=="catch") page[y,] <- N_at[,y] * S_a
        if(sample_type!="catch") page[y,] <- N_at[,y]
      page <- page / rowSums(page)

      ################################################
      ## Probability of sampling a given length bin
      ################################################
      plb <- matrix(ncol = length(highs), nrow = tyears)
      for (y in 1:tyears)
        plb[y, ] <- page[y, ] %*% plba
      plb <- plb / rowSums(plb)

      #######################
      ## Length frequencies
      #######################
      LF <- array(0, dim = dim(plb))
      rownames(LF) <- 1:tyears
      for (y in 1:tyears) {
        if(is.na(sum(plb[y,]))==FALSE){
          LF[y, ] <- rmultinom(n = 1,
                             size = comp_sample[y],
                             prob = plb[y, ])
        }
        if(is.na(sum(plb[y,]))){
          LF[y,] <- NA
        }
      }

      Outs <- NULL
      Outs$plba <- plba
      Outs$plb <- plb
      Outs$page <- page
      Outs$LF <- LF
      return(Outs)
    })


  }