# GRIWM implementation
# Adapted from https://github.com/sarakahanamoku/GRIWM/blob/master/GRIWM.R

library(dplyr)

GRIWM <- function(df, alpha, .iter=10000){
  
  dates <- df %>% dplyr::pull(dates)
  sd <- df %>% dplyr::pull(sd)
  k <- length(df$dates) # original returns 1
  iter <- .iter
  
  # initialize vectors
  w_time_mci <- rep(0,iter)
  
  for (c in 1:iter) {
    date_samp <- rep(0, k)
    ## resampling of the standard deviation of each date from a Gaussian distribution
    for (ii in 1:k) {
      
      date_samp[ii] <-
        round(rnorm(1, mean = as.numeric(dates[ii]), sd = as.numeric(sd[ii])))
    }
    date_samp <- (sort(date_samp))
    
    ## down-wighting procedure: calculate weighted McInerney date & confidence interval
    last.diff <- 1 / (date_samp - date_samp[1])[-1]
    weight <- last.diff / last.diff[1]
    
    if (last.diff[1] == Inf) {
      weight <- last.diff / last.diff[2]
      weight <- weight[-1]
    }
    
    ldate <- length(date_samp)
    T.mci.lst.vec <- rep(0, ldate - 1)
    
    for (m in 1:(ldate - 1)) {
      date.it <- date_samp[1:(1 + m)]
      date.age.it <- date_samp[1:(1 + m)]
      
      date.mci.it <- rev(max(date.it) + 1 - date.it)
      
      k <- length(date.it)
      t.n <- date.mci.it[k]
      n <- k
      T.rng <- t.n - date.mci.it[1]
      
      #probability of finding another record, estimated from the previous sighting rate and the time
      #since the last observation
      i <- t.n
      p.iter <- 1
      while (p.iter > alpha)
      {
        i <- i + 1
        p.iter <- (1 - (n / t.n)) ^ (i - t.n)
      }
      
      T.mci.lst.vec[m] <- max(date.it) + 1 - i
    }
    
    if (last.diff[1] == Inf) {
      w_time_mci[c] <-
        round((sum(weight * T.mci.lst.vec[-1])) / sum(weight), 0)
    }
    
    if (last.diff[1] != Inf) {
      w_time_mci[c] <- round((sum(weight * T.mci.lst.vec)) / sum(weight), 0)
    }
    
    
    # if (c %% iterdiv == 0) # this "progress output" slows down the code a lot!
    #   print(c)
  }
  
  #Confidence interval setting
  lower_ci <- round(quantile(na.omit(w_time_mci), probs = (alpha / 2)), 0)
  
  # Centroid
  centroid <- round(median(na.omit(w_time_mci)), 0)
  
  # Upper boundary of CI
  upper_ci <- round(quantile(na.omit(w_time_mci), probs = (1 - alpha / 2)), 0)
  
  estimate <- tibble(upper_ci, centroid, lower_ci)
  return(estimate)
}