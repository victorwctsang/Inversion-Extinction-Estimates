# GRIWM implementation
# Adapted from https://github.com/sarakahanamoku/GRIWM/blob/master/GRIWM.R

library(dplyr)

GRIWM = function(df, alpha, K, .n_iter=10000){
  
  # initialize vectors
  dates = df %>% dplyr::pull(dates)
  sd = df %>% dplyr::pull(sd)
  n = length(df$dates)
  estimates.griwm = rep(0,.n_iter)
  
  for (iter in 1:.n_iter) {
    date_samp = rep(0, n)
    ## resampling of the standard deviation of each date from a Gaussian distribution
    for (i in 1:n) {
      date_samp[i] = round(rnorm(1, mean = as.numeric(dates[i]), sd = as.numeric(sd[i])))
    }
    date_samp = (sort(date_samp))
    
    ## down-weighting procedure
    
    ### Calculate weights
    last.diff = 1 / (date_samp - date_samp[1])[-1]
    weight = last.diff / last.diff[1]
    
    if (last.diff[1] == Inf) {
      weight = last.diff / last.diff[2]
      weight = weight[-1]
    }
    
    ### Calculate McInerny dates
    estimates.mcinerny = rep(0, n-1)
    lambda = (n-1)/(K-date_samp[1])
    
    for (k in 2:n) {
      date.recent_k = date_samp[1:k]
      estimates.mcinerny[k-1] = date_samp[1] - log(alpha)/log(1-lambda)
    }
    
    # Get weighted extinction date
    if (last.diff[1] == Inf) {
      estimates.griwm[iter] =
        round((sum(weight * estimates.mcinerny)) / sum(weight), 0)
    }
    
    if (last.diff[1] != Inf) {
      estimates.griwm[iter] = round((sum(weight * estimates.mcinerny)) / sum(weight), 0)
    }
  }
  
  #Confidence interval setting
  lower_ci = round(quantile(na.omit(estimates.griwm), probs = (alpha / 2)), 0)
  
  # Centroid
  centroid = round(median(na.omit(estimates.griwm)), 0)
  
  # Upper boundary of CI
  upper_ci = round(quantile(na.omit(estimates.griwm), probs = (1 - alpha / 2)), 0)
  
  estimate = tibble(lower_ci, centroid, upper_ci)
  return(estimate)
}