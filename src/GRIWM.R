# GRIWM implementation
# Adapted from https://github.com/sarakahanamoku/GRIWM/blob/master/GRIWM.R

GRIWM = function(df, alpha, K, p_t, .n_iter=10000, bias_adjusted=F){
  
  # initialize vectors
  dates = dplyr::pull(df, dates)
  sd = dplyr::pull(df, sd)
  n = length(df$dates)
  estimates.griwm = rep(0,.n_iter)
  
  for (iter in 1:.n_iter) {
    dates.resampled = rep(0, n)
    ## Gaussian resampling step
    for (i in 1:n) {
      dates.resampled[i] = dates[i] + rnorm(1, 0, mean(sd))
    }
    dates.resampled = (sort(dates.resampled))
    
    ## down-weighting procedure
    
    ### Calculate weights
    last.diff = 1 / (dates.resampled - dates.resampled[1])[-1]
    weight = last.diff / last.diff[1]
    
    if (last.diff[1] == Inf) {
      weight = last.diff / last.diff[2]
      weight = weight[-1]
    }
    
    ### Calculate McInerny dates
    estimates.mcinerny = rep(0, n-1)
    
    lambda = n/(K-dates.resampled[1])
    
    for (k in 2:n) {
      dates.recent_k = dates.resampled[1:k]
      theta = dates.recent_k[1] - 1
      if (bias_adjusted == T) {
        theta = uniroot(f=function(theta) {(1 - n/(K-theta))^(dates.recent_k[1] - theta) - p_t}, interval=c(min(dates) - (K-min(dates))/2, min(dates) + (K-min(dates))/2), extendInt="yes")$root
      } else {
        theta = dates.recent_k[1] - log(p_t)/log(1 - n/(K-dates.recent_k[1]))
      }
      estimates.mcinerny[k-1] = theta
    }

    # Get weighted extinction date
    if (last.diff[1] == Inf) {
      estimates.griwm[iter] = round((sum(weight * estimates.mcinerny[-1])) / sum(weight), 0)
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
  
  estimate = dplyr::tibble(lower_ci, centroid, upper_ci)
  return(estimate)
}