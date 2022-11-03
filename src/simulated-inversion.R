
library(mgcv)

simulated_inversion = function (alpha, dates, sd, K, theta.test_vec, dating_error.mean=0) {
  n.fossils = length(dates)
  n.theta.test_vec = length(theta.test_vec)
  dating_error.sd = mean(sd)
  
  # Calculate S_star
  S_star = rep(NA, n.theta.test_vec)
  for (i in 1:n.theta.test_vec) {
    fossil.resamples = runif(n.fossils, theta.test_vec[i], K) + rnorm(n.fossils, dating_error.mean, dating_error.sd)
    S_star[i] = min(fossil.resamples)
  }
  
  # Create indicators
  m = min(dates)
  indicator_vec = sapply(S_star, function(x) ifelse(x > m, 1, 0))
  
  # Estimate model using smooth regression
  model_gam = gam(indicator_vec ~ s(theta.test_vec),
                  data=data.frame(theta.test_vec, indicator_vec),
                  family=binomial)

  # Use inversion
  lower = get_min_est(theta.test_vec, model_gam, q=alpha/2)
  upper = get_min_est(theta.test_vec, model_gam, q=1-alpha/2)
  point = get_min_est(theta.test_vec, model_gam, q=0.5)
  if (abs(upper - lower) > 1000000) {
    browser()
  }
  return(list(lower=lower, upper=upper, point=point))
}

get_min_est <- function (theta.test_vec , model_gam, q){
  prediction <- function (theta.test_vec , model_gam ){
    yval <- predict ( model_gam , newdata = data.frame ( theta.test_vec = theta.test_vec ), type = "response")
    return ( yval - q )
  }
  theta_est <- uniroot ( prediction , lower = 5000 , upper = 15000, extendInt = "yes", model_gam = model_gam, maxiter = 2000)$root
  return ( theta_est )
}