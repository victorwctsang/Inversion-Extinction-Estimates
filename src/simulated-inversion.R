
library(mgcv)

simulated_inversion = function (alpha, dates, sd, K, theta.test_vec=NULL, dating_error.mean=0, return_model=F) {
  n.fossils = length(dates)
  n.theta.test_vec = length(theta.test_vec)
  dating_error.sd = mean(sd)
  
  flag.return_thetas = F
  if (is.null(theta.test_vec)) {
    theta.test_vec = seq(from = ba_mle(dates, K) - 10*sd(dates),
                         to = min(dates) +2*sd(dates),
                         length.out = 5000)
    flag.return_thetas = T
  }
  
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
  results = list(lower=lower, upper=upper, point=point)
  
  # Return results
  if (return_model) {
    results$model = model_gam
  }
  if (flag.return_thetas) {
    results$theta.test_vec = theta.test_vec
  }
  return(results)
}

get_min_est <- function (theta.test_vec , model_gam, q){
  prediction <- function (theta.test_vec , model_gam ){
    yval <- predict ( model_gam , newdata = data.frame ( theta.test_vec = theta.test_vec ), type = "response")
    return ( yval - q )
  }
  theta_est <- uniroot ( prediction , lower = min(theta.test_vec) , upper = max(theta.test_vec), extendInt = "yes", model_gam = model_gam, maxiter = 2000)$root
  return ( theta_est )
}

ba_mle = function (X, K) {
  n.X = length(X)
  return(min(X) * (n.X + 1) / n.X - K / n.X)
}