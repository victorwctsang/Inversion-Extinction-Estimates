
library(readxl)
source("minmi-functions.R")

getFossilData = function (path, sheet, range, col_names, col_types) {
  read_excel(path, sheet, range, col_names, col_types)
}

boundFossilData = function (fossil.df, K) {
  return(fossil.df[fossil.df$age < K, ])
}

getStdDevFromFossilData = function (path, K, ...) {
  df = getFossilData(path, ...)
  df = boundFossilData(df, K)
  return(df$sd)
}

create_result_df = function (n.sims, method, error_factor) {
  n.methods = length(method)
  n.error_factors = length(error_factor)
  df.length = n.sims * n.methods * n.error_factors
  
  data.frame(sim_id=rep(1:n.sims, each=(n.methods*n.error_factors)),
             method=rep(method, each=n.error_factors, length.out=df.length),
             error_factor=rep(error_factor, length.out=df.length),
             lower=rep(NA, df.length),
             upper=rep(NA, df.length),
             mean=rep(NA, df.length),
             width=rep(NA, df.length),
             contains_theta=rep(NA, df.length),
             compute_time=rep(NA, df.length))
}

simulate_dataset = function (error_factor, theta.true, K, eps.mean, eps.sigma, n) {
  X = runif(n, min=theta.true, max=K)
  eps = rnorm(n, mean=eps.mean, sd=error_factor*eps.sigma)
  W = X + eps
  return(W)
}

simulate_datasets = function (n.sims, error_factor, theta.true, K, eps.mean, eps.sigma, n) {
  n.error_factors = length(error_factor)
  df.length = n.sims * n.error_factors
  
  dataset.df = data.frame(id=1:df.length,
                          error_factor=rep(error_factor, length.out=df.length))
  W = lapply(dataset.df$error_factor, function(x) simulate_dataset(x, theta.true, K, eps.mean, eps.sigma, n))
  dataset.df$W = W
  return(dataset.df)
}

estimate_quantile = function (W, sd, method, q, K, dating_error.mean) {
  estimate=NA
  runtime = NA
  start_time=Sys.time()
  estimate = switch(method,
    MLE = if (q != 0.5) NULL else mle(W),
    `BA-MLE` = if (q != 0.5) NULL else ba_mle(W, K),
    STRAUSS = if (q != 0.5) NULL else strauss(W, K),
    MINMI = minmi(q, W, sd, K, dating_error.mean)
  )
  runtime = as.numeric(difftime(Sys.time(), start_time), units="secs")
  return(list(q=q, estimate=estimate, runtime=runtime))
}

mle = function (W) {
  return(min(W))
}

ba_mle = function (W, K) {
  n = length(W)
  return(min(W) * (n+1)/n - K/n)
}

strauss = function (W) {
  n = length(W)
  return((n*min(W) - max(W))/(n-1))
}

minmi = function (q, W, sd, K, dating_error.mean=0, .B_init=500) {
  m = min(W)
  n = length(W)
  dating_error.sd = mean(sd)
  
  estimate = estimate_quantile.minmi(q=q, K=K, W=W)
  
  if (dating_error.sd != 0) {
    # Generate initial MC samples
    u.init = runif(.B_init, 0, 1)
    
    # Estimate optimal value for B
    B = find_optimal_B(max_var=(0.2*dating_error.sd)^2, q=q, K=K, m=m, u=u.init, eps.mean=dating_error.mean, eps.sigma=dating_error.sd)
    
    # Generate MC samples
    mc.samples = runif(B, 0, 1)
    estimate = estimate_quantile.minmi(q=q, K=K, W=W, u=mc.samples, eps.mean=dating_error.mean, eps.sigma=dating_error.sd)
  }
  return(estimate)
}