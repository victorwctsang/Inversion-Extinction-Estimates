
library(readxl)

getFossilData = function (path, sheet, range, col_names, col_types) {
  read_excel(path, sheet, range, col_names, col_types)
}

boundFossilData = function (fossil.df, K) {
  return(fossil.df[fossil.df$age < K, ])
}

getStdDevFromFossilData = function (path, K, ...) {
  df = getFossilData(path, ...)
  df = boundFossilData(df, K)
  return(mean(df$sd))
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
  
  dataset.df = data.frame(sim_id=rep(1:n.sims, each=n.error_factors),
                          error_factor=rep(error_factor, length.out=df.length))
  W = lapply(dataset.df$error_factor, function(x) simulate_dataset(x, theta.true, K, eps.mean, eps.sigma, n))
  dataset.df$W = W
  return(dataset.df)
}

