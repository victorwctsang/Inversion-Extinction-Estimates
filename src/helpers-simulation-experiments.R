
library(readxl)
source("../src/minmi-functions.R")
source("../src/GRIWM.R")
source("../src/simulated-inversion.R")
source("../src/garthwaite-robbins-munro-functions.R")

getFossilData = function (path, sheet, range, col_names, col_types) {
  read_excel(path, sheet, range, col_names, col_types)
}

boundFossilData = function (fossil.df, K) {
  return(fossil.df[fossil.df$age < K,])
}

getStdDevFromFossilData = function (path, K, n, ...) {
  df = getFossilData(path, ...)
  df = boundFossilData(df, K)
  return(df$sd[1:n])
}

create_result_df = function (n.sims, method, error_factor) {
  n.methods = length(method)
  n.error_factors = length(error_factor)
  df.length = n.sims * n.methods * n.error_factors
  
  data.frame(
    sim_id = rep(1:n.sims, each = (n.methods * n.error_factors)),
    method = rep(method, each = n.error_factors, length.out = df.length),
    error_factor = rep(error_factor, length.out = df.length),
    lower = rep(NA, df.length),
    upper = rep(NA, df.length),
    mean = rep(NA, df.length),
    width = rep(NA, df.length),
    contains_theta = rep(NA, df.length),
    compute_time = rep(NA, df.length)
  )
}

simulate_dataset = function (error_factor,
                             theta.true,
                             K,
                             eps.mean,
                             eps.sigma,
                             n) {
  eps = rnorm(n, mean = eps.mean, sd = error_factor * eps.sigma)
  X = runif(n, min = theta.true, max = K)
  W = X + eps
  return(W)
}

simulate_datasets = function (n.sims,
                              error_factor,
                              theta.true,
                              K,
                              eps.mean,
                              eps.sigma,
                              n) {
  n.error_factors = length(error_factor)
  df.length = n.sims * n.error_factors
  
  dataset.df = data.frame(
    id = 1:df.length,
    error_factor = rep(error_factor, length.out =
                         df.length)
  )
  W = lapply(dataset.df$error_factor, function(x)
    simulate_dataset(x, theta.true, K, eps.mean, eps.sigma, n))
  dataset.df$W = W
  return(dataset.df)
}

estimate_extinction = function (W, sd, method, K, dating_error.mean) {
  estimate = NA
  runtime = NA
  start_time = Sys.time()
  estimate = switch(
    method,
    MLE = mle(W),
    `BA-MLE` = ba_mle(W, K),
    Strauss = strauss(W),
    NULL
  )
  runtime = calculate_tdiff(start_time, Sys.time())
  return(list(point = estimate, point_runtime = runtime))
}

calculate_tdiff = function(start, end) {
  as.numeric(difftime(end, start), units = "secs")
}

estimate_conf_int = function (W,
                              sd,
                              method,
                              alpha,
                              K,
                              dating_error.mean,
                              B.point,
                              B.lower,
                              B.upper,
                              theta.test_vec) {
  estimate = switch(
    method,
    GRIWM = griwm(
      alpha = alpha,
      dates = W,
      sd = sd,
      K = K,
      p_t = alpha,
      bias_adjusted = F
    ),
    `GRIWM-corrected` = griwm(
      alpha = alpha,
      dates = W,
      sd = sd,
      K = K,
      p_t = 0.5,
      bias_adjusted = T
    ),
    MINMI = minmi(
      alpha = alpha,
      W = W,
      sd = sd,
      K = K,
      dating_error.mean = dating_error.mean,
      B.point = B.point,
      B.lower = B.lower,
      B.upper = B.upper
    ),
    `SI-UGM` = SI_UGM(
      alpha = alpha,
      dates = W,
      sd = sd,
      K = K,
      dating_error.mean = dating_error.mean,
      theta.test_vec = theta.test_vec
    ),
    `SI-RM` = SI_RM_process(
      alpha = alpha,
      dates = W,
      sd = sd,
      K = K,
      dating_error.mean = dating_error.mean,
      max_var = NULL
    ),
    `SI-RM-1000` = SI_RM_process(
      alpha = alpha,
      dates = W,
      sd = sd,
      K = K,
      dating_error.mean = dating_error.mean,
      max_var = -1
    )
  )
  return(estimate)
}

mle = function (W) {
  return(min(W))
}

ba_mle = function (W, K) {
  n = length(W)
  return(min(W) * (n + 1) / n - K / n)
}

strauss = function (W) {
  n = length(W)
  return((n * min(W) - max(W)) / (n - 1))
}

minmi = function (alpha,
                  W,
                  sd,
                  K,
                  dating_error.mean = 0,
                  .B_init = 500,
                  B.point,
                  B.lower,
                  B.upper) {
  m = min(W)
  n = length(W)
  dating_error.sd = mean(sd)
  
  lower = NULL
  point = NULL
  upper = NULL
  point_runtime = NULL
  conf_int_runtime = NULL
  
  if (dating_error.sd == 0) {
    # Confidence Interval
    conf_int.start_time = Sys.time()
    lower = estimate_quantile.minmi(q = alpha / 2, K = K, W = W)
    upper = estimate_quantile.minmi(q = 1 - alpha / 2,
                                    K = K,
                                    W = W)
    conf_int_runtime = calculate_tdiff(conf_int.start_time, Sys.time())
    
    # Point estimate
    point.start_time = Sys.time()
    point = estimate_quantile.minmi(q = 0.5, K = K, W = W)
    point_runtime = calculate_tdiff(point.start_time, Sys.time())
    
  } else {
    B.max = max(B.point, B.lower, B.upper)
    
    # Confidence interval
    conf_int.start_time = Sys.time()
    mc.samples = runif(B.max, 0, 1)
    lower = estimate_quantile.minmi(
      q = alpha / 2,
      K = K,
      W = W,
      u = mc.samples[1:B.lower],
      eps.mean = dating_error.mean,
      eps.sigma = dating_error.sd
    )
    upper = estimate_quantile.minmi(
      q = 1 - alpha / 2,
      K = K,
      W = W,
      u = mc.samples[1:B.upper],
      eps.mean = dating_error.mean,
      eps.sigma = dating_error.sd
    )
    conf_int_runtime = calculate_tdiff(conf_int.start_time, Sys.time())
    
    # Point estimate
    point.start_time = Sys.time()
    point = estimate_quantile.minmi(
      q = 0.5,
      K = K,
      W = W,
      u = mc.samples[1:B.point],
      eps.mean = dating_error.mean,
      eps.sigma = dating_error.sd
    )
    point_runtime = calculate_tdiff(point.start_time, Sys.time())
  }
  return(
    list(
      lower = lower,
      point = point,
      upper = upper,
      point_runtime = point_runtime,
      conf_int_runtime = conf_int_runtime,
      B.point = B.point,
      B.lower = B.lower,
      B.upper = B.upper
    )
  )
}

griwm = function(alpha,
                 dates,
                 sd,
                 K,
                 p_t,
                 .iter = 10000,
                 bias_adjusted) {
  df = data.frame(dates = dates, sd = sd)
  start_time = Sys.time()
  results = GRIWM(df, alpha, K, p_t, .iter, bias_adjusted = bias_adjusted)
  runtime = calculate_tdiff(start_time, Sys.time())
  return(
    list(
      lower = as.numeric(results$lower_ci),
      point = as.numeric(results$centroid),
      upper = as.numeric(results$upper_ci),
      point_runtime = runtime,
      conf_int_runtime = runtime
    )
  )
}

SI_UGM = function (alpha,
                   dates,
                   sd,
                   K,
                   dating_error.mean = 0,
                   theta.test_vec) {
  start_time = Sys.time()
  results = simulated_inversion(alpha, dates, sd, K, theta.test_vec, dating_error.mean = 0)
  runtime = calculate_tdiff(start_time, Sys.time())
  return(
    list(
      lower = results$lower,
      upper = results$upper,
      point = results$point,
      point_runtime = runtime,
      conf_int_runtime = runtime
    )
  )
}

SI_RM_process = function (alpha,
                          dates,
                          sd,
                          K,
                          dating_error.mean = 0,
                          max_var = NULL) {
  # use SI-UGM to estimate a model to get the gradient.
  SI_results = simulated_inversion(
    alpha,
    dates,
    sd,
    K,
    theta.test_vec,
    dating_error.mean = 0,
    return_model = T
  )
  
  CI_start_time = Sys.time()
  CI = estimate_CI.rm(
    W = dates,
    K = K,
    alpha = alpha,
    max_iter = 1000,
    eps.mean = dating_error.mean,
    eps.sigma = mean(sd),
    .model = SI_results$model,
    .CI_estimates = SI_results,
    return_iters = T,
    .starting_vals = NULL,
    max_var = max_var
  )
  conf_int_runtime = calculate_tdiff(CI_start_time, Sys.time())
  
  return(
    list(
      lower = CI$lower,
      upper = CI$upper,
      point = CI$point,
      point_runtime = conf_int_runtime,
      conf_int_runtime = conf_int_runtime,
      num_iters.lower = length(CI$lower.iters),
      num_iters.upper = length(CI$upper.iters)
    )
  )
}