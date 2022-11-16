

estimate_CI.rm = function (W,
                           K,
                           alpha,
                           max_iter,
                           eps.mean,
                           eps.sigma,
                           .model,
                           .CI_estimates,
                           return_iters = FALSE,
                           .starting_vals = NULL,
                           max_var = NULL) {
  n = length(W)
  # get point estimate of theta
  theta.hat = estimate_theta.rm(W)
  # calculate proportionality constant
  p = calc_prop_constant(alpha)
  
  # Get estimates of gradients at CI end points
  g = estimate_g(alpha, .model, .CI_estimates)
  
  # calculate m
  m = ceiling(min(50, 0.3 * (2 - alpha) / alpha))
  
  # Adjust number of iterations for m
  max_iter = max_iter + m - 1
  
  # initialize lower and upper vecs
  lower.iters = upper.iters = rep(NA, max_iter)
  
  # compute starting estimates
  if (is.null(.starting_vals)) {
    starting_ests = get_starting_vals(alpha,
                                      theta.hat,
                                      n,
                                      K,
                                      eps.mean,
                                      eps.sigma)
    lower.iters[m] = starting_ests[, "lower"]
    upper.iters[m] = starting_ests[, "upper"]
  } else {
    lower.iters[m] = .starting_vals[1]
    upper.iters[m] = .starting_vals[2]
  }
  
  # Apply RM process
  if (is.null(max_var)) {
    max_var = 0.2 * (eps.sigma) ^ 2
  }
  
  lower.iters = estimate_bound.rm(
    theta.iters = lower.iters,
    q = alpha / 2,
    theta.hat = theta.hat,
    n = n,
    K = K,
    p = p,
    m = m,
    g = g$lower,
    max_var = max_var,
    max_iter = max_iter,
    eps.mean = eps.mean,
    eps.sigma = eps.sigma
  )
  upper.iters = estimate_bound.rm(
    theta.iters = upper.iters,
    q = 1 - alpha / 2,
    theta.hat = theta.hat,
    n = n,
    K = K,
    p = p,
    m = m,
    g = g$upper,
    max_var = max_var,
    max_iter = max_iter,
    eps.mean = eps.mean,
    eps.sigma = eps.sigma
  )
  lower.iters = na.omit(lower.iters)
  upper.iters = na.omit(upper.iters)
  
  CI = list(lower = lower.iters[length(lower.iters)],
            point = theta.hat,
            upper = upper.iters[length(upper.iters)])
  if (return_iters == TRUE) {
    CI$lower.iters = lower.iters
    CI$upper.iters = upper.iters
  }
  return(CI)
}

estimate_theta.rm = function (W) {
  min(W)
}

calc_prop_constant = function (alpha) {
  # p := Step length proportionality constant
  2 / (qnorm(1 - alpha) * (2 * pi) ^ (-0.5) * exp(-(qnorm(1 - alpha) ^ 2) /
                                                    2))
}

get_starting_vals = function (alpha, theta, n, K, eps.mean, eps.sigma) {
  # Implements percentile method (Buckland, 1980; Efron, 1981)
  n.resamples = (2 - alpha) / alpha
  fossil.resamples = matrix(simulate_fossils(n * n.resamples,
                                             theta,
                                             K,
                                             eps.mean,
                                             eps.sigma),
                            ncol = n.resamples)
  theta.resamples = apply(fossil.resamples, 2, FUN = min)
  theta.resamples.sorted = sort(theta.resamples)
  return(
    cbind(
      lower = theta.resamples.sorted[2],
      point = median(theta.resamples.sorted),
      upper = theta.resamples.sorted[n.resamples - 1]
    )
  )
}

simulate_fossils = function (n,
                             theta,
                             K,
                             eps.mean = 0,
                             eps.sigma = 0) {
  # Simulate fossils assuming:
  # - Uniform deposition from theta to K
  # - Gaussian measurement error
  X = runif(n, min = theta, max = K)
  eps = rnorm(n, mean = eps.mean, sd = eps.sigma)
  W = X + eps
  return(W)
}

estimate_bound.rm = function (theta.iters,
                              q,
                              theta.hat,
                              n,
                              K,
                              p,
                              m,
                              g,
                              max_var,
                              max_iter,
                              eps.mean,
                              eps.sigma) {
  # Transformation to enforce theta < K
  eta.fun = function(theta)
    - log(K - theta) #log(K/(K-theta))
  theta.fun = function(eta)
    K - exp(-eta) #K*(1-exp(-et))
  
  eta_q.iters = eta.fun(theta.iters)
  
  # Start RM process
  i = m
  mc.var = Inf
  while (i < max_iter & mc.var > max_var) {
    # Generate resamples
    theta_q.hat = theta.fun(eta_q.iters[i])
    resamples = simulate_fossils(n,
                                 theta = theta_q.hat,
                                 K,
                                 eps.mean,
                                 eps.sigma)
    theta.hat.resample = estimate_theta.rm(resamples)
    
    # calculate step length
    c.eta = 2 * p * abs(eta.fun(theta_q.hat) - eta.fun(theta.hat))
    
    # Update
    if (theta.hat.resample <= theta.hat) {
      eta_q.iters[i + 1] = (eta_q.iters[i] + c.eta * q / i)
    } else {
      eta_q.iters[i + 1] = (eta_q.iters[i] - c.eta * (1 - q) / i)
    }
    
    # Calculate RM variance
    mc.var = calculate_rm_var(q, c.eta, i, g * exp(-eta_q.iters[i]))
    mc.var = mc.var / (K - theta.fun(eta_q.iters[i])) ^ 2
    
    # Iterate
    i = i + 1
  }
  theta_q.iters = theta.fun(eta_q.iters)
  return(theta_q.iters)
}

estimate_g = function (alpha, model, CI.estimates, .delta = 0.0001) {
  g = list(lower = NA, upper = NA)
  
  perturbs.lower = predict(model, newdata = data.frame(
    theta.test_vec = c(CI.estimates$lower - .delta, CI.estimates$lower + .delta)
  ))
  g$lower = (perturbs.lower[2] - perturbs.lower[1]) / (2 * .delta)
  
  perturbs.upper = predict(model, newdata = data.frame(
    theta.test_vec = c(CI.estimates$upper - .delta, CI.estimates$upper + .delta)
  ))
  g$upper = (perturbs.upper[2] - perturbs.upper[1]) / (2 * .delta)
  return(g)
}

calculate_rm_var = function (q, c, i, g) {
  rm_var = Inf
  if (g * c > 0.5) {
    rm_var = q * (1 - q) * c ^ 2 / (i * (2 * g * c - 1))
  }
  return(rm_var)
}