
estimate_CI.rm = function (W, K, alpha, max_iter, eps.mean, eps.sigma, return_iters=FALSE) {
  # estimate confidence interval using Garthwaite's Robbins-Munro process (1995)
  n=length(W)
  # get point estimate of theta
  theta.hat = estimate_theta.rm(W)
  # calculate p
  p = calc_prop_constant(alpha)
  # initialize lower and upper vecs
  lower.iters = upper.iters = rep(NA, max_iter)
  # calculate m
  m = ceiling(min(50, 0.3 * (2-alpha)/alpha))
  # compute starting estimates
  starting_ests = get_starting_vals(method="percentile", alpha, theta.hat, n, K, eps.mean, eps.sigma)
  lower.iters[m] = starting_ests[, "lower"]
  upper.iters[m] = starting_ests[, "upper"]
  # estimate lower
  lower.iters = estimate_bound.rm(lower.iters, alpha/2, theta.hat, n, K, p, m, max_var=0.2*eps.sigma, max_iter, eps.mean, eps.sigma)
  upper.iters = estimate_bound.rm(upper.iters, 1-alpha/2, theta.hat, n, K, p, m, max_var=0.2*eps.sigma, max_iter, eps.mean, eps.sigma)
  CI = list(CI.lower=lower.iters[length(lower.iters)], CI.upper=upper.iters[length(upper.iters)])
  if (return_iters == TRUE) {
    CI$CI.iters = cbind(lower=lower.iters, upper=upper.iters)
  }
  return(CI)
}

estimate_theta.rm = function (W) {
  min(W)
}

calc_prop_constant = function (alpha) {
  # p := Step length proportionality constant (garthwaite calls it k)
  # Garthwaite 1995 paper
  2 / (qnorm(1-alpha) * (2*pi)^(-0.5) * exp(- (qnorm(1-alpha)^2)/2))   # Normal distribution heuristic
}

get_starting_vals = function (method="percentile", ...) {
  switch(method,
         percentile=get_starting_est.percentile(...),
         analytic=get_starting_vals.analytic(...))
}

simulate_fossils = function (n, theta, K, eps.mean=0, eps.sigma=0) {
  # Simulate fossils assuming:
  # - Uniform deposition from theta to K
  # - Gaussian measurement error
  X = runif(n, min=theta, max=K)
  eps = rnorm(n, mean=eps.mean, sd=eps.sigma)
  W = X + eps
  return(W)
}

get_starting_est.percentile = function (alpha, theta, n, K, eps.mean, eps.sigma) {
  # Implements percentile method (Buckland, 1980; Efron, 1981)
  n.resamples = (2-alpha)/alpha
  fossil.resamples = matrix(simulate_fossils(n*n.resamples, theta, K, eps.mean, eps.sigma),
                            ncol=n.resamples)
  theta.resamples = apply(fossil.resamples, 2, FUN=min)
  theta.resamples.sorted = sort(theta.resamples)
  return(cbind(lower=theta.resamples.sorted[2],
               upper=theta.resamples.sorted[n.resamples-1]))
}

get_starting_est.analytic = function () {} # TODO

estimate_bound.rm = function (theta_q.iters, alpha.star, theta.hat, n, K, p, m, max_var=1000, max_iter, eps.mean, eps.sigma) {
  i = m
  mc.var = Inf
  while (i < max_iter && mc.var > max_var) {
    # print(c(theta_q.iters[i], theta.hat, K))
    # Generate resamples
    resamples = simulate_fossils(n, theta=theta_q.iters[i], K, eps.mean, eps.sigma)
    theta.hat.resample = estimate_theta.rm(resamples)
    # calculate step length
    c = 2*p*abs(theta_q.iters[i] - theta.hat)
    # print(paste("Step:", c))
    # Update
    if (theta.hat.resample <= theta.hat) {
      theta_q.iters[i+1] = (theta_q.iters[i] + c*(alpha.star)/ i)
    } else {
      theta_q.iters[i+1] = (theta_q.iters[i] - c*(1-alpha.star) / i)
    }
    mc.var = var.rm(alpha, c, i)
    i = i+1
    if (i >= max_iter) {
      warning("max_iter exceeded.")
    }
  }
  return(theta_q.iters)
}

var.rm = function (alpha, c, i) {
  alpha * (1-alpha) * c^2 / i
}


estimate_CI.optimal_rm = function (W, K, alpha, max_iter, eps.mean, eps.sigma, return_iters=FALSE, theta_L, theta_U, B) {
  # estimate confidence interval using Garthwaite's Robbins-Munro process (1995)
  n=length(W)
  # get point estimate of theta
  theta.hat = estimate_theta.rm(W)
  # initialize lower and upper vecs
  lower.iters = upper.iters = rep(NA, max_iter)
  # calculate m
  m = ceiling(min(50, 0.3 * (2-alpha)/alpha))
  # TODO compute optimal starting estimates 
  starting_ests = get_starting_vals(method="percentile", alpha, theta.hat, n, K, eps.mean, eps.sigma)
  lower.iters[m] = starting_ests[, "lower"]
  upper.iters[m] = starting_ests[, "upper"]
  # estimate lower
  lower.iters = estimate_bound.optimal_rm(lower.iters, alpha/2, theta.hat, n, K, m, max_var=0.2*eps.sigma, max_iter, eps.mean, eps.sigma, theta_q.true=theta_L, B=B)
  upper.iters = estimate_bound.optimal_rm(upper.iters, 1-alpha/2, theta.hat, n, K, m, max_var=0.2*eps.sigma, max_iter, eps.mean, eps.sigma, theta_q.true=theta_U, B=B)
  CI = list(CI.lower=lower.iters[length(lower.iters)], CI.upper=upper.iters[length(upper.iters)])
  if (return_iters == TRUE) {
    CI$CI.iters = cbind(lower=lower.iters, upper=upper.iters)
  }
  return(CI)
}

estimate_bound.optimal_rm = function (theta_q.iters, alpha.star, theta.hat, n, K, m, max_var=1000, max_iter, eps.mean, eps.sigma, theta_q.true, B) {
  # optimal step length c
  c = 1/calculate_g(m, K, theta_q.true, B, eps.mean, eps.sigma)
  print(c)
  i = m
  mc.var = Inf
  while (i < max_iter && mc.var > max_var) {
    # Generate resamples
    resamples = simulate_fossils(n, theta=theta_q.iters[i], K, eps.mean, eps.sigma)
    theta.hat.resample = estimate_theta.rm(resamples)
    # Update
    if (theta.hat.resample <= theta.hat) {
      theta_q.iters[i+1] = (theta_q.iters[i] + c*(alpha.star)/ i)
    } else {
      theta_q.iters[i+1] = (theta_q.iters[i] - c*(1-alpha.star) / i)
    }
    mc.var = var.rm(alpha, c, i)
    i = i+1
    if (i >= max_iter) {
      warning("max_iter exceeded.")
    }
  }
  return(theta_q.iters)
}

calculate_g = function(m, K, theta, B, eps.mean, eps.sigma) {
  e = rtnorm(B, mean=eps.mean, sd=eps.sigma, a=-Inf, b=m-theta)
  (dnorm(K-theta, mean=eps.mean, sd=eps.sigma) * ( mean((K-m)/(K-e-theta)^2) ) + (mean((K-m)/(K-e-theta)) - dnorm(m-theta, mean=eps.mean, sd=eps.sigma))*pnorm(K-theta, mean=eps.mean, sd=eps.sigma))/(dnorm(K-theta, mean=eps.mean, sd=eps.sigma))^2
}