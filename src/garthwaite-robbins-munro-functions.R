
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
  lower.iters = estimate_bound.rm("lower", lower.iters, alpha, theta.hat, n, K, p, m, max_iter, eps.mean, eps.sigma)
  upper.iters = estimate_bound.rm("upper", upper.iters, alpha, theta.hat, n, K, p, m, max_iter, eps.mean, eps.sigma)
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

estimate_bound.rm = function (bound, bound.iters, alpha, theta.hat, n, K, p, m, max_iter, eps.mean, eps.sigma) {
  # TODO - simplify this (repeating functions because of calculation of c)
  switch(bound,
         lower=estimate_bound.rm.lower(bound.iters, alpha, theta.hat, n, K, p, m, max_iter, eps.mean, eps.sigma),
         upper=estimate_bound.rm.upper(bound.iters, 1-alpha, theta.hat, n, K, p, m, max_iter, eps.mean, eps.sigma))
}

estimate_bound.rm.lower = function (bound.iters, alpha, theta.hat, n, K, p, m, max_iter, eps.mean, eps.sigma) {
  for (i in m:max_iter) {
    # Generate resamples
    resamples = simulate_fossils(n, theta=bound.iters[i], K, eps.mean, eps.sigma)
    theta.hat.resample = estimate_theta.rm(resamples)
    # calculate step length
    c = min(2*p*abs(bound.iters[i] - theta.hat), 1000)
    # Update
    if (theta.hat.resample <= theta.hat) {
      bound.iters[i+1] = (bound.iters[i] - c*(alpha/2)/ i)
    } else {
      bound.iters[i+1] = (bound.iters[i] + c*(1-alpha/2) / i)
    }
  }
  return(bound.iters)
}

estimate_bound.rm.upper = function (bound.iters, alpha, theta.hat, n, K, p, m, max_iter, eps.mean, eps.sigma) {
  for (i in m:max_iter) {
    # Generate resamples
    resamples = simulate_fossils(n, theta=bound.iters[i], K, eps.mean, eps.sigma)
    theta.hat.resample = estimate_theta.rm(resamples)
    # calculate step length
    c = min(2*p*abs(theta.hat - bound.iters[i]), 1000)
    # Update
    if (theta.hat.resample <= theta.hat) {
      bound.iters[i+1] = (bound.iters[i] + c*(1-alpha/2)/ i)
    } else {
      bound.iters[i+1] = (bound.iters[i] - c*(alpha/2) / i)
    }
  }
  return(bound.iters)
}