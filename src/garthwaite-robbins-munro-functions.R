
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
  lower.iters = estimate_bound.rm(lower.iters, alpha/2, theta.hat, n, K, p, m, max_var=(0.2*eps.sigma)^2, max_iter, eps.mean, eps.sigma)
  upper.iters = estimate_bound.rm(upper.iters, 1-alpha/2, theta.hat, n, K, p, m, max_var=(0.2*eps.sigma)^2, max_iter, eps.mean, eps.sigma)
  lower.iters = na.omit(lower.iters)
  upper.iters = na.omit(upper.iters)

  CI = list(lower=lower.iters[length(lower.iters)],
            point=theta.hat,
            upper=upper.iters[length(upper.iters)])
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
               point=median(theta.resamples.sorted),
               upper=theta.resamples.sorted[n.resamples-1]))
}

get_starting_est.analytic = function () {
  # TODO 
}

estimate_bound.rm = function (theta.iters, q, theta.hat, n, K, p, m, max_var=1000, max_iter, eps.mean, eps.sigma) {
  # Transformation to force theta < K
  eta = function(theta) -log(K-theta)
  eta_q.iters = eta(theta.iters)
  i = m
  mc.var = Inf
  while (i < max_iter) {
    # Generate resamples
    theta_q.hat = K - exp(-eta_q.iters[i])
    resamples = simulate_fossils(n, theta=theta_q.hat, K, eps.mean, eps.sigma)
    theta.hat.resample = estimate_theta.rm(resamples)
    # calculate step length
    c = 2*p*abs(eta(theta_q.hat) - eta(theta.hat))
    
    # Update
    if (theta.hat.resample <= theta.hat) {
      eta_q.iters[i+1] = (eta_q.iters[i] + c*q/ i)
    } else {
      eta_q.iters[i+1] = (eta_q.iters[i] - c*(1-q) / i)
    }
    i = i+1
  }
  theta_q.iters = K - exp(-eta_q.iters)
  return(theta_q.iters)
}

var.rm = function (alpha, c, i) {
  # TODO - pretty sure this is wrong
  alpha * (1-alpha) * c^2 / i
}