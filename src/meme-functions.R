library(extraDistr)

estimate_CI.mc = function (
    alpha, n, K, W, u, eps.sigma, eps.mean, uniroot.interval
  ) {
  CI.lower = getThetaQuantile(q=alpha/2, K, W, u, eps.mean, eps.sigma, uniroot.interval)
  CI.upper = CI.lower
  if (alpha != 0.5) {
    CI.upper = getThetaQuantile(q=1-alpha/2, K, W, u, eps.mean, eps.sigma, uniroot.interval)
  }
  CI = list(CI.lower = CI.lower, CI.upper = CI.upper)
  return(CI)
}

estimate_quantile.mc = function (
  q, K, W, u=NA, eps.mean=NA, eps.sigma=NA, uniroot.interval=NA
) {
  n = length(W)
  m = min(W)
  theta_q.hat = NA
  if (all(is.na(c(u, eps.mean, eps.sigma, uniroot.interval)))) {
    # No Measurement Error Case
    theta_q.hat = K - q^(-1/n) * (K-m)
  } else if (any(is.na(c(u, eps.mean, eps.sigma, uniroot.interval)))) {
    # Edge case
    stop("Not all parameters for the measurement error estimation have been provided. Check your inputs and try again.")
  } else {
    # Measurement Error case
    res = uniroot(function(theta) prod(getP(theta, K, u, m, eps.mean, eps.sigma)) - q)
    theta_q.hat = res$root
  }
  return(theta_q.hat)
}

calc_P = function (theta, K, u, m, eps.mean, eps.sigma) {
  # P(W >= m)
  F.eps.m = pnorm(m-theta, mean=eps.mean, sd=eps.sigma)
  F.eps.K = pnorm(K-theta, mean=eps.mean, sd=eps.sigma)
  # Monte Carlo Integration to find phi.hat
  phi.hat = estimate_phi(u=u, mean=eps.mean, sd=eps.sigma, a=-Inf, b=m-theta, K=K, m=m, theta=theta)
  F.eps.K = pnorm(K-theta, mean=eps.mean, sd=eps.sigma)
  P = 1 - F.eps.m/F.eps.K + phi.hat/F.eps.K
  return(P)
}

estimate_phi = function (u, mean, sd, a, b, K, m, theta) {
  # Monte Carlo integral
  e = uniform_to_tnorm(u, mean, sd, a, b)
  norm.constant = pnorm(b, mean=mean, sd=sd)
  phi.hat = apply((K-m)/(K-e-theta), 1, mean) * norm.constant
  return(phi.hat)
}

uniform_to_tnorm = function (u, mean, sd, a, b) {
  # Transforms uniform mc samples to truncated normal
  n = nrow(u)
  B = ncol(u)
  mc.samples = matrix(qtnorm(p=u, mean=mean, sd=sd, a=a, b=b), ncol=B)
  return(mc.samples)
}