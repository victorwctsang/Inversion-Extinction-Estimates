library(extraDistr)

estimate_extinction.minmi = function (W, sd, K, level = NA) {
  result = list(extinction=NA)
  n = length(W)
  m = min(W)
  dating.sd = mean(sd)
  uniroot.interval = c(5000, 19500) #c(m-(K-m)/2, m+(K-m)/2)
  
  B = 100
  u = matrix(runif(n*B, min=0, max=1), ncol=B)
  
  result$extinction = estimate_quantile.minmi(q=0.5, K=K, W=W, u=u,
                                              eps.mean=0, eps.sigma=dating.sd,
                                              uniroot.interval=uniroot.interval)
  if (!is.na(level)) {
    alpha = (1 - level)/2
    result$lower = estimate_quantile.minmi(q=alpha, K=K, W=W, u=u,
                                           eps.mean=0, eps.sigma=dating.sd,
                                           uniroot.interval=uniroot.interval)
    result$upper = estimate_quantile.minmi(q=1-alpha, K=K, W=W, u=u,
                                           eps.mean=0, eps.sigma=dating.sd,
                                           uniroot.interval=uniroot.interval)
  }
  return(result)
}

estimate_CI.minmi = function (
    alpha, n, K, W, u, eps.mean, eps.sigma, uniroot.interval
  ) {
  CI.lower = estimate_quantile.mc(q=alpha/2, K, W, u, eps.mean, eps.sigma, uniroot.interval)
  CI.upper = CI.lower
  if (alpha != 0.5) {
    CI.upper = estimate_quantile.mc(q=1-alpha/2, K, W, u, eps.mean, eps.sigma, uniroot.interval)
  }
  CI = list(CI.lower = CI.lower, CI.upper = CI.upper)
  return(CI)
}

estimate_quantile.minmi = function (
  q, K, W, u=NA, eps.mean=0, eps.sigma=0, uniroot.interval=NA
) {
  n = length(W)
  m = min(W)
  theta_q.hat = NA
  if (eps.mean==0 & eps.sigma==0) {
    # No Measurement Error Case
    theta_q.hat = K - q^(-1/n) * (K-m)
  } else  {
    # Measurement Error case
    theta_q.hat = tryCatch(
      {
        uniroot(function(theta) prod(calc_P(theta, K, u, m, eps.mean, eps.sigma)) - q,
                interval=uniroot.interval)$root
      },
      error=function(cond) {
        message("Error in uniroot")
        message(cond)
        return(NA)
      }
    )
  }
  return(theta_q.hat)
}

calc_P = function (theta, K, u, m, eps.mean, eps.sigma) {
  # P(W >= m)
  F.eps.m = pnorm(m-theta, mean=eps.mean, sd=eps.sigma)
  F.eps.K = pnorm(K-theta, mean=eps.mean, sd=eps.sigma)
  psi.hat = estimate_psi(u=u, mean=eps.mean, sd=eps.sigma, a=-Inf, b=m-theta, K=K, m=m, theta=theta)
  P = 1 - F.eps.m/F.eps.K * psi.hat
  return(P)
}

estimate_psi = function (u, mean, sd, a, b, K, m, theta) {
  # Monte Carlo integral
  e = uniform_to_tnorm(u, mean, sd, a, b)
  psi.hat = apply((m-e-theta)/(K-e-theta), 1, mean)
  return(psi.hat)
}

uniform_to_tnorm = function (u, mean, sd, a, b) {
  # Transforms uniform mc samples to truncated normal
  n = nrow(u)
  B = ncol(u)
  mc.samples = matrix(qtnorm(p=u, mean=mean, sd=sd, a=a, b=b), ncol=B)
  return(mc.samples)
}