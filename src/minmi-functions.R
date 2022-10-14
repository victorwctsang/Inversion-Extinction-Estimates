library(extraDistr)

estimate_extinction.minmi = function (W, sd, K, level = NULL, B = NULL, speed=F) {
  result = list(extinction=NULL)
  n = length(W)
  m = min(W)
  dating.sd = mean(sd)
  uniroot.interval = c(m-(K-m)/2, m+(K-m)/2)
  
  start_time = Sys.time()
  # Calculate number of monte carlo samples to use
  B.extinction = B.lower = B.upper = NULL
  if (!is.null(B)) {
    B.extinction = B.lower = B.upper = B
  } else {
    B.extinction = find_optimal_B(max_var = (0.2*dating.sd)^2,
                                  q=0.5,
                                  K=K,
                                  m=m,
                                  eps.mean=0,
                                  eps.sigma=dating.sd,
                                  B=500)
    if (!is.null(level)) {
      alpha = (1 - level)/2
      B.lower = find_optimal_B(max_var = (0.2*dating.sd)^2,
                               q=alpha,
                               K=K,
                               m=m,
                               eps.mean=0,
                               eps.sigma=dating.sd,
                               B=500)
      B.upper = find_optimal_B(max_var = (0.2*dating.sd)^2,
                               q=1-alpha,
                               K=K,
                               m=m,
                               eps.mean=0,
                               eps.sigma=dating.sd,
                               B=500)
    }
  }
  
  # Generate Monte Carlo Samples
  B.max = max(B.extinction, B.lower, B.upper)
  u = runif(B.max, min=0, max=1)

  # Perform estimates  
  if (!is.null(level)) {
    alpha = (1 - level)/2
    result$lower = estimate_quantile.minmi(q=alpha, K=K, W=W, u=u[1:B.lower],
                                           eps.mean=0, eps.sigma=dating.sd,
                                           uniroot.interval=uniroot.interval)
    result$upper = estimate_quantile.minmi(q=1-alpha, K=K, W=W, u=u[, 1:B.upper],
                                           eps.mean=0, eps.sigma=dating.sd,
                                           uniroot.interval=uniroot.interval)
  }
  result$extinction = estimate_quantile.minmi(q=0.5, K=K, W=W, u=u[1:B.extinction],
                                              eps.mean=0, eps.sigma=dating.sd,
                                              uniroot.interval=uniroot.interval)
  
  if (speed) {print(Sys.time() - start_time)}
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

estimate_quantile.minmi = function (q, K, W, u=NA, eps.mean=0, eps.sigma=0) {
  n = length(W)
  m = min(W)
  # No Measurement Error case
  theta_q.hat = K - q^(-1/n) * (K-m)
  if (eps.mean!=0 || eps.sigma!=0) {
    # Measurement Error case
    newton.res = pracma::newtonRaphson(fun=function(theta) estimating_eqn(theta, q, K, u, m, eps.mean, eps.sigma),
                                       x0=theta_q.hat,
                                       dfun=function(theta) estimating_eqn.deriv(theta, K, u, m, eps.mean, eps.sigma),
                                       maxiter = 200)
    theta_q.hat = newton.res$root
  }
  return(theta_q.hat)
}

estimating_eqn = function (theta, q, K, u, m, eps.mean, eps.sigma) {
  F.eps.m = pnorm(m-theta, mean=eps.mean, sd=eps.sigma)
  F.eps.K = pnorm(K-theta, mean=eps.mean, sd=eps.sigma)
  psi.hat = estimate_psi(u=u, mean=eps.mean, sd=eps.sigma, a=-Inf, b=m-theta, K=K, m=m, theta=theta)
  1 - F.eps.m/F.eps.K * psi.hat - q^(1/n)
}

estimating_eqn.deriv = function (theta, K, u, m, eps.mean, eps.sigma) {
  # pdfs and CDF evaluations (for convenience)
  f_eps.K = dnorm(K-theta, eps.mean, eps.sigma)
  F_eps.K = pnorm(K-theta, eps.mean, eps.sigma)
  f_eps.m = dnorm(m-theta, eps.mean, eps.sigma)
  F_eps.m = pnorm(m-theta, eps.mean, eps.sigma)
  
  e = uniform_to_tnorm(u, eps.mean, eps.sigma, a=-Inf, b=m-theta)
  psi_hat = mean((m-e-theta)/(K-e-theta))
  psi_hat.prime = - mean((K-m)/(K-e-theta)^2)

  - psi_hat * (F_eps.m*f_eps.K - f_eps.m*F_eps.K)/F_eps.K^2 - psi_hat.prime
}

estimate_psi = function (u, mean, sd, a, b, K, m, theta) {
  # Monte Carlo integral
  e = uniform_to_tnorm(u, mean, sd, a, b)
  psi.hat = mean((m-e-theta)/(K-e-theta))
  return(psi.hat)
}

uniform_to_tnorm = function (u, mean, sd, a, b) {
  # Transforms uniform mc samples to truncated normal
  mc.samples = qtnorm(p=u, mean=mean, sd=sd, a=a, b=b)
  return(mc.samples)
}

find_optimal_B = function (max_var, q, K, m, eps.mean, eps.sigma, B=500) {
  # Initial estimate of theta_q using no measurement error case
  theta_q.hat.init = K - q^(-1/n)*(K-m) 
  
  # Monte carlo samples
  e = rtnorm(B, mean=eps.mean, sd=eps.sigma, a=-Inf, b=m-theta_q.hat.init)
  
  f_eps.K = dnorm(K-theta_q.hat.init, eps.mean, eps.sigma)
  F_eps.K = pnorm(K-theta_q.hat.init, eps.mean, eps.sigma)
  f_eps.m = dnorm(m-theta_q.hat.init, eps.mean, eps.sigma)
  F_eps.m = pnorm(m-theta_q.hat.init, eps.mean, eps.sigma)
  
  sigma2.psi_hat = var((m-e-theta_q.hat.init)/(K-e-theta_q.hat.init))
  psi_hat = mean((m-e-theta_q.hat.init)/(K-e-theta_q.hat.init)) * F_eps.m
  psi_hat.prime = mean((m-K)/(K-e-theta_q.hat.init)^2) * F_eps.m
  u.prime = psi_hat * ( (F_eps.m * f_eps.K)/F_eps.K^2 - f_eps.m/F_eps.K) + F_eps.m/F_eps.K * psi_hat.prime
  
  ceiling(sigma2.psi_hat/max_var * (F_eps.m/F_eps.K)^2 * (u.prime)^(-2))
}