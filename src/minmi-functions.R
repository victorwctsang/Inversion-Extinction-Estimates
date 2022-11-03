library(extraDistr)

estimate_extinction.minmi = function (W, sd, K, level = NULL, B = NULL, time=F, .B_init=500) {
  result = list(point=NULL)
  n = length(W)
  m = min(W)
  dating.sd = mean(sd)
  
  start_time = Sys.time()
  # Calculate number of monte carlo samples to use
  B.extinction = B.lower = B.upper = NULL
  if (!is.null(B)) {
    B.extinction = B.lower = B.upper = B
  } else {
    u.init = runif(.B_init, 0, 1)
    B.extinction = find_optimal_B(max_var = (0.2*dating.sd)^2,
                                  q=0.5,
                                  K=K,
                                  m=m,
                                  u=u.init,
                                  eps.mean=0,
                                  eps.sigma=dating.sd)
    if (!is.null(level)) {
      alpha = (1 - level)/2
      B.lower = find_optimal_B(max_var = (0.2*dating.sd)^2,
                               q=alpha,
                               K=K,
                               m=m,
                               u=u.init,
                               eps.mean=0,
                               eps.sigma=dating.sd)
      B.upper = find_optimal_B(max_var = (0.2*dating.sd)^2,
                               q=1-alpha,
                               K=K,
                               m=m,
                               u=u.init,
                               eps.mean=0,
                               eps.sigma=dating.sd)
    }
  }
  
  # Generate Monte Carlo Samples
  B.max = max(B.extinction, B.lower, B.upper)
  mc.samples = runif(B.max, min=0, max=1)

  # Perform estimates  
  if (!is.null(level)) {
    alpha = (1 - level)/2
    result$lower = estimate_quantile.minmi(q=alpha, K=K, W=W, u=mc.samples[1:B.lower],
                                           eps.mean=0, eps.sigma=dating.sd)
    result$upper = estimate_quantile.minmi(q=1-alpha, K=K, W=W, u=mc.samples[1:B.upper],
                                           eps.mean=0, eps.sigma=dating.sd)
  }
  result$point = estimate_quantile.minmi(q=0.5, K=K, W=W, u=mc.samples[1:B.extinction],
                                              eps.mean=0, eps.sigma=dating.sd)
  
  if (time) {
    result$time = Sys.time() - start_time
    print(result$time)
  }
  return(result)
}

estimate_CI.minmi = function (
    alpha, n, K, W, u, eps.mean, eps.sigma
  ) {
  CI.lower = estimate_quantile.mc(q=alpha/2, K, W, u, eps.mean, eps.sigma)
  CI.upper = CI.lower
  if (alpha != 0.5) {
    CI.upper = estimate_quantile.mc(q=1-alpha/2, K, W, u, eps.mean, eps.sigma)
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
                                       x0=theta_q.hat)
    theta_q.hat = newton.res$root
  }
  return(theta_q.hat)
}

estimating_eqn = function (theta, q, K, u, m, eps.mean, eps.sigma) {
  F.eps.m = pnorm(m-theta, mean=eps.mean, sd=eps.sigma)
  F.eps.K = pnorm(K-theta, mean=eps.mean, sd=eps.sigma)
  psi.hat = estimate_psi(u=u, mean=eps.mean, sd=eps.sigma, a=-Inf, b=m-theta, K=K, m=m, theta=theta)
  return(1 - F.eps.m/F.eps.K * psi.hat - q^(1/n))
}

estimating_eqn_deriv = function (theta, K, u, m, eps.mean, eps.sigma) {
  F.eps.m = pnorm(m-theta, mean=eps.mean, sd=eps.sigma)
  f.eps.K = dnorm(K-theta, mean=eps.mean, sd=eps.sigma)
  F.eps.K = pnorm(K-theta, mean=eps.mean, sd=eps.sigma)
  psi.hat = estimate_psi(u=u, mean=eps.mean, sd=eps.sigma, a=-Inf, b=m-theta, K=K, m=m, theta=theta)
  psi.hat.prime = mean( (m-K)/(m-uniform_to_tnorm(u, eps.mean, eps.sigma, -Inf, b=m-theta)-theta)^2 )
  
  return(-F.eps.m/F.eps.K * (f.eps.K/F.eps.K * psi.hat + psi.hat.prime))
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

find_optimal_B = function (max_var, q, K, m, u, eps.mean, eps.sigma) {
  # Initial estimate of theta_q using no measurement error case
  theta_q.hat.init = K - q^(-1/n)*(K-m) 
  
  # pdfs and CDF evaluations (for convenience)
  f_eps.K = dnorm(K-theta_q.hat.init, eps.mean, eps.sigma)
  F_eps.K = pnorm(K-theta_q.hat.init, eps.mean, eps.sigma)
  
  # Estimate sigma^2_var
  B = length(u)
  e = uniform_to_tnorm(u, eps.mean, eps.sigma, a=-Inf, b=m-theta_q.hat.init)
  sample_var.psi_hat = var((m-e-theta_q.hat.init)/(K-e-theta_q.hat.init))
  
  # Estimate \hat\psi and \hat\psi prime
  psi_hat = mean((m-e-theta_q.hat.init)/(K-e-theta_q.hat.init))
  psi_hat.prime = - mean((K-m)/(K-e-theta_q.hat.init)^2)
  
  optimal_B = 1/max_var * sample_var.psi_hat * ( f_eps.K/F_eps.K * psi_hat + psi_hat.prime )^(-2)
  return(ceiling(optimal_B))
}