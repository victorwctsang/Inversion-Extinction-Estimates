library(extraDistr)
# Generic function for MINMI point estimates and confidence intervals
estimate_extinction.minmi = function (W,
                                      sd,
                                      K,
                                      level = NULL,
                                      B = NULL,
                                      time = F,
                                      .B_init = 500,
                                      .max_var = NULL,
                                      return_Bs = F) {
  result = list(point = NULL)
  n = length(W)
  m = min(W)
  dating.sd = mean(sd)
  
  start_time = Sys.time()
  
  # Set the maximum variance
  max_var = 0.2 * (dating.sd) ^ 2
  if (!is.null(.max_var)) {
    max_var = .max_var
  }
  
  # Calculate number of monte carlo samples to use
  B.point = B.lower = B.upper = NULL
  if (!is.null(B)) {
    B.point = B.lower = B.upper = B
  } else {
    u.init = runif(.B_init, 0, 1)
    B.point = find_optimal_B(
      max_var = max_var,
      q = 0.5,
      K = K,
      m = m,
      n = n,
      u = u.init,
      eps.mean = 0,
      eps.sigma = dating.sd
    )
    if (!is.null(level)) {
      alpha = (1 - level) / 2
      B.lower = find_optimal_B(
        max_var = max_var,
        q = alpha,
        K = K,
        m = m,
        n = n,
        u = u.init,
        eps.mean = 0,
        eps.sigma = dating.sd
      )
      B.upper = find_optimal_B(
        max_var = max_var,
        q = 1 - alpha,
        K = K,
        m = m,
        n = n,
        u = u.init,
        eps.mean = 0,
        eps.sigma = dating.sd
      )
    }
  }
  
  # Generate Monte Carlo Samples
  B.max = max(B.point, B.lower, B.upper)
  mc.samples = runif(B.max, min = 0, max = 1)
  
  # Calculate estimates
  if (!is.null(level)) {
    alpha = (1 - level) / 2
    result$lower = estimate_quantile.minmi(
      q = alpha,
      K = K,
      W = W,
      u = mc.samples[1:B.lower],
      eps.mean = 0,
      eps.sigma = dating.sd
    )
    result$upper = estimate_quantile.minmi(
      q = 1 - alpha,
      K = K,
      W = W,
      u = mc.samples[1:B.upper],
      eps.mean = 0,
      eps.sigma = dating.sd
    )
  }
  result$point = estimate_quantile.minmi(
    q = 0.5,
    K = K,
    W = W,
    u = mc.samples[1:B.point],
    eps.mean = 0,
    eps.sigma = dating.sd
  )
  
  if (time) {
    result$time = Sys.time() - start_time
    print(result$time)
  }
  
  if (return_Bs == T) {
    result$B.point = B.point
    result$B.lower = B.lower
    result$B.upper = B.upper
  }
  return(result)
}

# Function for getting a MINMI estimate of some quantile q
estimate_quantile.minmi = function (q,
                                    K,
                                    W,
                                    u = NA,
                                    eps.mean = 0,
                                    eps.sigma = 0) {
  n = length(W)
  m = min(W)
  # No Measurement Error case
  theta_q.hat = K - q ^ (-1 / n) * (K - m)
  if (eps.mean != 0 || eps.sigma != 0) {
    # Measurement Error case
    newton.res = pracma::newtonRaphson(
      fun = function(theta)
        estimating_eqn(theta, q, K, u, m, n, eps.mean, eps.sigma),
      x0 = theta_q.hat
    )
    theta_q.hat = newton.res$root
  }
  return(theta_q.hat)
}

# Function for estimating the number of B to use
find_optimal_B = function (max_var, q, K, m, n, u, eps.mean, eps.sigma) {
  # Initial estimate of theta_q using no measurement error case
  theta_q.hat.init = K - q ^ (-1 / n) * (K - m)
  
  # pdfs and CDF evaluations (for convenience)
  f_eps.K = dnorm(K - theta_q.hat.init, eps.mean, eps.sigma)
  F_eps.K = pnorm(K - theta_q.hat.init, eps.mean, eps.sigma)
  # Estimate sigma^2_var
  B = length(u)
  e = uniform_to_tnorm(u,
                       eps.mean,
                       eps.sigma,
                       a = -Inf,
                       b = m - theta_q.hat.init)
  sample_var.psi_hat = var((m - e - theta_q.hat.init) / (K - e - theta_q.hat.init))
  
  # Estimate \hat\psi and \hat\psi prime
  psi_hat = mean((m - e - theta_q.hat.init) / (K - e - theta_q.hat.init))
  psi_hat.prime = -mean((K - m) / (K - e - theta_q.hat.init) ^ 2)
  
  optimal_B = 1 / max_var * sample_var.psi_hat * (f_eps.K / F_eps.K * psi_hat + psi_hat.prime) ^
    (-2)
  
  return(max(ceiling(optimal_B), 100))
}

# Estimating equation for Newton-Raphson's
estimating_eqn = function (theta, q, K, u, m, n, eps.mean, eps.sigma) {
  F.eps.m = pnorm(m - theta, mean = eps.mean, sd = eps.sigma)
  F.eps.K = pnorm(K - theta, mean = eps.mean, sd = eps.sigma)
  psi.hat = estimate_psi(
    u = u,
    mean = eps.mean,
    sd = eps.sigma,
    a = -Inf,
    b = m - theta,
    K = K,
    m = m,
    theta = theta
  )
  return(1 - F.eps.m / F.eps.K * psi.hat - q ^ (1 / n))
}

estimating_eqn_deriv = function (theta, K, u, m, eps.mean, eps.sigma) {
  F.eps.m = pnorm(m - theta, mean = eps.mean, sd = eps.sigma)
  f.eps.K = dnorm(K - theta, mean = eps.mean, sd = eps.sigma)
  F.eps.K = pnorm(K - theta, mean = eps.mean, sd = eps.sigma)
  psi.hat = estimate_psi(
    u = u,
    mean = eps.mean,
    sd = eps.sigma,
    a = -Inf,
    b = m - theta,
    K = K,
    m = m,
    theta = theta
  )
  psi.hat.prime = 
    mean((m - K) / (m - uniform_to_tnorm(u, eps.mean, eps.sigma,-Inf, b = m - theta) - theta) ^ 2)
  
  return(-F.eps.m / F.eps.K * (f.eps.K / F.eps.K * psi.hat + psi.hat.prime))
}

# Monte Carlo approximation of psi
estimate_psi = function (u, mean, sd, a, b, K, m, theta) {
  e = uniform_to_tnorm(u, mean, sd, a, b)
  psi.hat = mean((m - e - theta) / (K - e - theta))
  return(psi.hat)
}
# Helper function for getting Monte Carlo samples
# using a density transformation
uniform_to_tnorm = function (u, mean, sd, a, b) {
  mc.samples = qtnorm(
    p = u,
    mean = mean,
    sd = sd,
    a = a,
    b = b)
  return(mc.samples)
}