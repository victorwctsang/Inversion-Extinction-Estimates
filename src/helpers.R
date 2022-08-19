library(extraDistr)
library(ggplot2)
library(latex2exp)

simulateFossils = function (
  n, theta, K, eps.mean=0, eps.sigma=0
) {
  # Simulate fossils assuming:
  # - Uniform deposition from theta to K
  # - Gaussian measurement error
  X = runif(n, min=theta, max=K)
  eps = rnorm(n, mean=eps.mean, sd=eps.sigma)
  W = X + eps
  return(list(X=X, W=W, eps=eps))
}

.generateMCSamples = function (u, mean, sd, a, b) {
  # Transforms uniform mc samples to truncated normal
  n = nrow(u)
  B = ncol(u)
  mc.samples = matrix(qtnorm(p=u, mean=mean, sd=sd, a=a, b=b), ncol=B)
  return(mc.samples)
}

.calcPhiHat = function (u, mean, sd, a, b, K, m, theta) {
  # Monte Carlo integral
  e = .generateMCSamples(u, mean, sd, a, b)
  norm.constant = pnorm(b, mean=mean, sd=sd)
  phi.hat.theta = apply((K-m)/(K-e-theta), 1, mean) * norm.constant
  return(phi.hat.theta)
}

.calcA1 = function (
    theta, K, m, eps.mean, eps.sigma
) {
  F.eps.m = pnorm(m-theta, mean=eps.mean, sd=eps.sigma)
  F.eps.K = pnorm(K-theta, mean=eps.mean, sd=eps.sigma)
  return(1 - F.eps.m/F.eps.K)
}

.calcA2 = function (
  theta, K, u, m, eps.mean, eps.sigma
) {
  # Monte Carlo Integration to find phi.hat.theta
  phi.hat.theta = .calcPhiHat(u=u, mean=eps.mean, sd=eps.sigma, 
                              a=-Inf, b=m-theta, K=K, m=m, 
                              theta=theta)
  F.eps.K = pnorm(K-theta, mean=eps.mean, sd=eps.sigma)
  return(phi.hat.theta/F.eps.K)
}

getP = function (
  theta, K, u, m, eps.mean, eps.sigma
) {
  A1 = .calcA1(theta, K, m, eps.mean, eps.sigma)
  A2 = .calcA2(theta, K, u, m, eps.mean, eps.sigma)
  return(A1 + A2)
}

solve.for.theta_q.hat = function (
    theta.root, K, u, m, eps.mean, eps.sigma, n, q
) {
  # Function to solve for \hat{\theta}_q
  P = getP(theta.root, K, u, m, eps.mean, eps.sigma)
  val = prod(P) - q
  return(val)
}

getThetaQuantile = function (
  q, K, W, u=NA, eps.mean=0, eps.sigma=0, uniroot.interval=NA
) {
  n = length(W)
  m = min(W)
  theta_q.hat = NA
  if (eps.mean == 0 && eps.sigma == 0) {
    # no measurement error
    theta_q.hat = K - q^(-1/n) * (K-m)
  }
  else {
    # measurement error
    res = uniroot(
      solve.for.theta_q.hat,
      interval=uniroot.interval,
      K = K,
      n = n,
      u = u,
      q = q,
      m = m,
      eps.mean = eps.mean,
      eps.sigma = eps.sigma
    )
    theta_q.hat = res$root
  }
  return(theta_q.hat)
}

simulateConfidenceIntervals = function (
  nSim, alpha, n, theta, K, B=0, eps.mean=0, eps.sigma=0, uniroot.interval=NA
) {
  CI.L = rep(NA, nSim)
  CI.U = rep(NA, nSim)
  u = matrix(runif(n*B, min=0, max=1), ncol=B)
  
  for (i in 1:nSim) {
    sim.data = simulateFossils(n, theta, K, eps.mean, eps.sigma)
    CI.L[i] = getThetaQuantile(q=alpha/2,
                               K=K,
                               W=sim.data$W,
                               u=u,
                               eps.mean=eps.mean,
                               eps.sigma=eps.sigma,
                               uniroot.interval=uniroot.interval)
    CI.U[i] = getThetaQuantile(q=1-alpha/2,
                               K=K,
                               W=sim.data$W,
                               u=u,
                               eps.mean=eps.mean,
                               eps.sigma=eps.sigma,
                               uniroot.interval=uniroot.interval)
  }
  CI.df = data.frame(i=1:nSim, L=CI.L, U=CI.U)
  CI.df$mean = rowMeans(CI.df[,2:3])
  CI.df$containsTheta = (CI.df$L < theta) & (theta < CI.df$U)
  return(CI.df)
}

plot.simCIs = function(df, theta) {
  ggplot(data=df, aes(x=i, y=mean)) +
    geom_errorbar(aes(ymax = U, ymin = L, colour=containsTheta)) +
    geom_hline(yintercept=theta) +
    labs(title=paste("Confidence interval;", TeX("\\theta"), "=", theta),
         y=TeX("\\theta"))
}

