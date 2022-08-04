library(extraDistr)
library(ggplot2)

simulateFossils = function (
  n, theta, K, eps.mean=0, eps.sigma=0
) {
  X = runif(n, min=theta, max=K)
  eps = rnorm(n, mean=eps.mean, sd=eps.sigma)
  W = X + eps
  return(list(X=X, W=W, eps=eps))
}

.generateMCSamples = function (n, mean, sd, a, b) {
  rtnorm(n, mean=mean, sd=sd, a=a, b=b)
}

.F.eps = function (q, mean, sd) {
  pnorm(q, mean, sd)
}

.a.theta = function (theta, K, m, mean, sd, phi) {
  F.eps.m = .F.eps(m-theta, mean, sd)
  F.eps.K = .F.eps(K-theta, mean, sd)
  
  a = (K-theta)/(K-m) * (1 - F.eps.m/F.eps.K) +
    F.eps.m/F.eps.K +
    F.eps.m/F.eps.K * phi
  return(a)
}

solve.for.theta_q.hat = function (
  theta.root, K, n, B, q, m, eps.mean, eps.sigma
) {
  # Monte Carlo Integration to find phi.hat.theta
  mc.samples = matrix(.generateMCSamples(n*B, eps.mean, eps.sigma, -Inf, m-theta.root), ncol=B)
  phi.hat.theta = apply(mc.samples/(K-theta.root-mc.samples), 1, mean)
  a.theta = .a.theta(theta.root, K, m, eps.mean, eps.sigma, phi.hat.theta)
  # Return the equation we're trying to solve
  return(theta.root - K + q^(-1/n) * (K-m) * exp(mean(log(a.theta))))
}

getP = function (
  theta, K, B, m, eps.mean, eps.sigma
) {
  # Monte Carlo Integration to find phi.hat.theta
  e = .generateMCSamples(B, eps.mean, eps.sigma, -Inf, m-theta)
  phi.hat.theta = mean(e/(K-theta-e))
  a.theta = .a.theta(theta, K, m, eps.mean, eps.sigma, phi.hat.theta)
  P = (K-m)/(K-theta) * a.theta
  return(P)
}

getThetaQuantile = function (
  q, K, W, B=0, eps.mean=0, eps.sigma=0, uniroot.interval=NA
) {
  n = length(W)
  m = min(W)
  theta_q.hat = NA
  if (eps.mean == 0 && eps.sigma == 0 && B == 0) {
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
      B = B,
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

  for (i in 1:nSim) {
    sim.data = simulateFossils(n, theta, K, eps.mean, eps.sigma)
    CI.L[i] = getThetaQuantile(q=alpha/2,
                               K=K,
                               W=sim.data$W,
                               B=B,
                               eps.mean=eps.mean,
                               eps.sigma=eps.sigma,
                               uniroot.interval=uniroot.interval)
    CI.U[i] = getThetaQuantile(q=1-alpha/2,
                               K=K,
                               W=sim.data$W,
                               B=B,
                               eps.mean=eps.mean,
                               eps.sigma=eps.sigma,
                               uniroot.interval=uniroot.interval)
  }
  CI.df = data.frame(i=1:nSim, L=L, U=U)
  CI.df$mean = rowMeans(CI.df[,2:3])
  CI.df$containsTheta = (CI.df$L < theta) & (theta < CI.df$U)
  return(CI.df)
}

plot.simCIs = function(df, theta) {
  ggplot(data=df, aes(x=i, y=mean)) +
    geom_errorbar(aes(ymax = U, ymin = L, colour=containsTheta)) +
    geom_hline(yintercept=theta) +
    ggtitle(paste("Confidence interval;", "theta", "=", theta)) +
    ylab("Theta")
}

