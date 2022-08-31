

calc_step_length_k = function (alpha) {
  # K := Step length proportionality constant
  # Garthwaite 1995 paper
  # Normal distribution heuristic
  2 / (qnorm(1-alpha) * (2*pi)^(-0.5) * exp(- (qnorm(1-alpha)^2)/2))
}

calc_unbiased_estimator = function (X, K, n) {
  min(X) - (K - min(X))/(n-1)
}

generate_starting_estimates.percentile = function (alpha, theta, n, K) {
  # Implements percentile method (Buckland, 1980; Efron, 1981)
  n.resamples = (2-alpha)/alpha
  fossil.resamples = matrix(runif(n * n.resamples, min=theta, max=K),
                            ncol=n.resamples)
  theta.resamples = apply(fossil.resamples,
                          2,
                          FUN=min)#function (x) calc_unbiased_estimator(x, b, n))
  theta.resamples.sorted = sort(theta.resamples)
  return(list(L=theta.resamples.sorted[2],
              U=theta.resamples.sorted[n.resamples-1],
              m=ceiling(min(50, 0.3 * n.resamples))))
}

rm_process = function (
    bound_vec, theta.hat, n, alpha, K, m, max_iter, k
) {
  for (i in m:max_iter) {
    # Generate resamples
    resamples = runif(n=n, min=bound_vec[i], max=K)
    
    # calculate step length
    c = min(2*k*abs(bound_vec[i] - theta.hat), 1000) # TODO
    
    # Update upper bound
    if (min(resamples) <= theta.hat) {
      bound_vec[i+1] = (bound_vec[i] + c*(1-alpha)/ i)
    } 
    else {
      bound_vec[i+1] = (bound_vec[i] - c*alpha / i)
    }
  }
  return(bound_vec[m:max_iter])
}

getConfidenceInterval.rm_process = function (theta.hat, n, alpha, K, max_iter, return_iters=FALSE) {
  # Step length proportionality constant
  k = calc_step_length_k(alpha)
  # Starting values
  starting_est = generate_starting_estimates.percentile(alpha, theta.hat, n, K)
  
  bound_lower = bound_upper = rep(NA, max_iter)
  m = starting_est$m
  bound_lower[m] = starting_est$L
  bound_upper[m] = starting_est$U
  
  lower = rm_process(bound_vec = bound_lower,
                     theta.hat = theta.hat,
                     n = n,
                     alpha = (1 - alpha),
                     K = K,
                     m = m,
                     max_iter = max_iter,
                     k = k)
  upper = rm_process(bound_vec = bound_upper,
                     theta.hat = theta.hat,
                     n = n,
                     alpha = alpha,
                     K = K,
                     m = m,
                     max_iter = max_iter,
                     k = k)
  CI.iters = cbind(lower, upper)
  return_vals = list(L=lower[length(lower)],
                     U=upper[length(upper)])
  return_vals$mean = mean(return_vals$L, return_vals$U)
  return_vals$width = return_vals$U - return_vals$L
  if (return_iters==TRUE) {
    return_vals$CI.iters=CI.iters
  }
  return(return_vals)
}