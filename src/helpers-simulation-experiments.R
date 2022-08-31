
library(readxl)

getFossilData = function (path, sheet, range, col_names, col_types) {
  read_excel(path, sheet, range, col_names, col_types)
}

boundFossilData = function (fossil.df, K) {
  return(fossil.df[fossil.df$age < K, ])
}

getStdDevFromFossilData = function (path, K, ...) {
  df = getFossilData(path, ...)
  df = boundFossilData(df, K)
  return(mean(df$sd))
}

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