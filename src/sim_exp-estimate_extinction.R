# Simulation Experiments

# Load synthetic data and configuration
source("helpers-simulation-experiments.R")

load("../data/synthetic-data.RData")

attach(synthetic.data.config)
set.seed(seed)

alpha = 0.05

methods.point_estimates = c("MLE", "BA-MLE", "MINMI")
methods.conf_int = c("MINMI")

# Run trials

results = data.frame(
  id=integer(),
  error_factor=double(),
  q=double(),
  method=character(),
  estimate=double(),
  time=double()
)

for (i in 1:nrow(datasets)) {
  iter = datasets[i,]
  W = iter$W[[1]]
  print(paste("Dataset ID:", iter$id))
  
  # point estimates
  for (method in methods.point_estimates) {
    estimation = estimate_quantile(W=W, sd=iter$error_factor*fossil.sd, method=method, q=0.5, K=K, dating_error.mean=dating_error.mean)
    results = tibble::add_row(results,
                              id=iter$id,
                              error_factor=iter$error_factor,
                              q=estimation$q,
                              method=method,
                              estimate=estimation$estimate,
                              time=estimation$runtime)
  }
  
  # upper/lower estimates
  for (method in methods.conf_int) {
    # lower
    estimation = estimate_quantile(W=W, sd=iter$error_factor*fossil.sd, method=method, q=alpha/2, K=K, dating_error.mean=dating_error.mean)
    results = tibble::add_row(results,
                              id=iter$id,
                              error_factor=iter$error_factor,
                              q=estimation$q,
                              method=method,
                              estimate=estimation$estimate,
                              time=estimation$runtime)
    # upper
    estimation = estimate_quantile(W=W, sd=iter$error_factor*fossil.sd, method=method, q=1-alpha/2, K=K, dating_error.mean=dating_error.mean)
    results = tibble::add_row(results,
                              id=iter$id,
                              error_factor=iter$error_factor,
                              q=estimation$q,
                              method=method,
                              estimate=estimation$estimate,
                              time=estimation$runtime)
  }
}

View(results)
saveRDS(results, "../data/sim_exp-estimate_extinction_results.RDS")
