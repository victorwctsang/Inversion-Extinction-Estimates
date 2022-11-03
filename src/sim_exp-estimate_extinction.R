# Simulation Experiments

# Load synthetic data and configuration
source("helpers-simulation-experiments.R")

load("../data/synthetic-data.RData")

attach(synthetic.data.config)
set.seed(seed)

alpha = 0.05

methods.point_estimates = c("MLE", "BA-MLE", "STRAUSS")
methods.conf_int = c("GRIWM", "MINMI", "SI-UGM", "GB-RM") #"GRIWM-corrected", 

# Run trials
# results = readRDS("../data/sim_exp-estimate_extinction_results.RDS")
results = data.frame(
  id=integer(),
  error_factor=double(),
  method=factor(),
  lower=double(),
  point=double(),
  upper=double(),
  point_runtime=double(),
  conf_int_runtime=double()
)

start_time = Sys.time()
for (i in 1:nrow(datasets)) {
  iter = datasets[i,]
  W = as.numeric(iter$W[[1]])
  sd = as.numeric(iter$error_factor*fossil.sd)
  print(paste("Dataset ID:", iter$id))
  
  print("- Point Estimates")
  for (method in methods.point_estimates) {
    print(paste("--", method))
    estimation = estimate_extinction(W=W, sd=sd, method=method, K=K, dating_error.mean=dating_error.mean)
    results = tibble::add_row(results, id=iter$id, error_factor=iter$error_factor, method=method, point=estimation$point, point_runtime=estimation$point_runtime, conf_int_runtime=estimation$conf_int_runtime)
  }
  
  print("- Confidence Intervals")
  for (method in methods.conf_int) {
    print(paste("--", method))
    estimation = estimate_conf_int(W=W, sd=sd, method=method, alpha=alpha, K=K, dating_error.mean=dating_error.mean)
    results = tibble::add_row(results, id=iter$id, error_factor=iter$error_factor, method=method, lower=estimation$lower, point=estimation$point, upper=estimation$upper, point_runtime=estimation$point_runtime, conf_int_runtime=estimation$conf_int_runtime)
  }
}
print(Sys.time() - start_time)

View(results)
saveRDS(results, "../data/sim_exp-estimate_extinction_results.RDS")
