# Simulation Experiments

# Load synthetic data and configuration
source("../src/helpers-simulation-experiments.R")
source("../src/minmi-functions.R")

load("../data/synthetic-data.RData")

attach(synthetic.data.config)
set.seed(seed)

alpha = 0.05

methods.point_estimates = c()#"STRAUSS", "MLE", "BA-MLE")
methods.conf_int = c("MINMI")#"GRIWM", "GRIWM-corrected", "GRIWM", "SI-RM", "SI-UGM", "SI-RM-corrected")

# Run trials
results = readRDS("../data/sim_exp-estimate_extinction_results.RDS")
results = results[results$method != "MINMI", ]
# results = data.frame(
#   id=integer(),
#   error_factor=double(),
#   method=factor(),
#   lower=double(),
#   point=double(),
#   upper=double(),
#   point_runtime=double(),
#   conf_int_runtime=double()
# )

# Generate initial MC samples
u.init = runif(500, 0, 1)

# Estimate optimal values for B
B.minmi = data.frame(
  error_factor = error_factors[-1],
  point = sapply(error_factors[-1],
                 FUN = function(x) {
                   find_optimal_B(
                     q = 0.5,
                     max_var = 0.2 * (mean(fossil.sd)) ^ 2,
                     K = K,
                     m = min(datasets[1, "W"][[1]]),
                     n = n,
                     u = u.init,
                     eps.mean = 0,
                     eps.sigma = x * mean(fossil.sd)
                   )
                 }),
  lower = sapply(error_factors[-1],
                 FUN = function(x) {
                   find_optimal_B(
                     q = 0.025,
                     max_var = 0.2 * (mean(fossil.sd)) ^ 2,
                     K = K,
                     m = min(datasets[1, "W"][[1]]),
                     n = n,
                     u = u.init,
                     eps.mean = 0,
                     eps.sigma = x * mean(fossil.sd)
                   )
                 }),
  upper = sapply(error_factors[-1],
                 FUN = function(x) {
                   find_optimal_B(
                     q = 0.975,
                     max_var = 0.2 * (mean(fossil.sd)) ^ 2,
                     K = K,
                     m = min(datasets[1, "W"][[1]]),
                     n = n,
                     u = u.init,
                     eps.mean = 0,
                     eps.sigma = x * mean(fossil.sd)
                   )
                 })
)
B.minmi = rbind(B.minmi, c(0, 2, 2, 2))

start_time = Sys.time()
for (i in 1:nrow(datasets)) {
  iter = datasets[i, ]
  W = as.numeric(iter$W[[1]])
  sd = as.numeric(iter$error_factor * fossil.sd)
  print(paste("Dataset ID:", iter$id))
  
  print("- Point Estimates")
  for (method in methods.point_estimates) {
    print(paste("--", method))
    estimation = estimate_extinction(
      W = W,
      sd = sd,
      method = method,
      K = K,
      dating_error.mean = dating_error.mean
    )
    results = tibble::add_row(
      results,
      id = iter$id,
      error_factor = iter$error_factor,
      method = method,
      point = estimation$point,
      point_runtime = estimation$point_runtime,
      conf_int_runtime = estimation$conf_int_runtime
    )
  }
  
  print("- Confidence Intervals")
  for (method in methods.conf_int) {
    print(paste("--", method))
    estimation = estimate_conf_int(
      W = W,
      sd = sd,
      method = method,
      alpha = alpha,
      K = K,
      dating_error.mean = dating_error.mean,
      B.point = B.minmi[B.minmi$error_factor == iter$error_factor, "point"],
      B.lower = B.minmi[B.minmi$error_factor == iter$error_factor, "lower"],
      B.upper = B.minmi[B.minmi$error_factor == iter$error_factor, "upper"]
    )
    results = tibble::add_row(
      results,
      id = iter$id,
      error_factor = iter$error_factor,
      method = method,
      lower = estimation$lower,
      point = estimation$point,
      upper = estimation$upper,
      point_runtime = estimation$point_runtime,
      conf_int_runtime = estimation$conf_int_runtime,
      B.point = estimation$B.point,
      B.lower = estimation$B.lower,
      B.upper = estimation$B.upper
    )
  }
}
sims.time = Sys.time() - start_time
print(sims.time)

View(results)
saveRDS(results, "../data/sim_exp-estimate_extinction_results.RDS")
