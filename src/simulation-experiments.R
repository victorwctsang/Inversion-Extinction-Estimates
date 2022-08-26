# Simulation Experiments
# Approaches:
# - Measurement Error Monte Carlo Estimate
# - GRIWM
# - GRIWM Corrected
# Metrics:
# - Point Estimates:
#   - MSE
#   - Bias
#   - Variance
#   - L1 Norm
#   - Computational efficiency
# - 95% confidence intervals
#   - Coverage with varying error
#   - Interval width with varying error
#   - Computational efficiency

################################################################################
# Experiment Setup
################################################################################

library(tidyverse)
source("helpers.R")

set.seed(2022)

n.sims = 250
n = 30
theta.true = 10000
K = 20000


dating_error.mean = 0
dating_error.sd = getStdDevFromFossilData(path='../data/fossildata.xlsx',
                                          K=K,
                                          sheet="MammothPrimEBer", 
                                          range="M3:N36", 
                                          col_names=c("age", "sd"), 
                                          col_types=c('numeric', 'numeric'))

alpha = 0.05
B = 1000
uniroot.interval = c(5000, 25000) # TODO

error_factor = c(0, 1, 2, 4)
col.names = c("Lower", "Upper", "Mean", "Width", "Contains_Theta", "Time")
results.colnames = paste0(paste(col.names, rep(error_factor, each=length(col.names)), sep="."), "sigma")

results.CI.df = data.frame(matrix(NA, 
                                  ncol=length(col.names)*length(error_factor), 
                                  nrow=n.sims,
                                  dimnames = list(1:n.sims,
                                                  results.colnames)))

################################################################################
# Run experiments
################################################################################

# Simulate datasets with varying amounts of measurement error
sim.datasets.0 = lapply(1:n.sims, function (x) simulateFossils(n, theta.true, K, dating_error.mean, 0*dating_error.sd))
sim.datasets.1 = lapply(1:n.sims, function (x) simulateFossils(n, theta.true, K, dating_error.mean, 1*dating_error.sd))
sim.datasets.2 = lapply(1:n.sims, function (x) simulateFossils(n, theta.true, K, dating_error.mean, 2*dating_error.sd))
sim.datasets.4 = lapply(1:n.sims, function (x) simulateFossils(n, theta.true, K, dating_error.mean, 4*dating_error.sd))

# Generate MC Samples
u = matrix(runif(n*B, min=0, max=1), ncol=B)

# Compute confidence intervals

for (sim in 1:n.sims) {
  W = sim.datasets.0[[sim]]$W
  for (error_case_idx in 1:4) {
    start_col_idx = (error_case_idx-1)*6+1
    end_col_idx = start_col_idx+3
    contains_theta_col_idx = end_col_idx + 1
    time_col_idx = contains_theta_col_idx + 1

    start_time = Sys.time()
    results.CI.df[sim, start_col_idx:end_col_idx] = getConfidenceInterval(alpha, n, K, W, u, error_factor[error_case_idx]*dating_error.sd, dating_error.mean, uniroot.interval)
    results.CI.df[sim, contains_theta_col_idx] = (results.CI.df[sim, start_col_idx] < theta.true) & (theta.true < results.CI.df[sim, start_col_idx+1])
    results.CI.df[sim, time_col_idx] = Sys.time() - start_time
  }
}

View(results.CI.df)

results.CI.df %>%
  select(starts_with(c("Contains", "Time"))) %>%
  colMeans() %>%
  as.data.frame() %>%
  View()

results.CI.df %>%
  select(starts_with("Width")) %>%
  colMeans() %>%
  as.data.frame() %>%
  View()

################################################################################
# Log results
################################################################################

