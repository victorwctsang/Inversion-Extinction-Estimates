# Simulation Experiments
# Approaches:
# - New Method
# - Garthwaite Robbins Munro
# Metrics:
# - 95% confidence intervals
#   - Coverage
#   - Interval width
#   - Computation time
#   - Above, but as error increases

library(tidyverse)
source("helpers-simulation-experiments.R")
source("new-method-functions.R")
source("garthwaite-robbins-munro-functions.R")

################################################################################
# Experiment Setup
################################################################################

set.seed(2022)

n.sims = 300                      # Number of simulated datasets
n = 30                            # Number of fossil samples in each dataset
theta.true = 10000                # True extinction date
K = 20000                         # Upper bound for fossil ages

dating_error.mean = 0             # True mean of radiocarbon dating error
dating_error.sd = getStdDevFromFossilData(path='../data/fossildata.xlsx',
                                          K=K,
                                          sheet="MammothPrimEBer", 
                                          range="M3:N36", 
                                          col_names=c("age", "sd"), 
                                          col_types=c('numeric', 'numeric'))

alpha = 0.05                      # 95% confidence interval
B = 500                           # Number of monte carlo samples
uniroot.interval = c(2000, 19500) # TODO uniroot search interval
max_iter = 2000                   # Number of iterations for Garthwaite RM

methods = c("new_method", "rm")           # Methods to use
error_factors = c(0, 1, 2, 4)     # Different multiples to apply to dating error sd

################################################################################
# Initialise variables
################################################################################

u = matrix(runif(n*B, min=0, max=1), ncol=B)                 # Base MC Samples

result.df = create_result_df(n.sims, methods, error_factors)
sim.datasets = simulate_datasets(n.sims, error_factors, theta.true, K, dating_error.mean, dating_error.sd, n)

#####################
###################################################################################
# Run experiments########################################################

for (i in 1:nrow(result.df)) {
  row = result.df[i, ]
  W = sim.datasets[(sim.datasets$sim_id == row$sim_id) & (sim.datasets$error_factor == row$error_factor), "W"][[1]]
  start_time = Sys.time()
  if (row$method == "new_method") {
    row[c("lower", "upper")] = estimate_CI.mc(alpha, n, K, W, u, dating_error.mean, (row$error_factor * dating_error.sd), uniroot.interval)
  } else if (row$method == "rm") {
    row[c("lower", "upper")] = estimate_CI.rm(W, K, alpha, max_iter, dating_error.mean, (row$error_factor * dating_error.sd), return_iters=FALSE)
  }
  row$compute_time = Sys.time() - start_time
  row$mean = mean(c(row$lower, row$upper))
  row$width = row$upper - row$lower
  row$contains_theta = (row$lower < theta.true) & (theta.true < row$upper)
  result.df[i, ] = row
}

View(result.df)

result.df %>%
  group_by(method, error_factor) %>%
  summarise(coverage = mean(contains_theta),
            avg.width = mean(width),
            avg.time = mean(compute_time))
