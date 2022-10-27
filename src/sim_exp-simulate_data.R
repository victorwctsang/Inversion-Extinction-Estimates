# Simulation Experiments
# Generate synthetic datasets for simulation experiments

library(tidyverse)
source("helpers-simulation-experiments.R")

options(nwarnings = 10000) 

# Configure synthetic data generation

synthetic.data.config = list(
  seed = 2022,
  theta.true = 10000,
  K = 20000,
  n.datasets = 1000,
  n = 20,
  error_factors = c(0, 1, 2, 4),
  n.error_factors = length(error_factors),
  dating_error.mean = 0,
  fossil.sd = NA
)

synthetic.data.config$fossil.sd = getStdDevFromFossilData(
  path='../data/fossildata.xlsx',
  sheet="MammothPrimEBer", 
  range="M3:N36", 
  col_names=c("age", "sd"), 
  col_types=c('numeric', 'numeric'),
  K=synthetic.data.config$K,
  n=synthetic.data.config$n
)

attach(synthetic.data.config)

# Generate synthetic data
set.seed(seed)
datasets = simulate_datasets(n.datasets, error_factors, theta.true, K, dating_error.mean, mean(fossil.sd), n)

# Export config and datasets
save(synthetic.data.config, datasets, file = "../data/synthetic-data.RData")