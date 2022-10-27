# Simulation Experiments

library(tidyverse)

load("../data/synthetic-data.RData")
attach(synthetic.data.config)

estimates = readRDS("../data/sim_exp-estimate_extinction_results.RDS")
View(estimates)

# Point estimates
performance.point_estimates = estimates %>%
  group_by(error_factor, method) %>%
  summarise(MSE_000 = mean((point - theta.true)^2)/1000,
            bias = mean(point)-theta.true,
            variance_000 = var(point)/1000,
            avg_runtime = mean(point_runtime))

View(performance.point_estimates)

# Confidence Intervals
performance.conf_int_estimates = estimates %>%
  filter(!is.na(conf_int_runtime)) %>%
  mutate(width = upper - lower,
         contains_theta = ifelse(theta.true > lower & theta.true < upper, 1, 0)) %>%
  group_by(error_factor, method) %>%
  summarise(coverage = mean(contains_theta),
            avg_width = mean(width),
            avg_runtime = mean(conf_int_runtime))

View(performance.conf_int_estimates)
