
source("new-method-functions.R")

set.seed(2022)

n = 30                            # Number of fossil samples in each dataset
theta.true = 10000                # True extinction date
K = 20000                         # Upper bound for fossil ages

dating_error.mean = 0             # True mean of radiocarbon dating error
dating_error.sd = 700             # True Standard deviation

alpha = 0.05                      # 95% confidence interval
B = 500                           # Number of monte carlo samples
uniroot.interval = c(2000, 19500) # TODO uniroot search interval

u=matrix(runif(n*B, 0, 1), ncol=B)# Base MC Samples

# Simulate some data
W = runif(n=n, min=theta.true, max=K) + rnorm(n=30, mean=dating_error.mean, sd=dating_error.sd)
m=min(W)

# Get an initial estimate of upper and lower bounds
init.CI = estimate_CI.mc(alpha, n, K, W, u, dating_error.mean, dating_error.sd, uniroot.interval)
init.CI

# Use this initial estimate to generate n samples of measurement error
e = list()
e$lower = rtnorm(n, mean=dating_error.mean, sd=dating_error.sd, a=-Inf, b=m-init.CI$CI.lower)
e$upper = rtnorm(n, mean=dating_error.mean, sd=dating_error.sd, a=-Inf, b=m-init.CI$CI.upper)
e

# Estimate sigma.phi
sigma.phi = list()
sigma.phi$lower = sd((K - m)/(K-e$lower-init.CI$CI.lower))
sigma.phi$upper = sd((K - m)/(K-e$upper-init.CI$CI.upper))
sigma.phi

# Find optimal B for a variance of 10% of our sampling error
optimal_B = list()
optimal_B$lower = find_optimal_B(target_sd=0.2*dating_error.sd,
                                 q=alpha/2,
                                 theta=init.CI$CI.lower,
                                 K, n, m,
                                 sd.phi=sigma.phi$lower,
                                 dating_error.mean,
                                 dating_error.sd)

optimal_B$upper = find_optimal_B(target_sd=0.2*dating_error.sd,
                                 q=(1 - alpha/2),
                                 theta=init.CI$CI.upper,
                                 K, n, m,
                                 sd.phi=sigma.phi$upper,
                                 dating_error.mean,
                                 dating_error.sd)
optimal_B

