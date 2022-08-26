
source("helpers.R")
library(readxl)


# Set parameters/assumptions
set.seed(2022)
K = 18000
B = 5000
alpha = 0.05
uniroot.interval = c(10000, 18000)

# Import mammoth data from fossildata.xlsx
fossil_dataframe = read_excel(
  path='../data/fossildata.xlsx', 
  sheet="MammothPrimEBer", 
  range="M3:N36", 
  col_names=c("age", "sd"), 
  col_types=c('numeric', 'numeric')
)
fossil_dataframe = fossil_dataframe[fossil_dataframe$age<K,]

# Define W, sigma, n, m
n = nrow(fossil_dataframe)
W = fossil_dataframe$age
sigma = mean(fossil_dataframe$sd)
m = min(W)

mc.samples = matrix(runif(n*B, min=0, max=1), ncol=B)

# Estimate confidence intervals

start.time = Sys.time()
getConfidenceInterval(alpha, n, K, W, u, sigma, 0, uniroot.interval)
Sys.time() - start.time
