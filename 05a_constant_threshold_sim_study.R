jobid <- 1
###
## Load required packages ----
###

library("purrr")
library("dplyr")
#library("ggplot2")
library("threshold")

# Set seeds for each catalogue
set.seed(1234)
seeds <- sample(1:1e7, size = 1e4, replace = FALSE)


seed <- seeds[jobid]

###
## Set true parameters and simulate rounded GPD data ----
###

sim_par <- list(
  u = 0,
  v = 0.32,
  scale_0 = 0.55,
  shape = -0.1,
  n = 1500,
  to_nearest = 0.1)



## Flat threshold - hard censoring



## Flat threshold - phased censoring
