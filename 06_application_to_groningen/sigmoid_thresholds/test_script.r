get_sigmoid_values <- function(tau, mu, sigma, v_1, v_2){
  v_1 + pnorm((tau - mu)/ sigma) * (v_2 - v_1)
}
#_______________________________________________________________________________
###
## Source required scripts ----
###


source("./src/_gpd.R")
source("./src/round_to_nearest.R")
source("./src/llh_gpd_rd_varu.R")
source("./src/mle_gpd_rd_varu.R")
source("./src/sample_mles_gpd_rd_stepped.R")
source('src/bootstrapping.R')
source("src/bootstrap_catalogues_sigmoid.r")
source("./src/thresholding/sample_latent_magnitudes.R")
source("./src/thresholding/standardise_gpd_sample.R")

source("./src/get_mle_and_bootstrap.R") # Need to adapt for sigmoid ?
source("./src/get_dq1.R")

print("functions sourced.")

#_______________________________________________________________________________
###
## Load required packages ----
###

library("purrr")
library("dplyr")
library("ggplot2")
library("ParBayesianOptimization")

## Set true parameters
#_______________________________________________________________________________
###
## Set true parameters and simulate rounded GPD data ----
###

# True parameters
sig_u <- 0.55
xi  <- - 0.1
u   <- 0
v <- c(0.83,0.42)
sig_0 <- sig_u + (0 - u) * xi
n   <- 4000
to_nearest <- 0.1
detection_rate <- 7

mu <- 1500
sigma <- 300
v_1 <- 0.83
v_2 <- 0.42
v_tau <- get_sigmoid_values(tau = 1:n, mu, sigma, v_1, v_2)

plot(1:n,v_tau, type = 'l')

## Save simulation paramters for use when plotting
#rds_path <- './Output/data/true_parameter_values.RData'
#save(sig_u, xi, u, v, sig_0, n, to_nearest, file =  rds_path)

set.seed(2345)

m_raw <- rgpd(n = n, scale = sig_u, shape = xi, mu = u)
cat_full <- data.frame(
  m_raw = m_raw,
  u = rep(u, n),
  v = v_tau,
  index = seq_along(m_raw))

# reduce to rounded, censored catalogue
cat_phased <- cat_full %>%
  mutate(distance_below = pmax(v - m_raw, 0),
         detection_prob = pexp(q = distance_below ,rate = detection_rate, lower.tail = FALSE),
         u = runif(n),
         keep = u < detection_prob) %>%
  filter(keep) %>%
  mutate(m_rounded = round_to_nearest(m_raw, to_nearest),
         index = seq_along(m_rounded)) %>%
  select(m = m_rounded, u, v, index)

cat_hard <- cat_full %>%
  filter(m_raw >= v) %>%
  mutate(m_rounded = round_to_nearest(m_raw, to_nearest),
         index = seq_along(m_rounded)) %>%
  select(m = m_rounded, u, v, index)

# Make some plots
plot(cat_full$index, cat_full$m_raw)
lines(cat_full$index, cat_full$v, col = 2)

plot(cat_hard$index, cat_hard$m)
lines(cat_hard$index, cat_hard$v, col = 2)

plot(cat_phased$index, cat_phased$m)
lines(cat_phased$index, cat_phased$v, col = 2)

# set cat to be tho one you want
cat <- cat_hard

# record "true" threshold parameters
n_tot <- NROW(cat)
mu_true <- min(which(cat$v < mean(c(v_1,v_2))))

catalogue_size <- tibble(cat_id = cat_id, n_total = n_tot, n_1 = n_1, n_2 = n_2)
rds_path <- paste0("./Output/data/catalogue_sizes/catalogue_size_",cat_id,".RDS")
#saveRDS(object = catalogue_size, file = rds_path)



#_______________________________________________________________________________
###
## Function that calculates reciporical of d(q,1) for a proposed threshold ----
###

FUN_4 <- function(mu, sigma, v_1, v_2, plot = FALSE){

  threshold <- get_sigmoid_values(tau = 1:NROW(cat), mu, sigma, v_1, v_2)
  #threshold <- get_sigmoid_values(tau = 1:NROW(cat), 700, 100, 1.1, 0.7)

  mles <- get_mle_and_bootstrap(
    cat = cat,
    threshold = threshold,
    to_nearest = 0.1,
    u = 0,
    n_bootstrap_mles = 750,
    verbose = TRUE)

  dq1_value <- get_dq1(
    cat = cat,
    threshold = threshold,
    to_nearest = 0.1,
    u =  0,
    mle_obj = mles,
    n_eval_pts = 501,
    y_samples_per_mle = 1)

  return(list(Score = 1/dq1_value))
}

FUN_4_visual <- function(mu, sigma, v_1, v_2){

  threshold <- get_sigmoid_values(tau = 1:NROW(cat), mu, sigma, v_1, v_2)
  #threshold <- get_sigmoid_values(tau = 1:NROW(cat), 700, 100, 1.1, 0.7)

  mles <- get_mle_and_bootstrap(
    cat = cat,
    threshold = threshold,
    to_nearest = 0.1,
    u = 0,
    n_bootstrap_mles = 750,
    verbose = TRUE)

  dq1_value <- get_dq1(
    cat = cat,
    threshold = threshold,
    to_nearest = 0.1,
    u =  0,
    mle_obj = mles,
    n_eval_pts = 501,
    y_samples_per_mle = 1)

    plot(cat$index, cat$m, main = paste0("d(q,1) = ", round(dq1_value,5)))
    lines(cat$index, threshold, col = 2, lwd = 2)

  return(list(Score = 1/dq1_value))
}
FUN_4(mu = 300, sigma = 50, v_1 = 0.85, v_2 = 0.45,plot = TRUE)
FUN_4_visual(mu = 150, sigma = 500, v_1 = 1.34, v_2 = 0.43)

#_______________________________________________________________________________
###
## Set up Bayesian optimization ----
###

# Bounds on each threshold parameter
bounds <- list(
  v_1 = c(0.3,1.4),
  v_2 = c(0.3,1.4),
  mu = c(150, 1250),
  sigma = c(1,500))

# Initial evaluation points (chosen randomly)
n_init <- 20
initGrid <- data.frame(
  v_1 = runif(n = n_init, min = bounds[[1]][1], max = bounds[[1]][2]),
  v_2 = runif(n = n_init, min = bounds[[2]][1], max = bounds[[2]][2]),
  mu  = runif(n = n_init, min = bounds[[3]][1], max = bounds[[3]][2]),
  sigma = runif(n = n_init, min = bounds[[4]][1], max = bounds[[4]][2]))

# Output container
opt_object <- bayesOpt(
  FUN = FUN_4
  , bounds = bounds
  , initGrid = initGrid
  , iters.n = 1
  , plotProgress = TRUE
  , verbose = 2
  , acq = "ei"
)

###
## Do Bayesian optimization ----
###
iterations <- 20
opt_object <- addIterations(opt_object, iters.n = iterations, verbose = 2, plotProgress = TRUE)

# save output
rds_path <- paste0("./opt_object.RDS")
#saveRDS(object = opt_object, file = rds_path)
opt_object <- readRDS(rds_path)

opt_pars <- unlist(ParBayesianOptimization::getBestPars(opt_object))
FUN_4_visual(mu = opt_pars[3], sigma = opt_pars[4], v_1 = opt_pars[1], v_2 = opt_pars[2])
lines(cat$index, cat$v, col = 4, lwd = 2)

FUN_4_visual(mu = 300, sigma = 90, v_1 = 0.83, v_2 = 0.42)
lines(cat$index, cat$v, col = 4, lwd = 2)
