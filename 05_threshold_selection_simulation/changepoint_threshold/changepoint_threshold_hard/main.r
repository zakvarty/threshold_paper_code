##  changepoint_threshold_hard/main.r
###
## Selecting a changepoint threshold when censoring increases exponentially below
## the modelling threshold.
##
## Can be run on locally for a single catalogue by manually setting the jobid,
## or on STORM for an array of jobids (simulated catalogues).
##
###
## Author: Zak Varty
##
START_TIME <- Sys.time()
#_______________________________________________________________________________
### Get the jobid from the command line argument
jobid  <- as.numeric(commandArgs(trailing=TRUE))
#jobid <- 1 #(for single job / running locally)
print(paste0("This is changepoint_threshold_hard/main.r, where jobid = ", jobid, "."))

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

source("./src/thresholding/sample_latent_magnitudes.R")
source("./src/thresholding/standardise_gpd_sample.R")

source("./src/get_mle_and_bootstrap.R")
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


#_______________________________________________________________________________
###
## Set seeds for each catalogue (will use 1-500) ----
###
set.seed(1234)
seeds <- sample(1:1e7, size = 1e3, replace = FALSE)
print("seeds 1:10 are:")
seeds[1:10]

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

# Save simulation paramters for use when plotting
rds_path <- './Output/data/true_parameter_values.RData'
save(sig_u, xi, u, v, sig_0, n, to_nearest, file =  rds_path)

# Range of catalogues to run in this job
cat_range <- (jobid - 1) * 50 + 1:50

#_______________________________________________________________________________
###
## Function that calculates reciporical of d(q,1) for a proposed threshold ----
###
FUN_3 <- function(threshold_1, threshold_2, change_location){
  step_length_1 <- floor(change_location)
  step_length_2 <- n_tot - step_length_1
  threshold <- c(rep(threshold_1, step_length_1), rep(threshold_2, step_length_2))

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

#_______________________________________________________________________________
###
##  ----
###
for (cat_id in cat_range){
  try(expr = {
    print("#####")
    print(paste0("### STARTING cat_id = ", cat_id ))
    print("#####")
    ###
    ## Simulate rounded GPD data ----
    ###
    # set the seed for this catalogue
    set.seed(seeds[cat_id])

    # Simulate catalogue (reproducibly)
    m_raw <- rgpd(n = n, scale = sig_u, shape = xi, mu = u)
    cat_full <- data.frame(
      m_raw = m_raw,
      u = rep(u, n),
      v = rep(v, each = n/2),
      index = seq_along(m_raw))

    # reduce to rounded, censored catalogue
    cat <- cat_full %>%
      filter(m_raw >= v) %>%
      mutate(m_rounded = round_to_nearest(m_raw, to_nearest),
             index = seq_along(m_rounded)) %>%
      select(m = m_rounded, u, v, index)

    ###
    ## Save catalouge size and  true step lengths ----
    ###
    n_tot <- NROW(cat)
    n_1 <- sum(cat$v == v[1])
    n_2 <- sum(cat$v == v[2])

    catalogue_size <- tibble(cat_id = cat_id, n_total = n_tot, n_1 = n_1, n_2 = n_2)
    rds_path <- paste0("./Output/data/catalogue_sizes/catalogue_size_",cat_id,".RDS")
    saveRDS(object = catalogue_size, file = rds_path)

    ###
    ## Set up Bayesian optimization ----
    ###
    # Bounds on each threshold parameter
    bounds <- list(
      threshold_1 = c(0,1.5),
      threshold_2 = c(0,1.5),
      change_location = c(100,1000))

    # Initial evaluation points (chosen randomly)
    n_init <- 20
    initGrid <- data.frame(
      threshold_1 = runif(n = n_init, min = bounds[[1]][1], max = bounds[[1]][2]),
      threshold_2 = runif(n = n_init, min = bounds[[2]][1], max = bounds[[2]][2]),
      change_location = runif(n = n_init, min = bounds[[3]][1], max = bounds[[3]][2]))

    # Output container
    opt_object <- bayesOpt(
      FUN = FUN_3
      , bounds = bounds
      , initGrid = initGrid
      , iters.n = 1
      , plotProgress = FALSE
      , verbose = 2
      , acq = "ei"
    )

    ###
    ## Do Bayesian optimization ----
    ###
    iterations <- 79
    opt_object <- addIterations(opt_object, iters.n = iterations, verbose = 2, plotProgress = FALSE)

    # save output
    rds_path <- paste0("./Output/data/opt_objects/opt_object_",cat_id,".RDS")
    saveRDS(object = opt_object, file = rds_path)

    # save best parameters found
    rds_path <- paste0("./Output/data/selected/selected_threshold_parameters_",cat_id,".RDS")
    saveRDS(object = as.data.frame(getBestPars(opt_object)), file = rds_path)

    print(paste0("DONE cat_id = ", cat_id))
  })
}

print(paste0("DONE jobid = ", jobid))
