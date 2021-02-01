##  sigmoid_Groningen/main.r
###
## Selecting a sigmiod threshold for Gronignen earthquake catalogue.
## Can be run on locally or STORM.
###
## Author: Zak Varty
##

#_______________________________________________________________________________
print("This is changepoint_threshold_Groningen/main.r")


################################################################################
################################################################################
###########         Part 0: Set up                             #################
################################################################################
################################################################################
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
source("./src/get_sigmoid_values.r")
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
## Set seed
###
set.seed(1234)
seeds <- sample(x = 1:1e7,size = 1e3,replace = FALSE)
#_______________________________________________________________________________
###
## Load Groningnen catalogue ----
###

gron_cat_raw <- readr:::read_csv("./Input/data/2020-05-14_14-46-48_cat.csv")

gron_cat <- gron_cat_raw %>%
  filter(date > as.Date("1995-04-01")) %>%
  filter(date < as.Date("2020-01-01")) %>%
  select(date, time, julian, mag) %>%
  arrange(julian) %>%
  mutate(index = seq_along(date))

rds_path <- "./Output/data/gron_cat.RDS"
saveRDS(gron_cat, rds_path)

cat <- gron_cat %>%
  mutate(m = mag) %>%
  select(index,m)

rds_path <- "./Output/data/cat.RDS"
saveRDS(cat, rds_path)

################################################################################
################################################################################
###########         Part 1: Joint parameter selection          #################
################################################################################
################################################################################
#_______________________________________________________________________________
###
## Function that calculates reciporical of d(q,1) for a proposed threshold ----
###

# Calculate d(q,1) for cat and 4-parameter sigmoid threshold
FUN_4 <- function(mu, sigma, v_1, v_2, plot = FALSE){

  threshold <- get_sigmoid_values(tau = 1:NROW(cat), mu, sigma, v_1, v_2)

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

# Copy of FUN_4, that also plots the threshold and shows metric value
FUN_4_visual <- function(mu, sigma, v_1, v_2){

  threshold <- get_sigmoid_values(tau = 1:NROW(cat), mu, sigma, v_1, v_2)

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

# Examples
#FUN_4(mu = 800, sigma = 50, v_1 = 0.85, v_2 = 0.45,plot = TRUE)
#FUN_4_visual(mu = 800, sigma = 200, v_1 = 1.7, v_2 = 1.7)

#_______________________________________________________________________________
###
## Bayesian optimization ----
###

# Bounds on each threshold parameter
bounds <- list(
  v_1 = c(0.4,1.7),
  v_2 = c(0.4,1.7),
  mu = c(200, 1100),
  sigma = c(1,500))

# Set number of initialisation point to run
n_inits = 5
saveRDS(n_inits, "./Output/data/n_inits.RDS")

# DO bayesian optimisation from n_inits random starting values
for(run_number in 1:n_inits){
  # Set initial evaluation points (chosen randomly)
  set.seed(seeds[run_number])
  n_init <- 5

  initGrid <- data.frame(
    v_1 = runif(n = n_init, min = bounds[[1]][1], max = bounds[[1]][2]),
    v_2 = runif(n = n_init, min = bounds[[2]][1], max = bounds[[2]][2]),
    mu  = runif(n = n_init, min = bounds[[3]][1], max = bounds[[3]][2]),
    sigma = runif(n = n_init, min = bounds[[4]][1], max = bounds[[4]][2])
    )

  # Create output container
  opt_object <- bayesOpt(
    FUN = FUN_4
    , bounds = bounds
    , initGrid = initGrid
    , iters.n = 1
    , plotProgress = TRUE
    , verbose = 2
    , acq = "ei"
  )

  # add iterations to the Bayesian optimization
  iterations <- 89
  opt_object <- addIterations(opt_object, iters.n = iterations, verbose = 2, plotProgress = TRUE)

  # Save  the output
  rds_path <- paste0("./Output/data/opt_object_",run_number,".RDS")
  saveRDS(object = opt_object, file = rds_path)
}
#_______________________________________________________________________________
###
## Get the top n_ranks threshold parmameters from each BO initialisation ----
###

n_ranks = 5
saveRDS(n_ranks, "./Output/data/n_ranks.RDS")

top_thresholds <- data.frame(
  init = rep(1:n_inits, each = n_ranks),
  rank = rep(1:n_ranks, n_inits),
  v_1 = rep(NA, n_inits * n_ranks),
  v_2 = rep(NA, n_inits * n_ranks),
  mu = rep(NA, n_inits * n_ranks),
  sigma = rep(NA, n_inits * n_ranks),
  dq1 = rep(NA, n_inits * n_ranks)
)
for( i in 1:5){
  opt_object_loaded <- readRDS(file = paste0('Output/data/opt_object_',i,'.RDS'))
  best_pars_loaded <- opt_object_loaded$scoreSummary %>%
    arrange(-Score) %>%
    mutate(dq1 = 1/ Score) %>%
    select(v_1, v_2,mu, sigma, dq1) %>%
    head(n_ranks)
  top_thresholds[(i-1)*n_ranks + 1:5, 3:7] <- best_pars_loaded
}
saveRDS(top_thresholds, "./Output/data/top_thresholds.RDS")

#_______________________________________________________________________________
###
##  Rerun the chosen thresholds from each initialisation point --------
###
n_reruns = 10
saveRDS(n_reruns, './Output/data/n_reruns.RDS')
set.seed(1234)
rerun_values <- data.frame(
  init = rep(1:n_inits, each = n_reruns),
  rerun = rep(1:n_reruns, n_inits),
  v_1 = rep(NA, n_inits * n_reruns),
  v_2 = rep(NA, n_inits * n_reruns),
  mu = rep(NA, n_inits * n_reruns),
  sigma = rep(NA, n_inits * n_reruns),
  dq1 = rep(NA, n_inits * n_reruns)
)

for( i in 1:NROW(rerun_values)){
  current_init <- rerun_values$init[i]
  current_row <- which(top_thresholds$init == current_init & top_thresholds$rank == 1)
  current_params <- top_thresholds[current_row,]
  rerun_value <- FUN_4(
    mu = current_params$mu,
    sigma = current_params$sigma,
    v_1 = current_params$v_1,
    v_2 = current_params$v_2)
  rerun_values[i,3:7] <-  c(current_params[c(3:6)], 1/rerun_value$Score)
  print(paste0(i,"/",NROW(rerun_values)," done."))
}

rds_path <- "./Output/data/rerun_values_selected_thresholds.RDS"
saveRDS(rerun_values, rds_path)

#_______________________________________________________________________________
###
##  Rerun the chosen thresholds and runners up from each init ----
###

# set up storage
temp <- expand.grid(init = 1:n_inits, rank = 1:n_ranks, rerun = 1:n_reruns)
rerun_values <- data.frame(
  init = rep(NA, NROW(temp)),
  rank = rep(NA, NROW(temp)),
  rerun = rep(NA, NROW(temp)),
  v_1 = rep(NA, NROW(temp)),
  v_2 = rep(NA, NROW(temp)),
  mu = rep(NA, NROW(temp)),
  sigma = rep(NA, NROW(temp)),
  dq1 = rep(NA, NROW(temp))
)
rerun_values[,1:3] <- temp

# Do re-evaluations
set.seed(4321)
for( i in 1:NROW(rerun_values)){
  current_init <- rerun_values$init[i]
  current_rank <- rerun_values$rank[i]
  current_id <- which(top_thresholds$init == current_init & top_thresholds$rank == current_rank)
  current_params <- top_thresholds[current_id,]
  rerun_value <- FUN_4(
    mu = current_params$mu,
    sigma = current_params$sigma,
    v_1 = current_params$v_1,
    v_2 = current_params$v_2)
  rerun_values[i,4:8] <-  c(current_params[c(3:6)], 1/rerun_value$Score)
  print(paste0(i,"/",NROW(rerun_values)," done."))
}

# Save the rerun values
rds_path <- "./Output/data/rerun_values_good_thresholds.RDS"
saveRDS(rerun_values, rds_path)


################################################################################
################################################################################
###########         Part 2: fixing end levels           ########################
################################################################################
################################################################################

## Since ending levels are stable across all selected thresholds, fix those
## and try to fix  mu and sigma.

## Part 2 can be re-run independently of Part 1 but needs the set-up in Part 0
## and the definiton of FUN_4.


#_______________________________________________________________________________
###
## Restrict the function to have only 2 parameters that are optimised over.
###

FUN_4 <- function(mu, sigma, v_1, v_2, plot = FALSE){

  threshold <- get_sigmoid_values(tau = 1:NROW(cat), mu, sigma, v_1, v_2)

  mles <- get_mle_and_bootstrap(
    cat = cat,
    threshold = threshold,
    to_nearest = 0.1,
    u = 0,
    n_bootstrap_mles = 1500,
    verbose = TRUE)

  dq1_value <- get_dq1(
    cat = cat,
    threshold = threshold,
    to_nearest = 0.1,
    u =  0,
    mle_obj = mles,
    n_eval_pts = 501,
    y_samples_per_mle = 10)

  return(list(Score = 1/dq1_value))
}


FUN_4_restricted <- function(mu, sigma){
  FUN_4(mu, sigma, v_1 = 1.15, v_2 = 0.76)
}


#_______________________________________________________________________________
###
## Optimise over mu, sigma with v_1 and v_2 fixed.
###

# Set bounds on each threshold parameter
bounds <- list(
  mu = c(400, 1100),
  sigma = c(1,300))
saveRDS(bounds, './Output/data/restricted_selection/bounds.RDS')

n_inits <- 5
saveRDS(n_inits, './Output/data/restricted_selection/n_inits.RDS')

for (run_number in 1:n_inits){

  # Choose random starting points
  set.seed(seeds[run_number])
  n_init <- 15
  initGrid <- data.frame(
    mu    = runif(n = n_init, min = bounds[[1]][1], max = bounds[[1]][2]),
    sigma = runif(n = n_init, min = bounds[[2]][1], max = bounds[[2]][2])
  )

  # Create output container
  opt_object <- bayesOpt(
    FUN = FUN_4_restricted
    , bounds = bounds
    , initGrid = initGrid
    , iters.n = 1
    , plotProgress = TRUE
    , verbose = 2
    , acq = "ei"
  )

  # Add iterations to Bayesian optimization
  iterations <- 50
  opt_object <- addIterations(opt_object, iters.n = iterations, verbose = 2, plotProgress = TRUE)

  # save output
  rds_path <- paste0("./Output/data/restricted_selection/opt_object_",run_number,".RDS")
  saveRDS(opt_object, file = rds_path)
}


good_mu_sigmas <- opt_object$scoreSummary %>%
  mutate(dq1 = 1/ Score) %>%
  arrange(dq1) %>%
  head(10) %>%
  select(mu, sigma, dq1)

plot(good_mu_sigmas$mu, good_mu_sigmas$sigma)

# save output
rds_path <- paste0("./Output/opt_object_",run_number,".RDS")
saveRDS(object = opt_object, file = rds_path)
opt_object <- readRDS(rds_path)

opt_pars <- unlist(ParBayesianOptimization::getBestPars(opt_object))
FUN_4_visual(mu = opt_pars[3], sigma = opt_pars[4], v_1 = opt_pars[1], v_2 = opt_pars[2])

#_______________________________________________________________________________
###
##  Rerun the chosen (restricted) thresholds and runners up from each init ----
###
n_reruns <- 10
# set up storage
temp <- expand.grid(init = 1:n_inits, rank = 1:n_ranks, rerun = 1:n_reruns)
rerun_values <- data.frame(
  init = rep(NA, NROW(temp)),
  rank = rep(NA, NROW(temp)),
  rerun = rep(NA, NROW(temp)),
  v_1 = rep(1.15, NROW(temp)),
  v_2 = rep(0.76, NROW(temp)),
  mu = rep(NA, NROW(temp)),
  sigma = rep(NA, NROW(temp)),
  dq1 = rep(NA, NROW(temp))
)
rerun_values[,1:3] <- temp

# Do re-evaluations
set.seed(4321)
for( i in 1:NROW(rerun_values)){
  current_init <- rerun_values$init[i]
  current_rank <- rerun_values$rank[i]
  current_id <- which(top_thresholds$init == current_init & top_thresholds$rank == current_rank)
  current_params <- top_thresholds[current_id,]
  rerun_value <- FUN_4_restricted(
    mu = current_params$mu,
    sigma = current_params$sigma
   )
  rerun_values[i,4:8] <-  c(current_params[c(3:6)], 1/rerun_value$Score)
  print(paste0(i,"/",NROW(rerun_values)," done."))
}

# Save the rerun values
rds_path <- "./Output/data/restricted_selection/rerun_values_good_thresholds_1500-10.RDS"
saveRDS(rerun_values, rds_path)

