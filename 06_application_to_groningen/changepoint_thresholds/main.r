##  changepoint_threshold_Groningen/main.r
###
## Selecting a changepoint threshold for Gronignen earthquake catalogue.
## Can be run on locally or STORM.
###
## Author: Zak Varty
##

#_______________________________________________________________________________
print("This is changepoint_threshold_Groningen/main.r")

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
## Set seed
###
set.seed(1234)


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

cat <- gron_cat %>%
  mutate(m = mag) %>%
  select(index,m)

plot(cat)

#_______________________________________________________________________________
###
## Function that calculates reciporical of d(q,1) for a proposed threshold ----
###
FUN_3 <- function(threshold_1, threshold_2, change_location){
  START_TIME <- Sys.time()

  step_length_1 <- floor(change_location)
  step_length_2 <- NROW(cat) - step_length_1
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
    y_samples_per_mle = 10)

  END_TIME <- Sys.time()
  print(END_TIME - START_TIME)
  return(list(Score = 1/dq1_value))
}

FUN_3_visual <- function(threshold_1, threshold_2, change_location, ...){
  START_TIME <- Sys.time()
  step_length_1 <- floor(change_location)
  step_length_2 <- NROW(cat) - step_length_1
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
    y_samples_per_mle = 10)

  plot(cat, main = paste0("d(q,1) =  ",round(dq1_value,5)))
  lines(x = cat$index, y = threshold, ...)

  END_TIME <- Sys.time()
  print(END_TIME - START_TIME)

  return(list(Score = 1/dq1_value))
}

# Sanity check that it works for previously problematic threshold
FUN_3_visual(threshold_1 = 1.41, threshold_2 = 1.41, change_location = 700, col = 2, lwd = 2)


## Selecting a flat threshold by grid-search:
fun3_eval_pts <- seq(0.6,1.9, by = 0.02)
fun3_vals <- rep(NA,length(fun3_eval_pts))
for(i in seq_along(fun3_eval_pts)){
  fun3_vals[i] <- FUN_3(threshold_1 = fun3_eval_pts[i],
                        threshold_2 = fun3_eval_pts[i],
                        change_location = 200)
}
plot(x = fun3_eval_pts, y = 1/unlist(fun3_vals))

rds_path <- "./Output/data/flat_threshold_grid_search/grid_points.RDS"
saveRDS(object = fun3_eval_pts, file = rds_path)
rds_path <- "./Output/data/flat_threshold_grid_search/dq1_evaluations.RDS"
saveRDS(object = fun3_vals, file = rds_path)

pdf_path <- "./Output/plots/flat_threshold_grid_search/dq1_evaluations.pdf"
pdf(pdf_path, width = 7, height = 5)
plot(x = fun3_eval_pts,
     y = 1/unlist(fun3_vals),
     xlab = 'threshold value',
     ylab = 'd(q,1)',
     main = "Groningen catalogue, flat threshold selection", ylim = c(0.05,0.1), pch = 16)
dev.off()

#_______________________________________________________________________________
###
##  Changepoint type threshold, selected by bayesian optimisaion ----
###

###
## Set up Bayesian optimization ----
###
# Bounds on each threshold parameter
bounds <- list(
  threshold_1 = c(0.4,1.8),
  threshold_2 = c(0.4,1.8),
  change_location = c(100,1100))

# Initial evaluation points (chosen randomly)
set.seed(1234)
seeds <- runif(n = 1000, min = 10000, max = 99999)
seed_id <- 5
set.seed(seeds[seed_id])

n_init <- 20

initGrid <- data.frame(
  threshold_1 = runif(n = n_init, min = bounds[[1]][1], max = bounds[[1]][2]),
  threshold_2 = runif(n = n_init, min = bounds[[2]][1], max = bounds[[2]][2]),
  change_location = runif(n = n_init, min = bounds[[3]][1], max = bounds[[3]][2]))

# Output container
opt_object <- bayesOpt(
  FUN = FUN_3
  , bounds = bounds
  #, initGrid = initGrid
  , initPoints = n_init
  , iters.n = 1
  , plotProgress = FALSE
  , verbose = 2
  , acq = "ei"
)

###
## Do Bayesian optimization ----
###
iterations <- 79
opt_object <- addIterations(opt_object, iters.n = iterations, verbose = 2, plotProgress = TRUE)

opt_object$scoreSummary %>%
  select(threshold_1, threshold_2, change_location, Score) %>%
  arrange(desc(Score)) %>%
  plot()

best_pars <- getBestPars(opt_object)
FUN_3_visual(
  threshold_1 = best_pars[[1]],
  threshold_2 =  best_pars[[2]],
  change_location = best_pars[[3]],
  col = 2,
  lwd = 2)
FUN_3_visual(
  threshold_1 = 1.45,
  threshold_2 =  1.45,
  change_location = best_pars[[3]],
  col = 2,
  lwd = 2)

# save output
rds_path <- paste0("./Output/data/opt_objects/opt_object_",seed_id,".RDS")
saveRDS(object = opt_object, file = rds_path)


## Plot all selected thresholds
best_values <- rep(NA,6)

pdf("./Output/plots/BO_selected_thresholds.pdf", width = 7, height = 5)
plot(cat, main = "5 BO runs with different initial sets", pch = 16, col = 'grey70', ylab = 'magnitude')
for( i in 1:6){
  opt_object_loaded <- readRDS(file = paste0('Output/data/opt_objects/opt_object_',i,'.RDS'))
  best_value_loaded <- 1/max(opt_object_loaded$scoreSummary$Score)
  best_values[i+1] <- best_value_loaded
  best_pars_loaded <- getBestPars(opt_object_loaded)
  best_value_loaded <-
    best_threshold <- c(
      rep(best_pars_loaded[[1]], floor(best_pars_loaded[[3]])),
      rep(best_pars_loaded[[2]], NROW(cat) -floor(best_pars_loaded[[3]]))
    )
  lines(cat$index, best_threshold, lwd = 2, col = i+2)
  abline(col = 9, h = bounds[[1]])
  abline(col = 9, v = bounds[[3]])
  legend("topleft", legend = 0:4, col = 2:6, pch = 16, title = "init")
}

plot(
  x = gron_cat$date,
  y = gron_cat$mag,
  main = "5 BO runs with different initial sets",
  xlab = "date",
  ylab = "magnitude",
  pch = 16,
  col = 'grey70')
for( i in 0:4){
  opt_object_loaded <- readRDS(file = paste0('Output/data/opt_objects/opt_object_',i,'.RDS'))
  best_value_loaded <- 1/max(opt_object_loaded$scoreSummary$Score)
  best_values[i+1] <- best_value_loaded
  best_pars_loaded <- getBestPars(opt_object_loaded)
  best_value_loaded <-
  best_threshold <- c(
    rep(best_pars_loaded[[1]], floor(best_pars_loaded[[3]])),
    rep(best_pars_loaded[[2]], NROW(cat) -floor(best_pars_loaded[[3]]))
    )
  lines(gron_cat$date, best_threshold, lwd = 2, col = i+2)
  abline(col = 9, h = bounds[[1]])
  abline(col = 9, v = gron_cat$date[bounds[[3]]])
}
legend("topleft", legend = 0:4, col = 2:6, pch = 16, title = "init")
dev.off()
best_values

rerun_values <- data.frame(
  opt_0 = rep(NA,10),
  opt_1 = rep(NA, 10),
  opt_2 = rep(NA, 10) ,
  opt_3 = rep(NA, 10) ,
  opt_4 = rep(NA, 10))

for( i in 0:4){
  opt_object_loaded <- readRDS(file = paste0('Output/data/opt_objects/opt_object_',i,'.RDS'))
  best_pars_loaded <- getBestPars(opt_object_loaded)
  for(j in 1:NROW(rerun_values)){
    rerun_value <- FUN_3(
        threshold_1 = best_pars_loaded[[1]],
        threshold_2 = best_pars_loaded[[2]],
        change_location = best_pars_loaded[[3]]
        )
    rerun_values[j,i+1] <-  rerun_value[[1]]
    print(paste0(i,j," done"))
  }
}

## Plot variability in each d(q,1) for optimal thresholds on re-evaluation
## (Did we just get lucky when we tested that threshold?)
pdf("./Output/plots/BO_optima_comparison.pdf", width = 7, height = 5)
plot(
  y = 1/reshape2::melt(data = rerun_values)$value,
  x = rep(0:4, each = 10),
  col = rep(2:6, each = 10),
  xlab = "Best threshold from BO #",
  ylab = "d(q,1)",
  main = "Investigating Groningen BO selections"
)
points(x = 0:4, colMeans(1/rerun_values), col = 1, pch = "-", cex = 2)
points(x = 0:4, y = best_values, col = 1, pch = 16)
legend(
  'bottomright',
  legend = c("BO value", "reevaluations", "reeval mean"),
  pch = c(16,1,95),
  col = c(1,1,1))
dev.off()
