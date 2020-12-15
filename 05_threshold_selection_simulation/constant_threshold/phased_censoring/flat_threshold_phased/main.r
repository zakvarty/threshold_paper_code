##  flat_threshold_phased/main.r
###
## selects flat threshold using four metrics for each of 500 simulated catalgoues
## Catalogues are simulted with flat thresholds and phased censoring.
## Author: Zak Varty
##

#_______________________________________________________________________________
### Get the jobid from the command line argument
jobid  <- as.numeric(commandArgs(trailing=TRUE))
print(paste0("This is flat_threshold_phased/main.r, where jobid = ", jobid, "."))

#_______________________________________________________________________________
###
## Source required scripts ----
###
source("./src/_gpd.R")
source("./src/round_to_nearest.R")
source("./src/llh_gpd_rd_varu.R")
source("./src/mle_gpd_rd_varu.R")
source('src/bootstrapping.R')
source("./src/sample_mles_gpd_rd_stepped.R")

source("./src/thresholding/sample_latent_magnitudes.R")
source("./src/thresholding/standardise_gpd_sample.R")
print("functions sourced.")

#_______________________________________________________________________________
###
## Load required packages ----
###

library("purrr")
library("dplyr")
library("ggplot2")
# Set seeds for each catalogue
set.seed(1234)
seeds <- sample(1:1e7, size = 1e4, replace = FALSE)
print("seeds 1:10 are:")
seeds[1:10]

#_______________________________________________________________________________
###
## ## Use jobid to retrieve and set the seed for this run ----
###
seed <- seeds[jobid]
seed_index <- jobid
set.seed(seed)
print(paste0("The seed for jobid ", jobid, " is ", seed, "."))

#_______________________________________________________________________________
###
## Set true parameters and simulate rounded GPD data ----
###

# True parameters
sig_u <- 0.6
xi  <- -0.03
u <- 0.0
v <- 0.32
sig_0 <- sig_u + (0 - u) * xi
n   <- 2040 # gives expected catalgue size of 1500
to_nearest <- 0.1
censoring_rate <- 7

# Simulate unrounded magnitudes as exceedances of u (reproducibly)
m_raw <- rgpd(n = n, scale = sig_u, shape = xi, mu = u)
cat_full <- data.frame(
  m_raw = m_raw,
  u = rep(u, n),
  v = rep(v, n),
  index = seq_along(m_raw))

# Apply exponential censoring below v_true
temp <- cat_full %>% mutate(
  distance_below_v = pmax(v - m_raw,0),
  keep_prob = pexp(q = distance_below_v, rate = 7, lower.tail = FALSE),
  U = runif(n),
  keep =  U < keep_prob)

temp2 <- temp %>%
  filter(keep == TRUE) %>%
  mutate(index = seq_along(m_raw))

# reduce to rounded, censored catalogue
cat <- temp2 %>%
  mutate(m_rounded = round_to_nearest(m_raw, to_nearest)) %>%
  select(m = m_rounded, u, v, index)
rm(temp, temp2)

print("head of cat")
head(cat)
print("cat has this many rows")
n_1 <- NROW(cat)
n_1

#_______________________________________________________________________________
###
## Set modelling thresholds to try ----
###

# Set thresholds
thresholds_vec <- c(seq(0.0, 1.0, by = 0.025))
thresholds <- map(.x = thresholds_vec, .f = function(x){rep(x,n_1)})

# pdf(file = paste0("./Output/plots/cat_plots/cat_",seed_index,".pdf"),
#      width = 7,
#      height =5)
#  # Plot catalogue and thresholds
# plot(
#   x = cat$index,
#   y = cat$m,
#   bty = 'n',
#   type = 'n',
#   ylim = c(0, 3),
#   main = "Flat threshold, phased censoring",
#   ylab = "magnitude",
#   xlab = "index"
# )
# abline(h = thresholds_vec, col = 'gray70', lwd = 0.3)
# points(cat$index, cat$m, pch = 16, cex = 0.8)
# abline(h = v, col = 1, lwd = 1)
# dev.off()

#rds_path <- paste0("./Output/data/cats/cat_",seed_index,".RDS")
#saveRDS(object = cat, file = rds_path)
#print("catalogue saved")

#_______________________________________________________________________________
###
## Subset data based on each threshold ----
###

filter_catalogue <- function(cat, threshold, to_nearest){
  cat %>%
    mutate( threshold = threshold) %>%
    filter(m  > (threshold - to_nearest/2 + 1e-8))
}

cats_filtered<- map(
  .x = thresholds,
  .f = filter_catalogue,
  cat = cat,
  to_nearest = to_nearest)

mags_filtered <- map(cats_filtered, function(df){df$m})

thresholds_filtered <- map(cats_filtered, function(df){df$threshold})

#_______________________________________________________________________________
###
## Calculate MLE point estimate at each threshold ----
###
print("finding mles")
mles <- pmap(
  .l = list(x = mags_filtered, v = thresholds_filtered),
  .f = mle_gpd_rd_varu,
  u = 0,
  sigxi = c(1,0),
  to_nearest = 0.1,
  llh_val = FALSE
)

#print("found mles, saving mles")
#rds_path <- paste0("./Output/data/mles/mles_",seed_index,".RDS")
#saveRDS(object = mles, file = rds_path)
#print("mles saved")

sigma_0_mles <- map_dbl(mles, function(vec){vec[1]})
xi_mles <- map_dbl(mles, function(vec){vec[2]})

#_______________________________________________________________________________
###
## Calculate expected number of exceedances at each threshold ----
###
print("calculating weights vectors")
w_vectors <- pmap(
  .f = get_w_vector,
  .l = list(
    x = mags_filtered,
    v = thresholds_filtered,
    gpd_mle = mles
  ),
  to_nearest = 0.1,
  u = 0)

#_______________________________________________________________________________
###
## Sample MLE values at each threshold ----
###
print("sampling MLE values at each threshold")
mle_samples <- vector(mode = "list", length = length(mles))
step_lengths <- map(.x = cats_filtered, NROW)
v_values <- c(v)
print('starting...')
for(i in seq_along(thresholds)){
  if(i %% 10 == 1){
  message("iteration: ", i, "/", length(thresholds))
  }
  mle_samples[[i]] <- sample_mles_gpd_rd_stepped(
    mle = mles[[i]],
    u = 0,
    step_lengths = step_lengths[[i]],
    step_values = thresholds_vec[[i]],
    w = w_vectors[[i]],
    to_nearest = 0.1,
    n_sim = 750,
    verbose = TRUE)
}
print("done.")
## Save / load samples MLE values
#print("Saving mle samples.")
#rds_path <- paste0("./Output/data/mle_samples/mle_samples_",seed_index,".RDS")
#saveRDS(mle_samples, file = rds_path)

## Create parameter stability plots over thresholds: set-up
sigma_0_samples <- map(mle_samples, function(vec){vec[1,]})
xi_samples      <- map(mle_samples, function(vec){vec[2,]})

#sigma_0_means <- map_dbl(sigma_0_mles, mean)
#sigma_0_high  <- map_dbl(sigma_0_samples, quantile, probs = 0.975)
#sigma_0_low   <- map_dbl(sigma_0_samples, quantile, probs = 0.025)

#xi_means <- map_dbl(xi_mles, mean)
#xi_high <-  map_dbl(xi_samples, quantile, probs = 0.975)
#xi_low  <-  map_dbl(xi_samples, quantile, probs = 0.025)

# # Paramter stability plots
# pdf_path <- paste0("output/qq_flat_threshold_sim_multiple/stability_plots/parameter_stability_",seed_index,".pdf")
# pdf(file = pdf_path, width = 8, height = 5)
# par(mfrow = c(1,2))
# # mean and CI
# plot(x = rep(thresholds_vec, 2),
#      y = c(sigma_0_high, sigma_0_low),
#      pch = "-",
#      xlab = "threshold value",
#      ylab = "sigma_0")
# points(x =thresholds_vec, y = sigma_0_means)
# plot(x = rep(thresholds_vec, 2),
#      y = c(xi_high, xi_low),
#      pch = "-",
#      xlab = "threshold value",
#      ylab = "xi")
# points(x =thresholds_vec, y = xi_means)
# dev.off()
# rm(xi_high, xi_low, xi_means, xi_mles, sigma_0_high, sigma_0_low, sigma_0_means, sigma_0_mles)

#_______________________________________________________________________________
###
## Calculate PP distance for each sampled mle ----
###
print("Calculating PP distance for each sampled mle.")
n_eval_pts <- 501
n_rep <- 10 ## latent samples per mle sample
pp_WMSE_samples <- vector(mode = "list", length = length(thresholds))
pp_WMAE_samples <- vector(mode = "list", length = length(thresholds))
CI_probs <- c(0.025,0.975)

# functions for in for loop
pp_distance_vector_std_exp <- function(x, eval_probs){
  eval_quantiles <- qexp(p = eval_probs, rate = 1)
  sample_probs <- ecdf(x)(eval_quantiles)
  return(sample_probs - eval_probs)
}
weighted_MSE <- function(d, w){mean(d^2 * w)}
weighted_MAE <- function(d, w){mean(abs(d) * w)}

print("starting...")
for(i in seq_along(thresholds)){

  ## Sample true magnitudes above that threshold
  latent_samples <- purrr::pmap(
    .l = list(sig_u = rep(mle_samples[[i]][1, ], n_rep),
              xi = rep(mle_samples[[i]][2, ], n_rep)),
    .f = sample_latent_magnitudes,
    x = cats_filtered[[i]]$m,
    v = cats_filtered[[i]]$threshold,
    u = 0,
    to_nearest = to_nearest
  )

  ## Transform to have common Exp(1) distribution
  latent_samples_exp <- purrr::pmap(
    .f = standardise_gpd_sample,
    .l = list(
      y = latent_samples,
      sig_u = rep(mle_samples[[i]][1,], n_rep),
      xi = rep(mle_samples[[i]][2,], n_rep)),
    v = cats_filtered[[i]]$threshold,
    u = 0
  )

  p <- ppoints(n_eval_pts)
  weights <- map(.x = latent_samples_exp,
                 .f = function(y, p){( p * (1-p) / sum(!is.na(y)) )^(-0.5)},
                 p = p)

  ## Calculate discrepancy from Exp(1) pp-plot
  pp_deviation_vectors <- map(
    .x = latent_samples_exp,
    .f = pp_distance_vector_std_exp,
    eval_probs = p
  )

  # Combine to weighted sum for each sample
  pp_WMSE <- pmap_dbl(
    .l = list(d = pp_deviation_vectors, w = weights),
    .f = weighted_MSE
  )

  pp_WMAE <- pmap_dbl(
    .l = list(d = pp_deviation_vectors, w = weights),
    .f = weighted_MAE
  )

  ## record mean squared difference for each MLE
  pp_WMSE_samples[[i]] <- pp_WMSE
  pp_WMAE_samples[[i]] <- pp_WMAE
  if(i %% 10 == 1) print(i)
}
print("done.")

## Save / load results
#print("Saving pp samples.")
#rds_path <- paste0("./Output/data/pp_WMSE_samples/pp_WMSE_samples_",seed_index,".RDS")
#saveRDS(pp_WMSE_samples, file = rds_path)
#rds_path <- paste0("./Output/data/pp_WMAE_samples/pp_WMAE_samples_",seed_index,".RDS")
#saveRDS(pp_WMAE_samples, file = rds_path)

pp_EWMSE <- map_dbl(pp_WMSE_samples, mean)
pp_EWMAE <- map_dbl(pp_WMAE_samples, mean)

#thresholds_vec[which.min(pp_EWMSE)]
#thresholds_vec[which.min(pp_EWMAE)]

# # plot metrics at each threshold
# pdf_path <- paste0("output/qq_flat_threshold_sim_multiple/pp_metric_plots/pp_metrics_",seed_index,".pdf")
# pdf(file = pdf_path, width = 7, height = 5)
#
# plot(x = thresholds_vec, y = pp_EWMSE, main = paste0('selected ',thresholds_vec[which.min(pp_EWMSE)], ' true ',v), xlab = 'threshold')
# abline(v = v, col = 'grey70')
# abline(v = thresholds_vec[which.min(pp_EWMSE)])
#
# plot(x = thresholds_vec, y = pp_EWMAE, main = paste0('selected ',thresholds_vec[which.min(pp_EWMAE)], ' true ',v), xlab = 'threshold')
# abline(v = v, col = 'grey70')
# abline(v = thresholds_vec[which.min(pp_EWMAE)])
# dev.off()

#_______________________________________________________________________________
###
## Calculate QQ distance for each sampled mle ----
###
print("Calculating QQ distances.")
n_eval_pts <- 501
n_rep <- 10 ## latent samples per mle sample
qq_MSE_samples <- vector(mode = "list", length = length(thresholds))
qq_MAE_samples <- vector(mode = "list", length = length(thresholds))
CI_probs <- c(0.025,0.975)
print("starting...")
for(i in seq_along(thresholds)){
  ## Sample true magnitudes above that threshold
  latent_samples <- purrr::pmap(
    .l = list(sig_u = rep(mle_samples[[i]][1, ], n_rep),
              xi = rep(mle_samples[[i]][2, ], n_rep)),
    .f = sample_latent_magnitudes,
    x = cats_filtered[[i]]$m,
    v = cats_filtered[[i]]$threshold,
    u = 0,
    to_nearest = to_nearest
  )

  ## Transform to have common Exp(1) distribution
  latent_samples_exp <- purrr::pmap(
    .f = standardise_gpd_sample,
    .l = list(
      y = latent_samples,
      sig_u = rep(mle_samples[[i]][1,], n_rep),
      xi = rep(mle_samples[[i]][2,], n_rep)),
    v = cats_filtered[[i]]$threshold,
    u = 0
  )

  # get evaluation points equally spaced on probability scale
  eval_probs <- ppoints(n_eval_pts)
  model_quants <- qexp(p = eval_probs, rate = 1)
  sample_quants <- map(.x = latent_samples_exp, .f = quantile, probs = eval_probs, na.rm = TRUE)

  ## Calculate discrepancy from Exp(1) qq-plot
  qq_dist <- function(sample_q,  model_q){ifelse(test = is.finite(sample_q), yes = model_q - sample_q, no =  NA)}
  qq_deviation_vectors <- map(.x = sample_quants, .f = qq_dist, model_q = model_quants)

  # Combine to weighted sum for each sample
  qq_MSE <- map_dbl(.x = qq_deviation_vectors, .f = function(d){mean(d^2, na.rm = TRUE)})
  qq_MAE <- map_dbl(.x = qq_deviation_vectors, .f = function(d){mean(abs(d), na.rm = TRUE)})


  ## record mean squared difference for each MLE
  qq_MSE_samples[[i]] <- qq_MSE
  qq_MAE_samples[[i]] <- qq_MAE
  if(i %%10 ==1) print(i)
}
print("done.")

## Save / load results
#print("Saving qq distances.")
#rds_path <- paste0("./Output/data/qq_MSE_samples/qq_MSE_samples_",seed_index,".RDS")
#saveRDS(qq_MSE_samples, file = rds_path)

#rds_path <- paste0("./Output/data/qq_MAE_samples/qq_MAE_samples_",seed_index,".RDS")
#saveRDS(qq_MAE_samples, file = rds_path)
#print("done.")

## Find expected metric values and the threshold that minimises them
qq_EMSE <- map_dbl(qq_MSE_samples, mean)
qq_EMAE <- map_dbl(qq_MAE_samples, mean)
#thresholds_vec[which.min(qq_EMSE)]
#thresholds_vec[which.min(qq_EMAE)]

# # Plot the expected metrics against treshold values
# pdf_path <- paste0("output/qq_flat_threshold_sim_multiple/qq_metric_plots/qq_metrics_",seed_index,".pdf")
#
# pdf(file = pdf_path, width = 7, height = 5)
# #par(mfrow = c(1,2))
# plot(x = thresholds_vec, y= qq_EMSE, main = paste0('selected ',thresholds_vec[which.min(qq_EMSE)], ' true ',v))
# abline(v = thresholds_vec[which.min(qq_EMSE)])
# abline(v = v, col = 'grey70')
#
# plot(x = thresholds_vec, y= qq_EMAE, main = paste0('selected ',thresholds_vec[which.min(qq_EMAE)], ' true ',v))
# abline(v = thresholds_vec[which.min(qq_EMAE)])
# abline(v = v, col = 'grey70')
# dev.off()

#_______________________________________________________________________________
###
## Plot all metrics together ----
###
# print("Plotting all metrics together.")
# pdf_path <- paste0("./Output/plots/all_metric_plots/all_metrics_",seed_index,".pdf")
#
# pdf(file = pdf_path, width = 9, height = 9)
# par(mfrow = c(2,2))
#
# plot(x = thresholds_vec, y = pp_EWMSE, main = paste0('selected ',thresholds_vec[which.min(pp_EWMSE)], ' true ',v))
# abline(v = v, col = 'grey70')
# abline(v = thresholds_vec[which.min(pp_EWMSE)])
#
# plot(x = thresholds_vec, y = pp_EWMAE, main = paste0('selected ',thresholds_vec[which.min(pp_EWMAE)], ' true ',v))
# abline(v = v, col = 'grey70')
# abline(v = thresholds_vec[which.min(pp_EWMAE)])
#
# plot(x = thresholds_vec, y= qq_EMSE, main = paste0('selected ',thresholds_vec[which.min(qq_EMSE)], ' true ',v))
# abline(v = thresholds_vec[which.min(qq_EMSE)])
# abline(v = v, col = 'grey70')
#
# plot(x = thresholds_vec, y= qq_EMAE, main = paste0('selected ',thresholds_vec[which.min(qq_EMAE)], ' true ',v))
# abline(v = thresholds_vec[which.min(qq_EMAE)])
# abline(v = v, col = 'grey70')
# dev.off()
# print("done.")

# Save expected error metrics at each threshold
saveRDS(object = pp_EWMSE, file = paste0("./Output/data/pp_EWMSE/pp_EWMSE_", jobid,".RDS"))
saveRDS(object = pp_EWMAE, file = paste0("./Output/data/pp_EWMAE/pp_EWMAE_", jobid,".RDS"))
saveRDS(object = qq_EMSE, file = paste0("./Output/data/qq_EMSE/qq_EMSE_", jobid,".RDS"))
saveRDS(object = qq_EMAE, file = paste0("./Output/data/qq_EMAE/qq_EMAE_", jobid,".RDS"))


print("/n")
print(paste("completed catalogue", seed_index, 'of', length(seeds)))
print('END.')

## Plot a histogram
#pdf(paste0("./Output/plots/",jobid,".pdf"))
#hist(x, main = paste0(jobid))
#dev.off()

## Save
#save(x,file=paste0("./Output/data/",jobid,".RData"))
