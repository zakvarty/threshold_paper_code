##  batch_step_thresh.r
###
## Adaptation of qq_step_threshold_sim.r for running on many simulated catalogues
## in parallel on storm.
###
## Author: Zak Varty
##
## _Note to future self:_ save the average metric values to speed up download if
## using more than 10 jobs.
START_TIME <- Sys.time()
#_______________________________________________________________________________
### Get the jobid from the command line argument
jobid  <- as.numeric(commandArgs(trailing=TRUE))
print(paste0("This is batch_step_thresh.r, where jobid = ", jobid, "."))

#_______________________________________________________________________________
###
## Source required scripts ----
###

source("./src/_gpd.R")
source("./src/round_to_nearest.R")
source("./src/mle_gpd_rd_varu.R")
source("./src/sample_mles_gpd_rd_stepped.R")
#source("src/sample_mles_gpd_rd_varu.R")
#source("src/sample_mles_gpd_rd_varu_II.R")
source("./src/bootstrap_stability_plot.R")
source('src/bootstrapping.R')

source("./src/empirical_cdf.R")
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
seeds <- sample(1:1e7, size = 1e3, replace = FALSE)
print("seeds 1:10 are:")
seeds[1:10]

#_______________________________________________________________________________
###
## Use jobid to retrieve and set the seed for this run ----
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
sig <- 0.5
xi  <- -0.03
u <- 0.4
n   <- 1000
to_nearest <- 0.1
v_values <- c(0.42,0.42)
v_true <- rep(x = v_values, each = n / 2)

# Simulate full catalogue (reproducibly)
m_raw <- rgpd(n = n, scale = sig, shape = xi, mu = u)
cat_full <- data.frame(
  m_raw = m_raw,
  u = rep(u, n),
  v = v_true,
  index = seq_along(m_raw))

# reduce to rounded, censored catalogue
cat <- cat_full %>%
  filter(m_raw >= v) %>%
  mutate(m_rounded = round_to_nearest(m_raw, to_nearest),
         index = seq_along(m_rounded)) %>%
  select(m = m_rounded, u, v, index)

# rds_path <- paste0("./Output/data/cats/cat_",seed_index,".RDS")
# saveRDS(object = cat, file = rds_path)
# print("catalogue saved")

# Get number of events on each side of change
n_1 <- floor(sum(cat$v == v_values[1])/2)
n_2 <- ceiling(sum(cat$v == v_values[2])/2)
change_point <- n_1 + 0.5

#_______________________________________________________________________________
###
## Set stepped modelling thresholds to try ----
###

# Set thresholds (assuming decreasing)
thresholds_vec <- c(seq(0.3, 1.0, by = 0.025))
n_thresh <- length(thresholds_vec)

threshold_pairs <-  data.frame(v_1 = NULL, v_2 = NULL)
for ( i in 1:n_thresh){
  for ( j in 1:i){
    new_row <- data.frame(v_1 = thresholds_vec[i], v_2 = thresholds_vec[j])
    threshold_pairs <- rbind(threshold_pairs, new_row)
  }
}

# ## plot theshold pair values and true threshold
# plot(threshold_pairs, pch = 16, col = 'grey70')
# points(x = max(v_true), y = min(v_true), pch = 16, col = 1, cex = 0.8)

thresholds <- map2(.x = threshold_pairs$v_1,
                   .y = threshold_pairs$v_2,
                   .f = function(x,y){c(rep(x,n_1), rep(y, n_2))})

# pdf_path <- paste0("/Output/plots/cat_plots/cat_plot_",seed_index,".pdf")
# pdf("./Output/plots/step_cat.pdf", width = 7, height = 5)
# # Plot catalogue and thresholds
# plot(
#   x = cat$index,
#   y = cat$m,
#   ylim = c(0,4),
#   bty = 'n',
#   type = 'n',
#   main = "simulated catalogue for stepped threshold selection",
#   ylab = "magnitude",
#   xlab = "index"
# )
# abline(h = thresholds_vec, col = 'gray', lwd = 0.4)
# abline(v = n_1 + 0.5, col = 'gray')
# points(cat$index, cat$m, pch = 16, cex = 0.8)
# #abline(h = mu, col = 'lightgray')
# lines(x = c(1, n_1 + 0.5, n_1 + 0.5, n_1 + n_2), y = v_true %>% range %>% rev %>% rep(each = 2))
# dev.off()

#_______________________________________________________________________________
###
## Subset data based on each threshold ----
###

filter_catalogue <- function(cat, threshold, to_nearest){
  cat %>%
    mutate( threshold = threshold) %>%
    filter(m  > (threshold - to_nearest/2 + 1e-8))
}

cats_filtered <- map(
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
# rds_path <- paste0("./Output/data/mles/mles_",seed_index,".RDS")
# saveRDS(object = mles, file = rds_path)
# print("mles saved")

#_______________________________________________________________________________
###
## Calculate expected number of exceedances at each threshold ----
###
print("calculating expected exceedances")

w_vectors <- pmap(
  .f = get_w_vector,
  .l = list(
    x = mags_filtered,
    v = thresholds_filtered,
    gpd_mle = mles
  ),
  to_nearest = 0.1,
  u = 0)

# #Plot expeceted vs observed exceedances of stepped threshold: diagnostic only
# plot(
#   y = unlist(expected_exceedances),
#   x = map_dbl(.x = mags_filtered,.f =  length),
#   xlab = 'catalogue length',
#   ylab = 'expected number of exceedances',
#   main = 'at all thresholds, expected values <= observed count')
# abline(coef = c(0,1))


#_______________________________________________________________________________
###
## Sample MLE values for each pair of thresholds ----
###
print("sampling MLE values at each pair of thresholds")
mle_samples <- vector(mode = "list", length = length(mles))
print('starting...')
# this will take a while...
for(i in seq_along(thresholds)){
  if(i %% 10 == 1){
    message("iteration: ", i, "/", length(thresholds))
  }

  mle_samples[[i]] <- sample_mles_gpd_rd_stepped(
    mle = mles[[i]],
    u = 0,
    step_lengths = c(n_1, n_2),
    step_values = as.numeric(threshold_pairs[i,]),
    w = w_vectors[[i]],
    to_nearest = 0.1,
    n_sim = 750,
    verbose = TRUE)
}
print('done.')

# print('saving MLE samples')
# rds_path <- paste0("./Output/data/mle_samples/mle_samples_",seed_index,".RDS")
# saveRDS(object = mle_samples, file = rds_path)
# print('done.')

#_______________________________________________________________________________
###
## Calculate PP distances for each sampled mle ----
###
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

print("calculating PP-distances...")
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
  print(i)
}
print("done.")

# print("saving PP-distances")
# rds_path <- paste0("./Output/data/pp_WMSE_samples/pp_WMSE_samples_",seed_index,".RDS")
# saveRDS(object = pp_WMSE_samples, file = rds_path)
# rds_path <- paste0("./Output/data/pp_WMAE_samples/pp_WMAE_samples_",seed_index,".RDS")
# saveRDS(object = pp_WMAE_samples, file = rds_path)
# print("done.")

#_______________________________________________________________________________
###
## Calculate QQ distances for each sampled mle ----
###
n_eval_pts <- 501
n_rep <- 10 ## latent samples per mle sample
qq_MSE_samples <- vector(mode = "list", length = length(thresholds))
qq_MAE_samples <- vector(mode = "list", length = length(thresholds))
CI_probs <- c(0.025,0.975)

print("calculating  qq-distances...")
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


  ## record metric for each MLE
  qq_MSE_samples[[i]] <- qq_MSE
  qq_MAE_samples[[i]] <- qq_MAE
  print(i)
}
print('done.')

# Save distance measures
# print("saving QQ-distances")
# rds_path <- paste0("./Output/data/qq_MSE_samples/qq_MSE_samples_",seed_index,".RDS")
# saveRDS(object = qq_MSE_samples, file = rds_path)
# rds_path <- paste0("./Output/data/qq_MAE_samples/qq_MAE_samples_",seed_index,".RDS")
# saveRDS(object = qq_MAE_samples, file = rds_path)
# print("done.")

#_______________________________________________________________________________
###
## Calculate expected PP and QQ metrics for each threshold pair ----
###
print("Calculating expected metric values")
metrics <- data.frame(
  v_1 = threshold_pairs$v_1,
  v_2 = threshold_pairs$v_2,
  pp_EWMSE = map_dbl(pp_WMSE_samples, mean),
  pp_EWMAE = map_dbl(pp_WMAE_samples, mean),
  qq_EMSE  = map_dbl(qq_MSE_samples, mean),
  qq_EMAE  = map_dbl(qq_MAE_samples, mean))
print("done.")

print("Saving expected metric values")
rds_path <- paste0("./Output/data/metrics/metrics_",seed_index,".RDS")
saveRDS(metrics, rds_path)

print(paste0("true thresholds are:" , v_values, "."))
print(paste0("pp_EWMSE chose:" , metrics[which.min(metrics$pp_EWMSE),1:2], "."))
print(paste0("pp_EWMAE chose:" , metrics[which.min(metrics$pp_EWMAE),1:2], "."))
print(paste0("qq_EMSE  chose:" , metrics[which.min(metrics$qq_EMSE),1:2] , "."))
print(paste0("qq_EMAE  chose:" , metrics[which.min(metrics$qq_EMAE),1:2] , "."))


#_______________________________________________________________________________
###
## Threshold selection plots: PP and QQ distances ----
###
#
# library(ggplot2)
# pdf_path <- paste0("./Output/plots/all_metric_plots/all_metrics_",seed_index,".pdf")
# pdf(file = pdf_path, width = 7, height = 5)
# # PP squared error
# ggplot(data = NULL) +
#   geom_point(data = metrics, aes(x = v_1, y = v_2), color = 'grey70') +
#   # geom_point(data = metrics, size = 2.5, aes(x = v_1, y = v_2) , color = 'grey80')
#   geom_vline(xintercept = max(v_true), color = 'grey50') +
#   geom_hline(yintercept = min(v_true), color = 'grey50') +
#   geom_point(data = metrics %>%
#                mutate(is_minimum = as.factor (pp_EWMSE==min(pp_EWMSE))) %>%
#                filter(pp_EWMSE < 0.03),
#              size = 3,
#              aes(x = v_1, y = v_2, color =   pp_EWMSE , shape = is_minimum)) +
#   scale_color_continuous(type = "viridis") +
#   theme_minimal()
# # PP absolute error
# ggplot(data = NULL) +
#   geom_point(data = metrics, aes(x = v_1, y = v_2), color = 'grey70') +
#   # geom_point(data = metrics, size = 2.5, aes(x = v_1, y = v_2) , color = 'grey80')
#   geom_vline(xintercept = max(v_true), color = 'grey50') +
#   geom_hline(yintercept = min(v_true), color = 'grey50') +
#   geom_point(data = metrics %>%
#                mutate(is_minimum = as.factor (pp_EWMAE==min(pp_EWMAE))) %>%
#                filter(pp_EWMAE < 1.2),
#              size = 3,
#              aes(x = v_1, y = v_2, color = pp_EWMAE, shape = is_minimum)) +
#   scale_color_continuous(type = "viridis") +
#   theme_minimal()
# # QQ squared error
# ggplot(data = NULL) +
#   geom_point(data = metrics, aes(x = v_1, y = v_2), color = 'grey70') +
#   # geom_point(data = metrics, size = 2.5, aes(x = v_1, y = v_2) , color = 'grey80')
#   geom_vline(xintercept = max(v_true), color = 'grey50') +
#   geom_hline(yintercept = min(v_true), color = 'grey50') +
#   geom_point(data = metrics %>%
#                mutate(is_minimum = as.factor (qq_EMSE==min(qq_EMSE))) %>%
#                filter(qq_EMSE < 0.03),
#              size = 3,
#              aes(x = v_1, y = v_2, color =   qq_EMSE , shape = is_minimum)) +
#   scale_color_continuous(type = "viridis") +
#   theme_minimal()
# # QQ absolute error
# ggplot(data = NULL) +
#   geom_point(data = metrics, aes(x = v_1, y = v_2), color = 'grey70') +
#   # geom_point(data = metrics, size = 2.5, aes(x = v_1, y = v_2) , color = 'grey80')
#   geom_vline(xintercept = max(v_true), color = 'grey50') +
#   geom_hline(yintercept = min(v_true), color = 'grey50') +
#   geom_point(data = metrics %>%
#                mutate(is_minimum = as.factor (qq_EMAE==min(qq_EMAE))) %>%
#                filter(qq_EMAE < 0.09),
#              size = 3,
#              aes(x = v_1, y = v_2, color = qq_EMAE , shape = is_minimum)) +
#   scale_color_continuous(type = "viridis") +
#   theme_minimal()
# dev.off()

END_TIME <- Sys.time()
RUN_TIME <- END_TIME - START_TIME

print(paste("completed catalogue", seed_index))
print(paste("time taken:", RUN_TIME))
print('END.')

## Plot a histogram
#pdf(paste0("./Output/plots/",jobid,".pdf"))
#hist(x, main = paste0(jobid))
#dev.off()

## Save
#save(x,file=paste0("./Output/data/",jobid,".RData"))
