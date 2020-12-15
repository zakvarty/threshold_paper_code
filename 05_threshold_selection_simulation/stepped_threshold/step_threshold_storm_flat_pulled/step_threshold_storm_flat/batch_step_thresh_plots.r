##  batch_step_thresh.r
###
## Adaptation of qq_step_threshold_sim.r for running on many simulated catalogues
## in parallel on storm.
## This script uses storm output to construct plots of metric values.
###
## Author: Zak Varty
##

#_______________________________________________________________________________
###
## Load required packages ----
###

library("purrr")
library("dplyr")
library("ggplot2")


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
v_values <- c(0.83,0.42)
v_true <- rep(x = v_values, each = n / 2)


## For each seed_index
seeds <- c(1)

for (seed_index in seeds) {
  # Load catalgoue
  rds_path <- paste0("./Output/data/cats/cat_",seed_index,".RDS")
  cat <- readRDS(rds_path)

  # Get number of events on each side of change
  n_1 <- sum(cat$v == v_values[1])
  n_2 <- sum(cat$v == v_values[2])
  change_point <- n_1 + 0.5

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

  thresholds <- map2(.x = threshold_pairs$v_1,
                     .y = threshold_pairs$v_2,
                     .f = function(x,y){c(rep(x,n_1), rep(y, n_2))})

  ###
  ## plot theshold pair values and true threshold
  ###
  # plot(threshold_pairs, pch = 16, col = 'grey70')
  # points(x = max(v_true), y = min(v_true), pch = 16, col = 1, cex = 0.8)

  ###
  # Plot catalogue and thresholds
  ###
  pdf_path <- paste0("./Output/plots/cat_plots/cat_plot_",seed_index,".pdf")
  pdf(pdf_path, width = 7, height = 5)
  plot(
    x = cat$index,
    y = cat$m,
    ylim = c(0,4),
    bty = 'n',
    type = 'n',
    main = "simulated catalogue for stepped threshold selection",
    ylab = "magnitude",
    xlab = "index"
  )
  abline(h = thresholds_vec, col = 'gray', lwd = 0.4)
  abline(v = n_1 + 0.5, col = 'gray')
  points(cat$index, cat$m, pch = 16, cex = 0.8)
  #abline(h = mu, col = 'lightgray')
  lines(x = c(1, n_1 + 0.5, n_1 + 0.5, n_1 + n_2), y = v_true %>% range %>% rev %>% rep(each = 2))
  dev.off()

  ###
  # Load expected metric values or calculate from samples
  ###

  ## if saved:
  #rds_path <- paste0("./Output/data/metrics/metrics_",seed_index,".RDS")
  #metrics <- readRDS(rds_path)

  ## otherwise:

  rds_path <- paste0("./Output/data/pp_WMSE_samples/pp_WMSE_samples_",seed_index,".RDS")
  pp_WMSE_samples <- readRDS(rds_path)

  rds_path <- paste0("./Output/data/pp_WMAE_samples/pp_WMAE_samples_",seed_index,".RDS")
  pp_WMAE_samples <- readRDS(rds_path)

  rds_path <- paste0("./Output/data/qq_MSE_samples/qq_MSE_samples_",seed_index,".RDS")
  qq_MSE_samples <- readRDS(rds_path)

  rds_path <- paste0("./Output/data/qq_MAE_samples/qq_MAE_samples_",seed_index,".RDS")
  qq_MAE_samples <- readRDS(rds_path)

  metrics <- data.frame(
    v_1 = threshold_pairs$v_1,
    v_2 = threshold_pairs$v_2,
    pp_EWMSE = map_dbl(pp_WMSE_samples, mean),
    pp_EWMAE = map_dbl(pp_WMAE_samples, mean),
    qq_EMSE  = map_dbl(qq_MSE_samples, mean),
    qq_EMAE  = map_dbl(qq_MAE_samples, mean))

  print("Saving expected metric values")
  rds_path <- paste0("./Output/data/metrics/metrics_",seed_index,".RDS")
  saveRDS(metrics, rds_path)

  # Print to console the selected threshold for that catalogue
  print(paste0("--- Catalogue ",seed_index,"---"))
  print(paste0("true thresholds are:" , v_values, "."))
  print(paste0("pp_EWMSE chose:" , metrics[which.min(metrics$pp_EWMSE),1:2], "."))
  print(paste0("pp_EWMAE chose:" , metrics[which.min(metrics$pp_EWMAE),1:2], "."))
  print(paste0("qq_EMSE  chose:" , metrics[which.min(metrics$qq_EMSE),1:2] , "."))
  print(paste0("qq_EMAE  chose:" , metrics[which.min(metrics$qq_EMAE),1:2] , "."))

  #_______________________________________________________________________________
  ###
  ## Threshold selection plots: PP and QQ distances ----
  ###

  pdf_path <- paste0("./Output/plots/all_metric_plots/all_metrics_",seed_index,".pdf")
  pdf(file = pdf_path, width = 7, height = 5)
  # PP squared error
  ggplot(data = NULL) +
    geom_point(data = metrics, aes(x = v_1, y = v_2), color = 'grey70') +
    # geom_point(data = metrics, size = 2.5, aes(x = v_1, y = v_2) , color = 'grey80')
    geom_vline(xintercept = max(v_true), color = 'grey50') +
    geom_hline(yintercept = min(v_true), color = 'grey50') +
    geom_point(data = metrics %>%
                 mutate(is_minimum = as.factor (pp_EWMSE==min(pp_EWMSE))) %>%
                 filter(pp_EWMSE < 0.03),
               size = 3,
               aes(x = v_1, y = v_2, color =   pp_EWMSE , shape = is_minimum)) +
    scale_color_continuous(type = "viridis") +
    theme_minimal()
  # PP absolute error
  ggplot(data = NULL) +
    geom_point(data = metrics, aes(x = v_1, y = v_2), color = 'grey70') +
    # geom_point(data = metrics, size = 2.5, aes(x = v_1, y = v_2) , color = 'grey80')
    geom_vline(xintercept = max(v_true), color = 'grey50') +
    geom_hline(yintercept = min(v_true), color = 'grey50') +
    geom_point(data = metrics %>%
                 mutate(is_minimum = as.factor (pp_EWMAE==min(pp_EWMAE))) %>%
                 filter(pp_EWMAE < 1.2),
               size = 3,
               aes(x = v_1, y = v_2, color = pp_EWMAE, shape = is_minimum)) +
    scale_color_continuous(type = "viridis") +
    theme_minimal()
  # QQ squared error
  ggplot(data = NULL) +
    geom_point(data = metrics, aes(x = v_1, y = v_2), color = 'grey70') +
    # geom_point(data = metrics, size = 2.5, aes(x = v_1, y = v_2) , color = 'grey80')
    geom_vline(xintercept = max(v_true), color = 'grey50') +
    geom_hline(yintercept = min(v_true), color = 'grey50') +
    geom_point(data = metrics %>%
                 mutate(is_minimum = as.factor (qq_EMSE==min(qq_EMSE))) %>%
                 filter(qq_EMSE < 0.03),
               size = 3,
               aes(x = v_1, y = v_2, color =   qq_EMSE , shape = is_minimum)) +
    scale_color_continuous(type = "viridis") +
    theme_minimal()
  # QQ absolute error
  ggplot(data = NULL) +
    geom_point(data = metrics, aes(x = v_1, y = v_2), color = 'grey70') +
    # geom_point(data = metrics, size = 2.5, aes(x = v_1, y = v_2) , color = 'grey80')
    geom_vline(xintercept = max(v_true), color = 'grey50') +
    geom_hline(yintercept = min(v_true), color = 'grey50') +
    geom_point(data = metrics %>%
                 mutate(is_minimum = as.factor (qq_EMAE==min(qq_EMAE))) %>%
                 filter(qq_EMAE < 0.09),
               size = 3,
               aes(x = v_1, y = v_2, color = qq_EMAE , shape = is_minimum)) +
    scale_color_continuous(type = "viridis") +
    theme_minimal()
  dev.off()

}






