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


# Storage for selected thresholds
n_runs <- 10
selected <- tibble(
  jobid = rep(NA_real_, n_runs),
  true_v1 = rep(NA_real_, n_runs),
  true_v2 = rep(NA_real_, n_runs),
  pp_EWMSE_v1 = rep(NA_real_, n_runs),
  pp_EWMSE_v2 = rep(NA_real_, n_runs),
  pp_EWMAE_v1 = rep(NA_real_, n_runs),
  pp_EWMAE_v2 = rep(NA_real_, n_runs),
  qq_EMSE_v1 = rep(NA_real_, n_runs),
  qq_EMSE_v2 = rep(NA_real_, n_runs),
  qq_EMAE_v1 = rep(NA_real_, n_runs),
  qq_EMAE_v2 = rep(NA_real_, n_runs),
)

## For each seed_index
seeds <- 1:n_runs
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
    main = "Stepped thesehold, hard censoring",
    ylab = "magnitude",
    xlab = "index",
    cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
  )
  #abline(h = thresholds_vec, col = 'gray', lwd = 0.4)
  #abline(v = n_1 + 0.5, col = 'gray')
  lines(x = c(1, n_1 + 0.5, n_1 + 0.5, n_1 + n_2),
        y = v_true %>% range %>% rev %>% rep(each = 2),
        col = 2,
        lwd = 3)
  points(cat$index, cat$m, pch = 16, cex = 1)
  #abline(h = mu, col = 'lightgray')
  #lines(x = c(1, n_1 + 0.5, n_1 + 0.5, n_1 + n_2), y = v_true %>% range %>% rev %>% rep(each = 2),
        col = 2, lwd = 3, lty = 2)
  dev.off()

  ###
  # Load expected metric values or calculate from samples
  ###

  ## if saved:
  rds_path <- paste0("./Output/data/metrics/metrics_",seed_index,".RDS")
  metrics <- readRDS(rds_path)

  ## otherwise:

  # rds_path <- paste0("./Output/data/pp_WMSE_samples/pp_WMSE_samples_",seed_index,".RDS")
  # pp_WMSE_samples <- readRDS(rds_path)
  #
  # rds_path <- paste0("./Output/data/pp_WMAE_samples/pp_WMAE_samples_",seed_index,".RDS")
  # pp_WMAE_samples <- readRDS(rds_path)
  #
  # rds_path <- paste0("./Output/data/qq_MSE_samples/qq_MSE_samples_",seed_index,".RDS")
  # qq_MSE_samples <- readRDS(rds_path)
  #
  # rds_path <- paste0("./Output/data/qq_MAE_samples/qq_MAE_samples_",seed_index,".RDS")
  # qq_MAE_samples <- readRDS(rds_path)
  #
  # metrics <- data.frame(
  #   v_1 = threshold_pairs$v_1,
  #   v_2 = threshold_pairs$v_2,
  #   pp_EWMSE = map_dbl(pp_WMSE_samples, mean),
  #   pp_EWMAE = map_dbl(pp_WMAE_samples, mean),
  #   qq_EMSE  = map_dbl(qq_MSE_samples, mean),
  #   qq_EMAE  = map_dbl(qq_MAE_samples, mean))
  #
  # print("Saving expected metric values")
  # rds_path <- paste0("./Output/data/metrics/metrics_",seed_index,".RDS")
  # saveRDS(metrics, rds_path)

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

  ##
  # Create plot objects
  ##
  # PP squared error
  pp_squared_error_plot <- ggplot(data = NULL) +
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
  pp_absolute_error_plot <- ggplot(data = NULL) +
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
  qq_squared_error_plot <- ggplot(data = NULL) +
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
  qq_absolute_error_plot <- ggplot(data = NULL) +
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

  # All metrics in one file
  pdf_path <- paste0("./Output/plots/all_metric_plots/all_metrics_",seed_index,".pdf")
  pdf(file = pdf_path, width = 7, height = 5)
  print(pp_squared_error_plot)
  print(pp_absolute_error_plot)
  print(qq_squared_error_plot)
  print(qq_absolute_error_plot)
  dev.off()

  # pp_EWMSE
  pdf_path <- paste0("./Output/plots/pp_EWMSE_plots/pp_EWMSE_",seed_index,".pdf")
  pdf(file = pdf_path, width = 7, height = 5)
    print(pp_squared_error_plot)
  dev.off()

  # pp_EWMAE
  pdf_path <- paste0("./Output/plots/pp_EWMAE_plots/pp_EWMAE_",seed_index,".pdf")
  pdf(file = pdf_path, width = 7, height = 5)
    print(pp_absolute_error_plot)
  dev.off()

  # qq_EMSE
  pdf_path <- paste0("./Output/plots/qq_EMSE_plots/qq_EMSE_",seed_index,".pdf")
  pdf(file = pdf_path, width = 7, height = 5)
    print(qq_squared_error_plot)
  dev.off()

  # qq_EMAE
  pdf_path <- paste0("./Output/plots/qq_EMAE_plots/qq_EMAE_",seed_index,".pdf")
  pdf(file = pdf_path, width = 7, height = 5)
    print(qq_absolute_error_plot)
  dev.off()

  selected_this_run <- tibble(
    jobid = seed_index,
    true_v1 = v_values[1],
    true_v2 = v_values[2],
    pp_EWMSE_v1 = as.numeric(metrics[which.min(metrics$pp_EWMSE),1]),
    pp_EWMSE_v2 = as.numeric(metrics[which.min(metrics$pp_EWMSE),2]),
    pp_EWMAE_v1 = as.numeric(metrics[which.min(metrics$pp_EWMAE),1]),
    pp_EWMAE_v2 = as.numeric(metrics[which.min(metrics$pp_EWMAE),2]),
    qq_EMSE_v1 = as.numeric(metrics[which.min(metrics$qq_EMSE),1]),
    qq_EMSE_v2 = as.numeric(metrics[which.min(metrics$qq_EMSE),2]),
    qq_EMAE_v1 = as.numeric(metrics[which.min(metrics$qq_EMAE),1]),
    qq_EMAE_v2 = as.numeric(metrics[which.min(metrics$qq_EMAE),2]),
  )

  selected[seed_index,] <- selected_this_run[1,]

}

###
## Save threshold values chosen at each run
###
rds_path <- "./Output/data/selected.RDS"
saveRDS(selected, rds_path)

###
## Save selected thresholds for each run
###
opar <- par()
## v_1 selections
pdf_path <- "./Output/plots/selected_v1.pdf"
pdf(pdf_path, width = 7 , height = 5)
layout(rbind(1,2), heights=c(7,1))
plot(x = rep(selected$jobid, 4) + rep(0.1*c(-1, -0.5, 0.5, 1), each =10),
     y = c(
       selected$pp_EWMSE_v1,
       selected$pp_EWMAE_v1,
       selected$qq_EMSE_v1,
       selected$qq_EMAE_v1
     ),
     col = rep(1:4, each = 10),
     ylab = 'selected v_1',
     xlab = 'catalogue')
abline(h = selected$true_v1[1])
par(mar=c(0, 0, 0, 0))
plot.new()
legend('center',
       legend = c("pp_EWMSE","pp_EWMAE","qq_EMSE","qq_EMAE"),
       col = 1:4,
       pch = 1,
       ncol= 4,
       bty = "n")
par(mar = opar$mar)
dev.off()

## v_2 selections
pdf_path <- "./Output/plots/selected_v2.pdf"
pdf(pdf_path, width = 7 , height = 5)
layout(rbind(1,2), heights=c(7,1))
plot(x = rep(selected$jobid, 4) + rep(0.1*c(-1, -0.5, 0.5, 1), each =10),
     y = c(
       selected$pp_EWMSE_v2,
       selected$pp_EWMAE_v2,
       selected$qq_EMSE_v2,
       selected$qq_EMAE_v2
     ),
     col = rep(1:4, each = 10),
     ylab = 'selected v_2',
     xlab = 'catalogue')
abline(h = selected$true_v2[1])
par(mar=c(0, 0, 0, 0))
plot.new()
legend('center',
       legend = c("pp_EWMSE","pp_EWMAE","qq_EMSE","qq_EMAE"),
       col = 1:4,
       pch = 1,
       ncol= 4,
       bty = "n")
par(mar = opar$mar)
dev.off()


###
# plot the total error over all runs by each metric on each threshold
###

## v_1 errors
errors <- selected
summed_abs_errors <- errors %>%
  transmute(
   # jobid = jobid,
    pp_EWMSE_v1 = abs(pp_EWMSE_v1 - true_v1),
    pp_EWMAE_v1 = abs(pp_EWMAE_v1 - true_v1),
    qq_EMSE_v1  = abs(qq_EMSE_v1  - true_v1),
    qq_EMAE_v1  = abs(qq_EMAE_v1  - true_v1),
    pp_EWMSE_v2 = abs(pp_EWMSE_v2 - true_v2),
    pp_EWMAE_v2 = abs(pp_EWMAE_v2 - true_v2),
    qq_EMSE_v2  = abs(qq_EMSE_v2 - true_v2),
    qq_EMAE_v2  = abs(qq_EMAE_v2  - true_v2)
  ) %>%
  summarise_all(mean)

##
# Mean absolute error in threshold choice by each metric
##
optimal_v1 <- thresholds_vec[which.max(thresholds_vec > v_values[1])]
optimal_v2 <- thresholds_vec[which.max(thresholds_vec > v_values[2])]
min_error_v1 <- (optimal_v1 - v_values[1])
min_error_v2 <- (optimal_v2 - v_values[2])

pdf_path <- "./Output/plots/mean_absolute_selection_errors.pdf"
pdf(pdf_path, width = 7 , height = 5)
layout(cbind(1,2), widths = c(6.5,1.5))
plot(x = rep(c(1,2), each = 4),
     y = as.numeric(summed_abs_errors[1,]),
     xlab = 'threshold portion',
     ylab = 'mean absolute error',
     main = 'threshold selection error over 10 replications',
     ylim = c(0, max(as.numeric(summed_abs_errors[1,]))),
     col = rep(1:4,2),
     pch = "-")
points(x = c(1,2), y = c(min_error_v1, min_error_v2), pch = '-', col = 5)
par(mar=c(0, 0, 0, 0))
plot.new()
legend('center',
        title = 'metric',
       legend = c("pp_EWMSE","pp_EWMAE","qq_EMSE","qq_EMAE","optimal"),
       col = c(1,2,3,4,5),
       pch = c(rep(95,5)),
       ncol= 1,
       bty = "n")
par(mar = opar$mar)
dev.off()

