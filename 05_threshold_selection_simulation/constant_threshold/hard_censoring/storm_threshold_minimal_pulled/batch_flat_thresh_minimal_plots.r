##  batch_flat_thresh_minimal_plots.r
###
## Create plots from qq_flat_threshold_sim_multiple using storm data output
## (downloading pdfs would really slow things down.)
###
## Author: Zak Varty

#_______________________________________________________________________________
###
## Source required scripts ----
###
getwd()


#_______________________________________________________________________________
###
## Load required packages ----
###
library("purrr")
library("dplyr")
library("ggplot2")


#_______________________________________________________________________________
###
## Set true parameters ----
###

# True parameters
sig_u <- 0.55
xi  <- - 0.1
u   <- 0
v <- 0.32
sig_0 <- sig_u + (0 - u) * xi
n   <- 1500
to_nearest <- 0.1

thresholds_vec <- c(seq(0.0, 1.0, by = 0.025))

n_runs <- 500
selected <- tibble(
  jobid = rep(NA_real_, n_runs),
  true_threshold = rep(NA_real_, n_runs),
  qq_EMSE = rep(NA_real_, n_runs),
  qq_EMAE = rep(NA_real_, n_runs),
  pp_EWMSE = rep(NA_real_, n_runs),
  pp_EWMAE = rep(NA_real_, n_runs)
)

#_______________________________________________________________________________
###
## For each simulated catalogue... ----
###

###
## Get selected threshold by each metric
###
for(run_number in 1:n_runs){
   ###
   ## Get selected threshold by each metric
   ###
   ## Load the associated catalogue
   rds_path <- paste0('./Output/data/pp_EWMAE/pp_EWMAE_',run_number,'.RDS')
   pp_EWMAE <- readRDS(rds_path)

   rds_path <- paste0('./Output/data/pp_EWMSE/pp_EWMSE_',run_number,'.RDS')
   pp_EWMSE <- readRDS(rds_path)

   rds_path <- paste0('./Output/data/qq_EMAE/qq_EMAE_',run_number,'.RDS')
   qq_EMAE <- readRDS(rds_path)

   rds_path <- paste0('./Output/data/qq_EMSE/qq_EMSE_',run_number,'.RDS')
   qq_EMSE <- readRDS(rds_path)

   selected_this_run <- tibble(
      jobid = run_number,
      true_threshold = v,
      qq_EMSE = thresholds_vec[which.min(qq_EMSE)],
      qq_EMAE = thresholds_vec[which.min(qq_EMAE)],
      pp_EWMSE = thresholds_vec[which.min(pp_EWMSE)],
      pp_EWMAE = thresholds_vec[which.min(pp_EWMAE)]
   )

  selected[run_number,] <- selected_this_run[1,]
  print(run_number)
}


plot(selected$jobid, selected$qq_EMSE, xlab = 'simulated catalogue ID', ylab = 'qq_EMSE selection', ylim = c(0,1))
lines(selected$jobid, selected$true_threshold)
plot(selected$jobid, selected$qq_EMAE, xlab = 'simulated catalogue ID', ylab = 'qq_EMAE selection',ylim = c(0,1))
lines(selected$jobid, selected$true_threshold)
plot(selected$jobid, selected$pp_EWMSE, xlab = 'simulated catalogue ID', ylab = 'pp_EWMSE selection', ylim = c(0,1))
lines(selected$jobid, selected$true_threshold)
plot(selected$jobid, selected$pp_EWMAE, xlab = 'simulated catalogue ID', ylab = 'pp_EWMAE selection', ylim = c(0,1))
lines(selected$jobid, selected$true_threshold)


selection_probs <- tibble(
   threshold = thresholds_vec,
   qq_EMSE = map_dbl(.x = thresholds_vec, .f = function(thresh){mean(selected$qq_EMSE == thresh)}),
   qq_EMAE = map_dbl(.x = thresholds_vec, .f = function(thresh){mean(selected$qq_EMAE == thresh)}),
   pp_EWMSE = map_dbl(.x = thresholds_vec, .f = function(thresh){mean(selected$pp_EWMSE == thresh)}),
   pp_EWMAE = map_dbl(.x = thresholds_vec, .f = function(thresh){mean(selected$pp_EWMAE == thresh)})
)

###
## Summarise threshold selection over 500 simulated catalgoues
###
pdf(file = "./Output/plots/selected_thresholds.pdf", width = 7, height = 5)
plot(
   x = rep(selection_probs$threshold, 4),
   y = as.numeric(unlist(selection_probs[,2:5])),
   type = 'n',
   ylab = 'selection proportion',
   xlab = 'threshold')
lines(x = thresholds_vec, y = selection_probs$pp_EWMSE, col = 1)
lines(x = thresholds_vec, y = selection_probs$pp_EWMAE, col = 2)
lines(x = thresholds_vec, y = selection_probs$qq_EMSE, col = 3)
lines(x = thresholds_vec, y = selection_probs$qq_EMAE, col = 4)
abline(v = v, lty =2)

legend('topright',
       legend = c("pp_EWMSE", "pp_EWMAE","qq_EMSE","qq_EMAE"),
       lty = 1,
       col = 1:4,
       bty = 'n')
dev.off()

## Get RMSE of threshold selections
RMSE_threshold_selection <- function(selection_probs){
   sqrt(mean( rep((thresholds_vec - v)^2, times = selection_probs * n_runs)))
}

RMSE_dq1 <- RMSE_threshold_selection(selection_probs$qq_EMAE)
RMSE_dq2 <- RMSE_threshold_selection(selection_probs$qq_EMSE)
RMSE_dp1 <- RMSE_threshold_selection(selection_probs$pp_EWMAE)
RMSE_dp2 <- RMSE_threshold_selection(selection_probs$pp_EWMSE)


pdf(file = './Output/plots/selected_threshold_summary.pdf', width = 7, height = 5)
par(mfrow = c(2,2))
plot(
   x = selection_probs$threshold,
   y = selection_probs$qq_EMAE,
   type = 'h',
   ylab = 'selection proportion',
   xlab = 'threshold',
   ylim = c(0,1),
   lwd = 2,
   main = paste0("d(q,1), RMSE = ", round(RMSE_dq1, 2)))
abline(v=v , col = 2, lty = 2)
plot(
   x = selection_probs$threshold,
   y = selection_probs$qq_EMSE,
   type = 'h',
   ylab = 'selection proportion',
   xlab = 'threshold',
   ylim = c(0,1),
   lwd = 2,
   main = paste0("d(q,2), RMSE = ", round(RMSE_dq2, 2)))
abline(v=v , col = 2, lty = 2)
plot(
   x = selection_probs$threshold,
   y = selection_probs$pp_EWMAE,
   type = 'h',
   ylab = 'selection proportion',
   xlab = 'threshold',
   ylim = c(0,1),
   lwd = 2,
   main = paste0("d(p,1), RMSE = ", round(RMSE_dp1, 2)))
abline(v = v , col = 2, lty = 2)
plot(
   x = selection_probs$threshold,
   y = selection_probs$pp_EWMSE,
   type = 'h',
   ylab = 'selection proportion',
   xlab = 'threshold',
   ylim = c(0,1),
   lwd = 2,
   main = paste0("d(p,2), RMSE = ", round(RMSE_dp2, 2)))
abline(v = v,col = 2, lty = 2)
dev.off()


## Same plots but on individual
pdf(file = './Output/plots/selected_threshold_summary_individual.pdf', width = 4, height = 3)
par(mar = c(4.1,4.1,3.1,1.1))
plot(
   x = selection_probs$threshold,
   y = selection_probs$qq_EMAE,
   type = 'h',
   ylab = 'selection proportion',
   xlab = 'threshold',
   ylim = c(0,1),
   lwd = 2,
   main = paste0("d(q,1), RMSE = ", round(RMSE_dq1, 2)))
abline(v=v , col = 2, lty = 2, lwd = 2)
plot(
   x = selection_probs$threshold,
   y = selection_probs$qq_EMSE,
   type = 'h',
   ylab = 'selection proportion',
   xlab = 'threshold',
   ylim = c(0,1),
   lwd = 2,
   main = paste0("d(q,2), RMSE = ", round(RMSE_dq2, 2)))
abline(v=v , col = 2, lty = 2, lwd = 2)
plot(
   x = selection_probs$threshold,
   y = selection_probs$pp_EWMAE,
   type = 'h',
   ylab = 'selection proportion',
   xlab = 'threshold',
   ylim = c(0,1),
   lwd = 2,
   main = paste0("d(p,1), RMSE = ", round(RMSE_dp1, 2)))
abline(v = v , col = 2, lty = 2, lwd = 2)
plot(
   x = selection_probs$threshold,
   y = selection_probs$pp_EWMSE,
   type = 'h',
   ylab = 'selection proportion',
   xlab = 'threshold',
   ylim = c(0,1),
   lwd = 2,
   main = paste0("d(p,2), RMSE = ", round(RMSE_dp2, 2)))
abline(v = v,col = 2, lty = 2, lwd = 2)
dev.off()



