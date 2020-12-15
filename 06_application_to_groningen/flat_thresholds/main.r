##  flat_thresholds/main.r
##  Author: Zak Varty
##
##  Selecting a flat threshold for the Groningen catalogue by grid search
##  and comparing to the flat, conservative KNMI threshold of 1.45ML.

print("This is flat_thresholds/main.r")

###########         Part 0: Set up                             #################

## 0.1: Source required scripts ------------------------------------------------
source("../../00_src/_gpd.R")
source("../../00_src/round_to_nearest.R")
source("../../00_src/llh_gpd_rd_varu.R")
source("../../00_src/mle_gpd_rd_varu.R")
source("../../00_src/sample_mles_gpd_rd_stepped.R")
source('../../00_src/bootstrapping.R')

source("../../00_src/thresholding/sample_latent_magnitudes.R")
source("../../00_src/thresholding/standardise_gpd_sample.R")

source("../../00_src/get_mle_and_bootstrap.R")
source("../../00_src/get_dq1.R")
source("../../00_src/get_sigmoid_values.r")

source("../../00_src/fitting_exp_submodel/get_mle_and_bootstrap_exp.R")
source("../../00_src/fitting_exp_submodel/sample_mles_exp_rd_stepped.R")
source("../../00_src/fitting_exp_submodel/mle_exp_rd_varu.R")
source("../../00_src/fitting_exp_submodel/llh_exp_rd_varu.R")

source("../../00_src/plot_return_levels.r")
print("functions sourced.")

## 0.2: Load required packages -------------------------------------------------

library("purrr")
library("dplyr")
library("ggplot2")
library("furrr")
library("ismev")
#library("ParBayesianOptimization")

## 0.3: Set seed ---------------------------------------------------------------

set.seed(1234)
seeds <- sample(x = 1:1e7,size = 1e3,replace = FALSE)

## 0.4: Load Groningnen catalogue ----------------------------------------------

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


#############       Part 1: Fitting above 1.45                 #################

## 1.1: Fit abve 1.45ML using rounded GPD --------------------------------------
conservative_mle <- get_mle_and_bootstrap(
  cat = cat,
  threshold = rep(1.45, NROW(cat)),
  to_nearest = 0.1,
  u =-1,
  n_bootstrap_mles = 10000,
  verbose = TRUE)

# bootstrap mles
plot(conservative_mle$boot_mles)
# CI on sigma
quantile(conservative_mle$boot_mles[,1], prob = c(0.025,0.975))
# CI on xi
quantile(conservative_mle$boot_mles[,2], prob = c(0.025,0.975))


## 1.3:  Fit above 1.45ML using rounded exponential ----------------------------
conservative_exp_mle <- get_mle_and_bootstrap_exp(
  cat = cat,
  threshold = rep(1.45, NROW(cat)),
  to_nearest = 0.1,
  u = -1,
  n_bootstrap_mles = 10000,
  verbose = TRUE
)

# 1.4: Plot bootstrap mles from rounded exponential and rounded GPD models -----
pdf_path <- "./Output/plots/conservative_groningen_fitting/groningen_conservative_boot_mles.pdf"
pdf(pdf_path, width = 7, height = 5)
    par(mar = c(5.1,5.1,2.1,1.1))
    plot(x = conservative_mle$boot_mles[,1] + 2.45 * conservative_mle$boot_mles[,2],
         y = conservative_mle$boot_mles[,2],
         ylab = expression(hat(xi)),
         xlab = expression(hat(sigma)[1.45]),
         bty = 'n',
         cex.axis = 1.5,
         cex.lab = 1.5,
         #main = 'bootstrap mle values',
         xlim = c(0.3,0.6))
    points(x = conservative_exp_mle$boot_mles[,1] + 2.45 * conservative_exp_mle$boot_mles[,2],
           y = conservative_exp_mle$boot_mles[,2],
           col = 2)
    # legend('topright',
    #        legend = c("GPD","Exp"),
    #        pch = 1,
    #        col = c(1,2))
dev.off()


## 1.5: Likelihood ratio test for which model is best --------------------------
## p-value of LR test
pchisq(q = 2 * (conservative_mle$llh - conservative_exp_mle$llh), df = 1)


## 1.6: Return-level plots for each model --------------------------------------
pdf_path <-"./Output/plots/conservative_groningen_fitting/return_levels_conservative.pdf"
pdf(pdf_path,width = 7, height = 5)
    plot_return_levels(
      mle_and_bootstraps = conservative_mle,
      u = -1,
      v = 1.45,
      alpha = 0.05, log = 'x', lwd = 2,
      cex.axis = 1.4,
      cex.lab = 1.4)
    plot_return_levels(
      mle_and_bootstraps = conservative_exp_mle,
      u = -1,
      v = 1.45,
      alpha = 0.05, log = 'x',add = TRUE, col = 2, lwd = 2)
dev.off()


## 1.7: Modified QQ and PP plots above 1.45 ------------------------------------
source("src/modified_pp_plot.r")
source("src/modified_qq_plot.r")
pdf_path <- "./Output/plots/conservative_groningen_fitting/groningen_conservative_mod_PP_QQ.pdf"
pdf(pdf_path, width = 6, height = 6)
  modified_pp_plot(
    m = gron_cat$mag,
    threshold = rep(1.45,NROW(gron_cat)),
    mle_obj = conservative_mle,
    u = -1,
    to_nearest = 0.1,
    n_rep = 10,
    n_eval_pts = 101,
    cex.axis = 1.5,
    cex.lab = 1.5)

  modified_qq_plot(
    m = gron_cat$mag,
    threshold = rep(1.45,NROW(gron_cat)),
    mle_obj = conservative_mle,
    u = -1,
    to_nearest = 0.1,
    n_rep = 10,
    n_eval_pts =  101,
    cex.axis = 1.5,
    cex.lab = 1.5)
dev.off()




###########     Part 2: Fitting unrounded  above 1.45      #####################

## 2.1: Fit above 1.45ML using standard GPD ------------------------------------
ismev <- gpd.fit(xdat = gron_cat$mag, threshold = 1.5, npy = 1)
ismev$mles <- cbind(rnorm(n = 10000,mean = ismev$mle[1], ismev$se[1]),
                    rnorm(n = 10000,mean = ismev$mle[2], ismev$se[2]))

## 2.2: Compare parameter estimates with rounded GPD ---------------------------
hist(ismev$mles[,1] - conservative_mle$boot_mles[,1])
mean((ismev$mles[,2] - conservative_mle$boot_mles[,2])>0)
# difference is not significant relative to parameter uncertainty.



###########    Part 3: Selecting a flat threshold          #####################

## NB: clear workspace and run section 0

# Function to evaluate d(q,1) for a flat threshold
FUN_1 <- function(threshold_1, visual = FALSE){
  START_TIME <- Sys.time()
  threshold <- rep(threshold_1, NROW(cat))

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

  if(visual){
    plot(cat, main = paste0("d(q,1) =  ",round(dq1_value,5)))
    lines(x = cat$index, y = threshold, col = 2)
  }

  END_TIME <- Sys.time()
  print(END_TIME - START_TIME)
  return(dq1_value)
}
# example use:
# FUN_1(1.4,visual = TRUE)

# 3.1: Grid search for best flat threshold. ------------------------------------

# NB: If you want to run in series, not parallell, comment out plan() and
#     switch furrr:: to purrr:: equivalent.

plan(multisession, workers = 4)
potential_thresholds <- seq(0.3,1.8, by = 0.01)
dq1_values <- furrr:::future_map_dbl(.x = potential_thresholds, .f = FUN_1,.progress = TRUE)

metric_values <- data.frame(threshold = potential_thresholds, dq1_value = dq1_values)
#saveRDS(object = metric_values, file = './Output/data/flat_threshold_selection/metric_values.RDS')

## 3.2: Plot metric evaluations ------------------------------------------------

metric_values <- readRDS('./Output/data/flat_threshold_selection/metric_values.RDS')

pdf("./Output/plots/flat_threshold_selection/groningen_flat_metric_values.pdf", width = 7, height = 5)
    par(mar = c(5.1,5.1,2.1,1.1))
    plot(x = metric_values$threshold,
         y = (metric_values$dq1_value),
         cex = 0.8,
         pch = 16,
         log = 'y',
         xlab = "threshold value, v",
         ylab = "d(q,1)",
         cex.axis = 2,
         cex.lab = 2)
    abline(v =  seq(0.15,2.05,by = .1), lty = 2)
dev.off()

# Get best threshold and metric value for that threshold
min(dq1_values)
potential_thresholds[which.min(dq1_values)]

# 3.3: Get MLE and bootstrap for conservative & selected thresholds ------------

conservative_mle <- get_mle_and_bootstrap(
  cat = cat,
  threshold = rep(1.45, NROW(cat)),
  to_nearest = 0.1,
  u =-1,
  n_bootstrap_mles = 10000,
  verbose = TRUE)

chosen_mle <- get_mle_and_bootstrap(
  cat = cat,
  threshold = rep(1.07, NROW(cat)),
  to_nearest = 0.1,
  u =-1,
  n_bootstrap_mles = 10000,
  verbose = TRUE)

## 3.4: Get expected exceedance count for each threshold -----------------------
ENA_cons <- cat %>%
  filter(m >= 1.45 - 0.1 / 2 + 1e-8) %>%
  pull(m) %>%
  get_w_vector(v = 1.45,
               to_nearest = 0.1,
               gpd_mle = conservative_mle$mle,
               u =  -1) %>%
  sum()

ENA_chosen <- cat %>%
  filter(m >= 1.07 - 0.1 / 2 + 1e-8) %>%
  pull(m) %>%
  get_w_vector(v = 1.07,
               to_nearest = 0.1,
               gpd_mle = chosen_mle$mle,
               u =  -1) %>%
  sum()

## 3.5: Plot bootstrap mles using each threshold -------------------------------
pdf_path <- "./Output/plots/flat_threshold_selection/groningen_selected_boot_mles.pdf"

pdf(pdf_path, width = 7, height = 5)
    par(mar = c(5.1,5.1,2.1,1.1))
    plot(x = conservative_mle$boot_mles[,1] + 2.45 * conservative_mle$boot_mles[,2],
         y = conservative_mle$boot_mles[,2],
         ylab = expression(hat(xi)),
         xlab = expression(hat(sigma)[1.45]),
         bty = 'n',
         cex.axis = 1.5,
         cex.lab = 1.5,
         #main = 'bootstrap mle values',
         xlim = c(0.3,0.6))
    points(x = chosen_mle$boot_mles[,1] + 2.45 * chosen_mle$boot_mles[,2],
           y = chosen_mle$boot_mles[,2],
           col = 2)
    # legend('topright',
    #        legend = c("GPD","Exp"),
    #        pch = 1,
    #        col = c(1,2))
dev.off()

## 3.6: Plot estimated return levels above 1.45ML using each threshold ---------
pdf_path <- "./Output/plots/flat_threshold_selection/groningen_flat_selected_return_levels.pdf"

pdf(pdf_path, width = 7, height = 5)
    par(mar = c(5.1,5.1,2.1,2.1))
    rl_conservative <- plot_return_levels(
      mle_and_bootstraps = conservative_mle,
      u = -1,
      v = 1.45,
      alpha = 0.05, log = 'x', lwd = 2,
      cex.axis = 2,
      cex.lab = 2, col = 1)
    rl_chosen <- plot_return_levels(
      mle_and_bootstraps = chosen_mle,
      u = -1,
      v = 1.45,
      alpha = 0.05, log = 'x',add = TRUE, col = 2, lwd = 2)
    return_periods = 1/(0.1^(seq(0,2.5,length.out = 51)))
    points(x = return_periods,
          y = quantile(x = cat$m[cat$m > 1.45], probs = 1 - 1/ return_periods),
          col = 4, cex = 1.5, pch = 16)
dev.off()


## 3.7: Calculate relative standard deviations ---------------------------------
sds_chosen <- apply(chosen_mle$boot_mles, MARGIN = 2,sd)
sds_conservative <- apply(conservative_mle$boot_mles, MARGIN = 2, sd)
rel_sds <- sds_chosen / sds_conservative


## 3.8: Calculate the effective # extra, indept datapts above 1.45 -------------
Effective_extra_observations <- rel_sds^(-2) *  ENA_cons - ENA_cons

rl_CI_widths_chosen <- rl_chosen$CI_high - rl_chosen$CI_low
rl_CI_widths_conservative <- rl_conservative$CI_high - rl_conservative$CI_low

plot(
  x = rl_chosen$return_levels$return_period,
  y = rl_CI_widths_chosen - rl_CI_widths_conservative)

