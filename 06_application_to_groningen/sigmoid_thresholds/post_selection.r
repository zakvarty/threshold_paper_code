#_______________________________________________________________________________
print("This is sigmioid_Groningen/post_selection.r")


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

source("./src/fitting_exp_submodel/llh_exp_rd_varu.R")
source("./src/fitting_exp_submodel/mle_exp_rd_varu.R")
source("./src/fitting_exp_submodel/sample_mles_exp_rd_stepped.R")
source("./src/fitting_exp_submodel/get_mle_and_bootstrap_exp.R")

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

#_______________________________________________________________________________
###
## Set seed
###
set.seed(1234)
seeds <- sample(x = 1:1e7,size = 1e3,replace = FALSE)

#_______________________________________________________________________________
###
## Load Groningen catalogue ----
###
rds_path <- "./Output/data/gron_cat.RDS"
gron_cat <- readRDS(rds_path) %>% mutate(m = mag) %>% select(-mag)

#_______________________________________________________________________________
###
## Load top threshold thresholds ----
###

#rds_path <- "./Output/data/restricted_selection/top_thresholds.RDS"

selected_threshold_parameters <-  list(v_1 = 1.15, v_2 = 0.76, mu = 746, sigma = 1)

selected_threshold <- get_sigmoid_values(
  tau = gron_cat$index,
  mu = selected_threshold_parameters$mu,
  sigma = selected_threshold_parameters$sigma,
  v_1 = selected_threshold_parameters$v_1,
  v_2 = selected_threshold_parameters$v_2)

selected_threshold_flat <- rep(1.07, NROW(gron_cat))

conservative_threshold <- rep(1.45, NROW(gron_cat))
gron_cat <- gron_cat %>%
  mutate(v_cons = conservative_threshold,
         v_select = selected_threshold,
         v_flat = selected_threshold_flat)

################################################################################
################################################################################
###########         Part 1: Compare metric values              #################
################################################################################
################################################################################
n_repeats <- 100
metric_repeats <- data.frame(
  rep = 1:n_repeats,
  cons = rep(NA, n_repeats),
  flat = rep(NA, n_repeats),
  sigmiod = rep(NA, n_repeats)
)

FUN_threshold <- function(threshold, plot = FALSE){

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

  return(dq1_value)
}

cat <- gron_cat
# for(i in 1:n_repeats){
#   metric_repeats$cons[i] <- FUN_threshold(threshold = conservative_threshold)
#   print(c(i,1))
#   metric_repeats$flat[i] <- FUN_threshold(threshold = selected_threshold_flat)
#   print(c(i,2))
#   metric_repeats$sigmiod[i] <- FUN_threshold(threshold = selected_threshold)
#   print(c(i,3))
# }
# print("done")

#saveRDS(metric_repeats, "./Output/data/post_selection_threshold_repeats.RDS")
metric_repeats <- readRDS("./Output/data/post_selection_threshold_repeats.RDS")

colMeans(metric_repeats)
summarise_all(metric_repeats, quantile, probs = c(0.025,0.975))

pdf("Output/plots/post_selection/metric_repeats_post_selection.pdf", width = 7, height = 5)
plot(x = as.factor(rep(c("cons", "flat", "sigmoid"), each = n_repeats)),
     y = c(metric_repeats$cons, metric_repeats$flat, metric_repeats$sigmiod),
     log = 'y',
     xlab = 'threshold type',
     ylab = "d(q,1)")
dev.off()
qqnorm(metric_repeats$flat)
qqline(metric_repeats$flat)

################################################################################
################################################################################
###########         Part 2: Fit rounded GPD models             #################
################################################################################
################################################################################

## Get maximum likelihood estimates and bootstrap CI using each threshold
selected_mle <- get_mle_and_bootstrap(
  cat = gron_cat,
  threshold = selected_threshold,
  to_nearest = 0.1,
  u = -1,
  n_bootstrap_mles = 10000,
  verbose = TRUE)

selected_flat_mle <- get_mle_and_bootstrap(
  cat = gron_cat,
  threshold = selected_threshold_flat,
  to_nearest = 0.1,
  u = -1,
  n_bootstrap_mles = 10000,
  verbose = TRUE)

conservative_mle <- get_mle_and_bootstrap(
  cat = gron_cat,
  threshold = conservative_threshold,
  to_nearest = 0.1,
  u = -1,
  n_bootstrap_mles = 10000,
  verbose = TRUE)

selected_mle_exp <- get_mle_and_bootstrap_exp(
  cat = gron_cat,
  threshold = selected_threshold,
  to_nearest = 0.1,
  u = -1,
  n_bootstrap_mles = 10000,
  verbose = TRUE
)

selected_flat_mle_exp <- get_mle_and_bootstrap_exp(
  cat = gron_cat,
  threshold = selected_threshold_flat,
  to_nearest = 0.1,
  u = -1,
  n_bootstrap_mles = 10000,
  verbose = TRUE)

conservative_mle_exp <- get_mle_and_bootstrap_exp(
  cat = gron_cat,
  threshold = conservative_threshold,
  to_nearest = 0.1,
  u = -1,
  n_bootstrap_mles = 10000,
  verbose = TRUE
)

# # Save MLE objects
# saveRDS(selected_mle, "./Output/data/post_selection/selected_mle.RDS")
#saveRDS(selected_flat_mle,"./Output/data/post_selection/selected_flat_mle.RDS")
# saveRDS(conservative_mle,"./Output/data/post_selection/conservative_mle.RDS")
# saveRDS(selected_mle_exp,"./Output/data/post_selection/selected_mle_exp.RDS")
# saveRDS(selected_flat_mle_exp,"./Output/data/post_selection/selected_flat_mle_exp.RDS")
# saveRDS(conservative_mle_exp,"./Output/data/post_selection/conservative_mle_exp.RDS")

# Load MLE objects
selected_mle <- readRDS("./Output/data/post_selection/selected_mle.RDS")
selected_flat_mle <- readRDS("./Output/data/post_selection/selected_flat_mle.RDS")
conservative_mle <- readRDS("./Output/data/post_selection/conservative_mle.RDS")
selected_mle_exp <- readRDS("./Output/data/post_selection/selected_mle_exp.RDS")
selected_flat_mle_exp <- readRDS("./Output/data/post_selection/selected_flat_mle_exp.RDS")
conservative_mle_exp <- readRDS("./Output/data/post_selection/conservative_mle_exp.RDS")

# Xi bootstrap confidence intervals under each model
quantile(conservative_mle$boot_mles[,2], prob = c(0.025,0.975))
quantile(selected_flat_mle$boot_mles[,2], prob = c(0.025,0.975))
quantile(selected_mle$boot_mles[,2], prob = c(0.025,0.975))
quantile(selected_mle_exp$boot_mles[,2], prob = c(0.025,0.975))

# bootstrap confidence interval widths (xi)
diff(quantile(conservative_mle$boot_mles[,2], prob = c(0.025,0.975)))
diff(quantile(selected_flat_mle$boot_mles[,2], prob = c(0.025,0.975)))
diff(quantile(selected_mle$boot_mles[,2], prob = c(0.025,0.975)))

###
## Compare exp and GPD above sigmoid and conservative
###
# LR test sigmoid + GPD vs sigmoid + EXP
LR_test_statistic <- 2 * (selected_mle$llh - selected_mle_exp$llh)
LR_test_p_value <- pchisq(LR_test_statistic, df = 1,lower.tail = FALSE)
## 0.064
mean(selected_mle$boot_mles[,2] > 0) ## 0.015

# LR test flat + GPD vs flat + Exp
LR_test_statistic <- 2 * (selected_flat_mle$llh - selected_flat_mle_exp$llh)
LR_test_p_value <- pchisq(LR_test_statistic, df = 1,lower.tail = FALSE)
## 0.046
mean(selected_flat_mle$boot_mles[,2] > 0) ## 0.0055

# LR test conservative + GPD vs conservative + Exp
LR_test_statistic <- 2 * (conservative_mle$llh - conservative_mle_exp$llh)
LR_test_p_value <- pchisq(LR_test_statistic, df = 1,lower.tail = FALSE)
## 0.786
mean(conservative_mle$boot_mles[,2] > 0) ## 0.3305

###
## Get the expected sample size using each threshold.
###
E_count_selected <- sum(
  get_w_vector(
    x = gron_cat %>% filter(m > (v_select - 0.1 / 2 + 1e-8)) %>% pull(m),
    v = gron_cat %>% filter(m > (v_select - 0.1 / 2 + 1e-8)) %>% pull(v_select),
    to_nearest = 0.1,
    gpd_mle = selected_mle$mle,
    u = -1))

E_count_selected_flat <- sum(
  get_w_vector(
    x = gron_cat %>% filter(m > (v_flat - 0.1 / 2 + 1e-8)) %>% pull(m),
    v = gron_cat %>% filter(m > (v_flat - 0.1 / 2 + 1e-8)) %>% pull(v_flat),
    to_nearest = 0.1,
    gpd_mle = selected_mle$mle,
    u = -1))

E_count_conservative <- sum(
  get_w_vector(
    x = gron_cat %>% filter(m > (v_cons - 0.1 / 2 + 1e-8)) %>% pull(m),
    v = gron_cat %>% filter(m > (v_cons - 0.1 / 2 + 1e-8)) %>% pull(v_cons),
    to_nearest = 0.1,
    gpd_mle = selected_mle$mle,
    u = -1))

E_count_selected_flat - E_count_conservative
E_count_selected - E_count_conservative
###
## Sigmoid is worth how many extra obs over each flat threshold?
###
equivalent_extra_obs <- function(boot_mles_1, boot_mles_2, n_1){
  sd_1 <- apply(boot_mles_1, 2, sd)
  sd_2 <- apply(boot_mles_2, 2, sd)
  n_extra <- n_1 * (sd_1/sd_2)^2 - n_1
  return(n_extra)
}

# flat vs cons (363)
equivalent_extra_obs(
  boot_mles_1 =  conservative_mle$boot_mles,
  boot_mles_2 = selected_flat_mle$boot_mles,
  n_1 = E_count_conservative)
# sigmoid vs cons (509)
equivalent_extra_obs(
  boot_mles_1 =  conservative_mle$boot_mles,
  boot_mles_2 = selected_mle$boot_mles,
  n_1 = E_count_conservative)
# sigmoid vs flat (53)
equivalent_extra_obs(
  boot_mles_1 = selected_flat_mle$boot_mles,
  boot_mles_2 =  selected_mle$boot_mles,
  n_1 = E_count_conservative)

###
## Plot bootstrap CIs
###
pdf("./Output/plots/post_selection/parameter_comparison_cons_flat_sigmoid.pdf", width = 7, height = 5)
#png("./Output/plots/post_selection/parameter_comparison_cons_flat_sigmoid_1.png", width = 7, height = 5, units = "in")
par(mar= c(5.1,5.1,2.1,2.1), las = 0)
plot(y = conservative_mle$boot_mles[,2],
     x = conservative_mle$boot_mles[,1] + (2.45)*conservative_mle$boot_mles[,2],
     xlab = "",
     ylab = "",
     xlim = c(0.3,0.7),
     cex.axis = 2,
     cex.lab = 2,
     las = 3)
    mtext(text=expression(hat(xi)), side=2, line=3, las=1,cex = 2)
    mtext(text=expression(hat(sigma)[1.45]), side=1, line=4, las=1,cex = 2)
points(y = selected_flat_mle$boot_mles[,2],
       x = selected_flat_mle$boot_mles[,1] + (2.45)*selected_flat_mle$boot_mles[,2],
       col = 2)
points(y = selected_mle$boot_mles[,2],
       x = selected_mle$boot_mles[,1] + (2.45)*selected_mle$boot_mles[,2],
       col = 4)

plot(y = conservative_mle$boot_mles[,2],
     x = conservative_mle$boot_mles[,1] ,
     xlab = expression(hat(sigma)[-1]),
     ylab = expression(hat(xi)),
     xlim = c(0.1,1.1),
     cex.axis = 2,
     cex.lab = 1.5)
points(y = selected_flat_mle$boot_mles[,2],
       x = selected_flat_mle$boot_mles[,1],
       col = 2)
points(y = selected_mle$boot_mles[,2],
       x = selected_mle$boot_mles[,1],
       col = 4)

dev.off()# legend('topright',
#        legend = c("conservative", "sigmoid"),
#        title = "Threshold",
#        pch = c(1,1),
#        col = c(1,2))


source("./src/plot_return_levels.r")
pdf("./Output/plots/post_selection/return_level_comparison_cons_flat_sigmoid.pdf", width = 7, height = 5)
par(mar= c(5.1,5.1,2.1,2.1))
plot_return_levels(conservative_mle, u = -1, v = 1.45, col = 1, log = 'x', lwd = 2, cex.axis = 2, cex.lab = 2)
plot_return_levels(selected_flat_mle, u = -1, v = 1.45, col = 2, add = TRUE, lwd = 2)
plot_return_levels(selected_mle, u = -1, v = 1.45, col = 4, add = TRUE, lwd = 2)
#return_periods = 1/(0.1^(seq(0,2.5,length.out = 51)))
#points(x = return_periods,
#       y = quantile(x = cat$m[cat$m > 1.45], probs = 1 - 1/ return_periods),
#       col = 4, cex = 1.5, pch = 16)
dev.off()

# Upper endpoint proxy
cons_ueps <- plot_return_levels(conservative_mle, u = -1, v = 1.45,v_return_periods = c(1000,10000) , alpha = 0.2)
selected_uep <- plot_return_levels(selected_mle, u = -1, v = 1.45,v_return_periods = c(1000,10000), add = TRUE, col = 2, alpha = 0.2)

calculate_EUP <- function(mle_obj, u, alpha = 0.05){
  out <- list()
  mle <- u - mle_obj$mle[1] / mle_obj$mle[2]
  mle[mle_obj$mle[2] >= 0] <- NA

  boot_mles <- u - mle_obj$boot_mles[,1] / mle_obj$boot_mles[,2]
  boot_mles[mle_obj$boot_mles[,2] >= 0] <- Inf
  boot_CI <- quantile(boot_mles,probs = c(alpha/2, 1- alpha/2))


  out <- list(mle = mle, boot_CI =  boot_CI,alpha = alpha, boot_mles = boot_mles)
  return(out)
}

conservative_uep <- calculate_EUP(conservative_mle,u = -1 , alpha = 0.05)
selected_uep <- calculate_EUP(selected_mle, u = -1 , alpha = 0.05)

conservative_uep$mle
conservative_uep$boot_CI

selected_uep$mle
selected_uep$boot_CI
