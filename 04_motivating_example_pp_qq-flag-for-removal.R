##
# Making transformed QQ and PP plots for motivating example
#
# Note 2025-11-03: Not sure if this is needed anymore. Flag for removal.
##

SEED = 1234
thresholds_to_try <- c(0.5, 1.15, 1.85)
## 1: Source necessary code ----------------------------------------------------
#source('../00_src/_gpd.R')
#source('../00_src/round_to_nearest.R')

#source('../00_src/llh_gpd_rd_varu.R')
#source('../00_src/mle_gpd_rd_varu.R')

#source('../00_src/sample_mles_gpd_rd_stepped.R')
#source("../00_src/thresholding/sample_latent_magnitudes.R")
#source("../00_src/thresholding/standardise_gpd_sample.R")

#source('00_src/bootstrapping.R')

library(dplyr)
library(purrr)
library(threshold)

# 2: Simulate example catalogue ------------------------------------------------

## 2.1: Set underlying parameters ----------------------------------------------
cat_size = 1000
sig_0 = 0.55
xi = -0.1
u = 1.05
sig_u = sig_0 + xi * u
to_nearest = 0.1
step_values = c(1.65, 1.05)


## 2.2: Create full catalogue --------------------------------------------------

set.seed(SEED)

cat <- tibble(
  index = seq_len(cat_size),
  mag = rgpd_rd(n = cat_size, sig_u, xi, shift = u, to_nearest = to_nearest),
  v_cons = rep(step_values[1], cat_size),
  v_step = rep(step_values, each = cat_size / 2),
  exceeds_cons = (mag >= v_cons),
  exceeds_step = (mag >= v_step),
  exceeds_extra = (mag >= step_values[1]))

# Get the number of events in conservative or stepped catalogue
n_cons <- sum(cat$exceeds_cons)
n_step <- sum(cat$exceeds_step)
n_extra <- n_step - n_cons

# simulate additional excesses of the conservative threshold to investigate
# whether improvements are from increased sample size or observing a wider range
# of magnitude values.

scale_extra <- sig_u + xi * (step_values[1] - u)

cat_extra <- tibble(
  index = seq(from = cat_size + 1, to = cat_size + n_extra),
  mag = rgpd_rd(n_extra, scale_extra, xi, step_values[1], to_nearest),
  v_cons = rep(step_values[1], n_extra),
  v_step = rep(step_values[2], n_extra),
  exceeds_cons = (mag >= v_cons),
  exceeds_step = (mag >= v_step),
  exceeds_extra = (mag >= step_values[1])
)

cat_combined <- rbind(cat, cat_extra)
cat_combined$is_extra <- cat_combined$index > cat_size
cat_combined <-  mutate(cat_combined,
  is_extra = index > cat_size,
  plotting_colour = case_when(
    mag > v_cons ~ "grey10",
    mag < v_cons & index <= cat_size/2 ~ "grey60",
    mag < v_cons & index > cat_size/2 ~ "grey85")
)

rm(scale_extra, cat_extra)

## 2.3: Thin by censoring and threshold type  ----------------------------------

## 2.4: Plot catalogues --------------------------------------------------------
pdf_path <- "./00_outputs/04_motivating_example/cat_plots_motivating_example.pdf"
pdf(file = pdf_path, width = 4, height = 4)
par(mar = c(4.1,4.1,2.1,2.1))
plot(
  #main = "conservative",
  x = cat$index,
  y = cat$mag,
  pch = 16,
  col = if_else(cat$mag > cat$v_cons, true = "grey10", false = "grey60") ,
  cex = 0.8,
  bty = 'n',
  ylim = c(1,4),
  xlab = 'event index',
  ylab = 'magnitude')
lines(x = c(0,cat_size), y = rep(step_values[1],2), col = "purple", lwd = 2)

plot(
  #main = "stepped",
  x = cat$index,
  y = cat$mag,
  pch = 16,
  col = if_else(cat$mag > cat$v_step, true = "grey10", false = "grey60") ,
  cex = 0.8,
  bty = 'n',
  ylim = c(1,4),
  xlab = 'event index',
  ylab = 'magnitude')
mid_pt <- (cat_size + 1)/2
lines(x = c(0, mid_pt, mid_pt, cat_size), y = rep(step_values, each = 2), lwd = 2, col = "darkorange")

plot(
  #main = "extended",
  x = cat_combined$index,
  y = cat_combined$mag,
  pch = 16,
  col = cat_combined$plotting_colour ,
  cex = 0.8,
  bty = 'n',
  ylim = c(1,4),
  xlab = 'event index',
  ylab = 'magnitude')
lines(x = c(0,nrow(cat_combined)), y = rep(step_values[1],2), col = "purple", lwd = 2)

plot(
  #main = "trial thresholds",
  x = cat$index,
  y = cat$mag,
  pch = 16,
  col = if_else(cat$mag > cat$v_step, true = "grey10", false = NA) ,
  cex = 0.8,
  bty = 'n',
  ylim = c(0.5,4),
  xlab = 'event index',
  ylab = 'magnitude')
mid_pt <- (cat_size + 1)/2
lines(x = c(0, mid_pt, mid_pt, cat_size), y = rep(step_values, each = 2), lwd = 2, col = "darkorange")
for (i in seq_along(thresholds_to_try)) {
lines(x = c(0,cat_size), y = rep(thresholds_to_try[i],2), col = "purple", lwd = 2)
}
dev.off()


# 3.0 Estimate GPD parameters above a range of thresholds ----------------------
## 3.1 Set flat threshold heights to try ---------------------------------------

mles <- vector(mode = "list", length = length(thresholds_to_try))
boot_mles <- vector(mode = "list", length = length(thresholds_to_try))


for (i in seq_along(thresholds_to_try)) {
  v_value <- thresholds_to_try[i]
  print(paste("starting bootstrap for threshold", i))
  boot_mles <- parametric_bootstrap_mles(x = cat$mag, v = rep(v_value, nrow(cat)))
}
print("done")

#sigma_0_mles <- map_dbl(mles, function(vec){vec[1]})
#xi_mles <- map_dbl(mles, function(vec){vec[2]})

##
## ZV TO HERE 13-11-2025 - next step: add sample_latent_magnitudes to {threshold}
##

## 4.0: Construct modified PP and QQ plots for each threshold ------------------

# Set simulation parameters and bootstrap CI quantiles
PP_QQ <- list(
  n_eval_pts = 500,          ## m = number of probs at which to eval distance?
  n_rep = 100,               ## B = number of bootstrap samples?
  CI_probs = c(0.025, 0.975) ## CIs on plot
)

## Do it for one bootstrap first

y_sims <- purrr::pmap_dfr(
  .f = sample_unrounded_gpd,
  .l = list(scale = boot_mles$bootstrap_mles$scale,
  shape = boot_mles$bootstrap_mles$shape),
  x = boot_mles$config$x,
  shift = boot_mles$config$shift,
  to_nearest = boot_mles$config$to_nearest)

temp <- threshold::unround_gpd(n = 1,
                       x = boot_mles$config$x,
                       scale = boot_mles$bootstrap_mles$scale[1],
                       shape = boot_mles$bootstrap_mles$shape[1],
                       shift = boot_mles$config$shift,
                       to_nearest = boot_mles$config$to_nearest)

y_sim <- vector(mode = "list", length = PP_QQ$n_rep)
z_sim <- vector(mode = "list", length = PP_QQ$n_rep)
u_sim <- vector(mode = "list", length = PP_QQ$n_rep)
g_sim <- vector(mode = "list", length = PP_QQ$n_rep)

for (b in 1:PP_QQ$n_rep) {
  ## Sample latenet magnitudes
  y_sim[[b]] <- threshold::sample_unrounded_gpd(
    x = boot_mles$config$x,
    scale = boot_mles$bootstrap_mles$scale[b],
    shape = boot_mles$bootstrap_mles$shape[b],
    shift = boot_mles$config$shift,
    to_nearest = boot_mles$config$to_nearest)
  ## filter to threshold excesses
  is_excess <- y_sim[[b]] > boot_mles$config$v
  ## Transform to Exp(1) margins
  z_sim[[b]] <- threshold::gpd_to_exp(
    y = y_sim[[b]][is_excess],
    scale = boot_mles$bootstrap_mles$scale[b],
    shape = boot_mles$bootstrap_mles$shape[b],
    shift = boot_mles$config$shift,
    v = boot_mles$config$v[is_excess])
  # Transform to Uniform margins
  u_sim[[b]] <- pexp(q = z_sim[[b]])
  # Transforem to Gaussian margins
  g_sim[[b]] <- qnorm(p = u_sim[[b]])
}

theo_samples_exp <- matrix(data = rexp(n = length(z_sim[[1]]) * length(z_sim)),
                                       nrow = length(z_sim),
                                       ncol = length(z_sim[[1]]))


# Make plots

pdf_path <- "./output/plots/qq_pp_gauss_margins_motivating_example.pdf"
pdf(file = pdf_path, width = 5, height = 5)
    par(mar = c(4.1,4.1,1.1,1.1))
    for(i in seq_along(thresholds)){

      ## Sample true magnitudes above that threshold
      latent_samples <- purrr::pmap(
        .l = list(sig_u = rep(boot_mles[[i]][1, ], n_rep),
                  xi = rep(boot_mles[[i]][2, ], n_rep)),
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
        u = 1.0
      )

#       latent_samples_gauss <- map(
#         .x = latent_samples_exp,
#         .f = function(x){qnorm(pexp(x))}
#       )
#
#       ## Simulate actual N(0,1) vectors of same lengths
#       theo_samples_gauss <- map(
#         .x = latent_samples_exp,
#         .f = function(x){rnorm(n = length(x))})
#
#       ###
#       ## Construct PP plot on Gaussian margins
#       ###
#
#       # vector of x-coords for pp plot
#       p <- ppoints(n_eval_pts)
#       std_norm_quants <- qnorm(p = p)
#
#       get_sample_ecdf <- function(x, q){ecdf(x)(q)}
#
#       observed_ecdf_values <- map(
#         .x = latent_samples_gauss,
#         .f = get_sample_ecdf,
#         q = std_norm_quants)
#       observed_ecdf_values <- do.call(rbind, observed_ecdf_values)
#       observed_ecdf_CI <- apply(X = observed_ecdf_values,
#                                 MARGIN = 2,
#                                 FUN = quantile,
#                                 probs = CI_probs)
#
#       theoretical_ecdf_values <- map(
#         .x = theo_samples_gauss,
#         .f = get_sample_ecdf,
#         q = std_norm_quants
#       )
#       theoretical_ecdf_values <- do.call(rbind, theoretical_ecdf_values)
#       theoretical_ecdf_CI <- apply(X = theoretical_ecdf_values,
#                                    MARGIN = 2,
#                                    FUN = quantile,
#                                    probs = CI_probs)
#
#       above_sim_int <- observed_ecdf_CI[1,] > theoretical_ecdf_CI[2,]
#       below_sim_int <- observed_ecdf_CI[2,] < theoretical_ecdf_CI[1,]
#       point_colours <- rep('black', length(p))
#       point_colours[above_sim_int] <- 'red'
#       point_colours[below_sim_int] <- 'blue'
#
#       plot(
#         x = rep(p,4),
#         y = c(theoretical_ecdf_CI[1,],
#               theoretical_ecdf_CI[2,],
#               observed_ecdf_CI[1,],
#               observed_ecdf_CI[2,]),
#         ylab = 'sample probability',
#         xlab = 'model probability',
#         type = 'n',
#         asp = 1)
#       polygon(
#         x = c(p, rev(p)),
#         y = c(theoretical_ecdf_CI[1,], rev(theoretical_ecdf_CI[2,])),
#         col = "grey60",
#         border = NA
#       )
#       for(j in seq_along(p)){
#         lines(x = rep(p[j],2), y = observed_ecdf_CI[,j], col = point_colours[j])
#       }
#
#       ###
#       ## Construct QQ plot on Gaussian margins
#       ###
#       p <- ppoints(n_eval_pts)
#       std_norm_quants <- qnorm(p)
#
#       observed_quantiles <- map(
#         .x = latent_samples_gauss,
#         .f = quantile,
#         probs = p,
#         na.rm = TRUE
#       )
#       observed_quantiles_CI <- do.call(rbind, observed_quantiles)
#       observed_quantiles_CI <- apply(observed_quantiles_CI, 2, quantile, probs = CI_probs)
#
#       theoretical_quantiles <- map(
#         .x = theo_samples_gauss,
#         .f = quantile,
#         probs = p,
#         na.rm = TRUE
#       )
#       theoretical_quantiles_CI <- do.call(rbind, theoretical_quantiles)
#       theoretical_quantiles_CI <- apply(theoretical_quantiles_CI, 2, quantile, probs = CI_probs)
#
#       above_sim_int <- observed_quantiles_CI[1,] > theoretical_quantiles_CI[2,]
#       below_sim_int <- observed_quantiles_CI[2,] < theoretical_quantiles_CI[1,]
#       point_colours <- rep('black', length(p))
#       point_colours[above_sim_int] <- 'red'
#       point_colours[below_sim_int] <- 'blue'
#
#       plot(
#         x = rep(std_norm_quants,4),
#         y = c(theoretical_quantiles_CI[1,], theoretical_quantiles_CI[2,], observed_quantiles_CI[1,], observed_quantiles_CI[2,]),
#         ylab = 'sample quantile',
#         xlab = 'model quantiles',
#         type = 'n',
#         asp = 1)
#       polygon(
#         x = c(std_norm_quants, rev(std_norm_quants)),
#         y = c(theoretical_quantiles_CI[1,], rev(theoretical_quantiles_CI[2,])),
#         col = "grey60",
#         border = NA
#       )
#       for(j in seq_along(p)){
#         lines(
#           x = rep(std_norm_quants[j],2),
#           y = observed_quantiles_CI[,j],
#           col = point_colours[j])
#       }
#     }
# dev.off()

## 5.0: Construct modified PP and QQ plots of exponential margins --------------
pdf_path <- "./output/plots/qq_pp_exp_margins_motivating_example.pdf"
pdf(file = pdf_path, width = 5, height = 5)
    par(mar = c(4.1,4.1,1.1,1.1))

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
        u = 1.0
      )

      # latent_samples_gauss <- map(
      #   .x = latent_samples_exp,
      #   .f = function(x){qnorm(pexp(x))}
      # )

      ## Simulate actual N(0,1) vectors of same lengths
      #theo_samples_gauss <- map(
      #  .x = latent_samples_exp,
      #  .f = function(x){rnorm(n = length(x))})

      ## Simulate actual Exp(1) vectors of same lengths
      theo_samples_exp <- map(
        .x = latent_samples_exp,
        .f = function(x){rexp(n = length(x))})


      ###
      ## Construct PP plot on exponential margins
      ###

      # vector of x-coords for pp plot
      p <- ppoints(n_eval_pts)
      #std_norm_quants <- qnorm(p = p)
      std_exp_quants <- qexp(p = p)

      get_sample_ecdf <- function(x, q){ecdf(x)(q)}

      observed_ecdf_values <- map(
        #.x = latent_samples_gauss,
        .x = latent_samples_exp,
        .f = get_sample_ecdf,
        #q = std_norm_quants
        q = std_exp_quants)
      observed_ecdf_values <- do.call(rbind, observed_ecdf_values)
      observed_ecdf_CI <- apply(observed_ecdf_values, 2, quantile, probs = CI_probs)

      theoretical_ecdf_values <- map(
        #.x = theo_samples_gauss,
        .x = theo_samples_exp,
        .f = get_sample_ecdf,
        #q = std_norm_quants
        q = std_exp_quants
      )
      theoretical_ecdf_values <- do.call(rbind, theoretical_ecdf_values)
      theoretical_ecdf_CI <- apply(X = theoretical_ecdf_values,
                                   MARGIN = 2,
                                   FUN = quantile,
                                   probs = CI_probs)

      above_sim_int <- observed_ecdf_CI[1,] > theoretical_ecdf_CI[2,]
      below_sim_int <- observed_ecdf_CI[2,] < theoretical_ecdf_CI[1,]
      point_colours <- rep('black', length(p))
      point_colours[above_sim_int] <- 'red'
      point_colours[below_sim_int] <- 'blue'

      plot(
        x = rep(p,4),
        y = c(theoretical_ecdf_CI[1,],
              theoretical_ecdf_CI[2,],
              observed_ecdf_CI[1,],
              observed_ecdf_CI[2,]),
        ylab = 'sample probability',
        xlab = 'model probability',
        type = 'n',
        asp = 1)
      polygon(
        x = c(p, rev(p)),
        y = c(theoretical_ecdf_CI[1,], rev(theoretical_ecdf_CI[2,])),
        col = "grey60",
        border = NA
      )
      for(j in seq_along(p)){
        lines(x = rep(p[j],2), y = observed_ecdf_CI[,j], col = point_colours[j])
      }

      ###
      ## Construct QQ plot on exponential margins
      ###
      p <- ppoints(n_eval_pts)
      #std_norm_quants <- qnorm(p)
      std_exp_quants <- qexp(p)

      observed_quantiles <- map(
        #.x = latent_samples_gauss,
        .x = latent_samples_exp,
        .f = quantile,
        probs = p,
        na.rm = TRUE
      )
      observed_quantiles_CI <- do.call(rbind, observed_quantiles)
      observed_quantiles_CI <- apply(X = observed_quantiles_CI,
                                     MARGIN = 2,
                                     FUN = quantile,
                                     probs = CI_probs)

      theoretical_quantiles <- map(
        #.x = theo_samples_gauss,
        .x = theo_samples_exp,
        .f = quantile,
        probs = p,
        na.rm = TRUE
      )
      theoretical_quantiles_CI <- do.call(rbind, theoretical_quantiles)
      theoretical_quantiles_CI <- apply(X = theoretical_quantiles_CI,
                                        MARGIN = 2,
                                        FUN = quantile,
                                        probs = CI_probs)

      above_sim_int <- observed_quantiles_CI[1,] > theoretical_quantiles_CI[2,]
      below_sim_int <- observed_quantiles_CI[2,] < theoretical_quantiles_CI[1,]
      point_colours <- rep('black', length(p))
      point_colours[above_sim_int] <- 'red'
      point_colours[below_sim_int] <- 'blue'

      plot(
        #x = rep(std_norm_quants,4),
        x = rep(std_exp_quants,4),
        y = c(theoretical_quantiles_CI[1,],
              theoretical_quantiles_CI[2,],
              observed_quantiles_CI[1,],
              observed_quantiles_CI[2,]),
        ylab = 'sample quantile',
        xlab = 'model quantiles',
        type = 'n',
        asp = 1)
      polygon(
        #x = c(std_norm_quants, rev(std_norm_quants)),
        x = c(std_exp_quants, rev(std_exp_quants)),
        y = c(theoretical_quantiles_CI[1,], rev(theoretical_quantiles_CI[2,])),
        col = "grey60",
        border = NA
      )
      for(j in seq_along(p)){
        lines(
          #x = rep(std_norm_quants[j],2),
          x = rep(std_exp_quants[j],2),
          y = observed_quantiles_CI[,j],
          col = point_colours[j])
      }
    }
dev.off()


