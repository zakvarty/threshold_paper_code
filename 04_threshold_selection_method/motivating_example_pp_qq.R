##
# Making transformed QQ and PP plots for motivating example
##

## 1: Source necessary code ----
source('src/_gpd.R')
source('src/round_to_nearest.R')

source('src/llh_gpd_rd_varu.R')
source('src/mle_gpd_rd_varu.R')

source('src/sample_mles_gpd_rd_stepped.R')
source("src/thresholding/sample_latent_magnitudes.R")
source("src/thresholding/standardise_gpd_sample.R")

library(dplyr)
library('purrr')
#source('src/split_MSE.R')
# source('src/gpd_rd_to_unif.R')
# source('src/gpd_rd_to_unif_multiple.R')
# source('src/qq_norm.R')
# source('src/QQ_interval_gpd_rd.R')
# source('src/PP_interval_gpd_rd.R')
# source('src/split_MSE.R')

# 2: Simulate example catalogue ----
## 2.1: Set underlying parameters ----
cat_size <- 1000
u <- 1.05
sig_u <- 0.3
xi <- 0.1
to_nearest <- 0.1
v_cons <- rep(1.65, cat_size)
v_step <- rep(c(1.65, 1.05), each = cat_size / 2)
## 2.2: Create full catagloue ----
set.seed(4321)
mags_full <- rgpd_rd(
  n = cat_size,
  mu_latent = u,
  sig = sig_u,
  xi =  xi,
  to_nearest = to_nearest
)
cat_full <- data.frame(
  index = seq_along(mags_full),
  mag = mags_full,
  v_cons = v_cons,
  v_step = v_step
)

## 2.3: Thin by censoring and model inclusion ----
cat_censored <- filter(cat_full, mag >= v_step)
cat_step <- cat_censored
cat_cons <- filter(cat_censored, mag >= v_cons)

## 2.4: Plot catalogues ----
pdf(file = "./output/plots/motivating_example/cat_plot_motivating_example.pdf",width = 8,height = 4)
par(mar = c(4.1,4.1,2.1,2.1))
color <- if_else(cat_full$mag >= v_step, true = 'black', false = "grey50")
plot(x = cat_full$index,
     y = cat_full$mag, pch = 16,
     cex = 0.8,
     col = color,
     bty = 'n',
     ylim = c(1,4),
     xlab = 'event index',
     ylab = 'magnitude')
lines(x = c(0,500.5,500.5,1000), y = c(1.65, 1.65,1.05,1.05), lwd = 2)
lines(x = c(0,1000), y = c(1.65,1.65), lty = 2, col = "grey50", lwd = 2)
dev.off()

## Get the number of events in conservative or stepped catalogue
n_1 <- sum(cat_full$mag >= cat_full$v_step)
n_2 <- sum(cat_full$mag >= cat_full$v_cons)


#_______________________________________________________________________________
###
## Get simulated catalogue into correct format ----
###
cat <- cat_censored %>%
  transmute(m = mag, u = u, v = v_step, index)

n_1 <- NROW(cat)
#_______________________________________________________________________________
###
## Set modelling thresholds to try ----
###

# Set thresholds
thresholds_vec <- c(0.5, 1.15, 1.85)
thresholds <- map(.x = thresholds_vec, .f = function(x){rep(x,n_1)})

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

mles <- pmap(
  .l = list(x = mags_filtered, v = thresholds_filtered),
  .f = mle_gpd_rd_varu,
  u = 0,
  sigxi = c(1,0),
  to_nearest = 0.1,
  llh_val = FALSE
)
sigma_0_mles <- map_dbl(mles, function(vec){vec[1]})
xi_mles <- map_dbl(mles, function(vec){vec[2]})

#_______________________________________________________________________________
###
## Calculate expected number of exceedances at each threshold ----
###
source('src/bootstrapping.R')
expected_exceedances <- pmap(
  .f = get_n_v_mle,
  .l = list(
    x = mags_filtered,
    v = thresholds_filtered,
    gpd_mle = mles
  ),
  to_nearest = 0.1,
  u = 0)

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

mle_samples <- vector(mode = "list", length = length(mles))
step_lengths <- map(.x = cats_filtered, NROW)
#v_values <- c(v)
for(i in seq_along(thresholds)){
  message("iteration: ", i, "/", length(thresholds))

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

#_______________________________________________________________________________
###
## Construct PP plot for each threshold mle ----
###

n_eval_pts <- 101
n_rep <- 10 ## latent samples per mle sample
CI_probs <- c(0.025,0.975)


pdf(file = "./output/plots/motivating_example/motivating_example_qq_pp_gauss_margins.pdf", width = 5, height = 5)
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

  latent_samples_gauss <- map(
    .x = latent_samples_exp,
    .f = function(x){qnorm(pexp(x))}
  )

  ## Simulate actual N(0,1) vectors of same lengths
  theo_samples_gauss <- map(
    .x = latent_samples_exp,
    .f = function(x){rnorm(n = length(x))})

  ###
  ## Construct PP plot on Gaussian margins
  ###

  # vector of x-coords for pp plot
  p <- ppoints(n_eval_pts)
  std_norm_quants <- qnorm(p = p)

  get_sample_ecdf <- function(x, q){ecdf(x)(q)}

  observed_ecdf_values <- map(
    .x = latent_samples_gauss,
    .f = get_sample_ecdf,
    q = std_norm_quants)
  observed_ecdf_values <- do.call(rbind, observed_ecdf_values)
  observed_ecdf_CI <- apply(observed_ecdf_values, 2, quantile, probs = CI_probs)

  theoretical_ecdf_values <- map(
    .x = theo_samples_gauss,
    .f = get_sample_ecdf,
    q = std_norm_quants
  )
  theoretical_ecdf_values <- do.call(rbind, theoretical_ecdf_values)
  theoretical_ecdf_CI <- apply(theoretical_ecdf_values, 2, quantile, probs = CI_probs)

  above_sim_int <- observed_ecdf_CI[1,] > theoretical_ecdf_CI[2,]
  below_sim_int <- observed_ecdf_CI[2,] < theoretical_ecdf_CI[1,]
  point_colours <- rep('black', length(p))
  point_colours[above_sim_int] <- 'red'
  point_colours[below_sim_int] <- 'blue'

  plot(
    x = rep(p,4),
    y = c(theoretical_ecdf_CI[1,], theoretical_ecdf_CI[2,], observed_ecdf_CI[1,], observed_ecdf_CI[2,]),
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
  std_norm_quants <- qnorm(p)

  observed_quantiles <- map(
    .x = latent_samples_gauss,
    .f = quantile,
    probs = p,
    na.rm = TRUE
  )
  observed_quantiles_CI <- do.call(rbind, observed_quantiles)
  observed_quantiles_CI <- apply(observed_quantiles_CI, 2, quantile, probs = CI_probs)

  theoretical_quantiles <- map(
    .x = theo_samples_gauss,
    .f = quantile,
    probs = p,
    na.rm = TRUE
  )
  theoretical_quantiles_CI <- do.call(rbind, theoretical_quantiles)
  theoretical_quantiles_CI <- apply(theoretical_quantiles_CI, 2, quantile, probs = CI_probs)

  above_sim_int <- observed_quantiles_CI[1,] > theoretical_quantiles_CI[2,]
  below_sim_int <- observed_quantiles_CI[2,] < theoretical_quantiles_CI[1,]
  point_colours <- rep('black', length(p))
  point_colours[above_sim_int] <- 'red'
  point_colours[below_sim_int] <- 'blue'

  plot(
    x = rep(std_norm_quants,4),
    y = c(theoretical_quantiles_CI[1,], theoretical_quantiles_CI[2,], observed_quantiles_CI[1,], observed_quantiles_CI[2,]),
    ylab = 'sample quantile',
    xlab = 'model quantiles',
    type = 'n',
    asp = 1)
  polygon(
    x = c(std_norm_quants, rev(std_norm_quants)),
    y = c(theoretical_quantiles_CI[1,], rev(theoretical_quantiles_CI[2,])),
    col = "grey60",
    border = NA
  )
  for(j in seq_along(p)){
    lines(
      x = rep(std_norm_quants[j],2),
      y = observed_quantiles_CI[,j],
      col = point_colours[j])
  }
}
dev.off()

pdf(file = "./output/plots/motivating_example/motivating_example_qq_pp_exp_margins.pdf", width = 5, height = 5)
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
  theo_samples_gauss <- map(
    .x = latent_samples_exp,
    .f = function(x){rnorm(n = length(x))})

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
  theoretical_ecdf_CI <- apply(theoretical_ecdf_values, 2, quantile, probs = CI_probs)

  above_sim_int <- observed_ecdf_CI[1,] > theoretical_ecdf_CI[2,]
  below_sim_int <- observed_ecdf_CI[2,] < theoretical_ecdf_CI[1,]
  point_colours <- rep('black', length(p))
  point_colours[above_sim_int] <- 'red'
  point_colours[below_sim_int] <- 'blue'

  plot(
    x = rep(p,4),
    y = c(theoretical_ecdf_CI[1,], theoretical_ecdf_CI[2,], observed_ecdf_CI[1,], observed_ecdf_CI[2,]),
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
  observed_quantiles_CI <- apply(observed_quantiles_CI, 2, quantile, probs = CI_probs)

  theoretical_quantiles <- map(
    #.x = theo_samples_gauss,
    .x = theo_samples_exp,
    .f = quantile,
    probs = p,
    na.rm = TRUE
  )
  theoretical_quantiles_CI <- do.call(rbind, theoretical_quantiles)
  theoretical_quantiles_CI <- apply(theoretical_quantiles_CI, 2, quantile, probs = CI_probs)

  above_sim_int <- observed_quantiles_CI[1,] > theoretical_quantiles_CI[2,]
  below_sim_int <- observed_quantiles_CI[2,] < theoretical_quantiles_CI[1,]
  point_colours <- rep('black', length(p))
  point_colours[above_sim_int] <- 'red'
  point_colours[below_sim_int] <- 'blue'

  plot(
    #x = rep(std_norm_quants,4),
    x = rep(std_exp_quants,4),
    y = c(theoretical_quantiles_CI[1,], theoretical_quantiles_CI[2,], observed_quantiles_CI[1,], observed_quantiles_CI[2,]),
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


