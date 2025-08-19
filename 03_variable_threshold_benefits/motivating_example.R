##
# Small simulation study to show the benefits of stepped threshold
# Forms section 3 of paper on threshold selection
#
# For conservative, stepped and extended catalogues show:
#  - parameter recovery,
#  - mse decomposition,
#  - conditional return level estimation
##

SEED = 1234
DELTA = 0.1
N_SIMS = 1000

OUTPUT_PATH <- here("00_outputs", "03_motivating_example")

# 1: Source necessary code -----------------------------------------------------
library(dplyr)
library(here)
library(threshold)
source("03_helpers.R")

# 2: Simulate example catalogue ------------------------------------------------
set.seed(SEED)
## 2.1: Set underlying parameters ----------------------------------------------
cat_size <- 1000
u <- 1.05
sig_u <- 0.3
xi <- 0.1
to_nearest <- DELTA
v_cons <- rep(1.65, cat_size)
v_step <- rep(c(1.65, 1.05), each = cat_size / 2)


## 2.2: Create full catalogue --------------------------------------------------
mags_full <- rgpd_rd(
  n = cat_size,
  shift_latent = u,
  scale = sig_u,
  shape =  xi,
  to_nearest = to_nearest
)

cat_full <- data.frame(
  index = seq_along(mags_full),
  mag = mags_full,
  v_cons = v_cons,
  v_step = v_step
)

cat_step <- filter(cat_full, mag >= v_step)
cat_cons <- filter(cat_full, mag >= v_cons)

## 2.4: Plot catalogues --------------------------------------------------------
plot_catalogue <- function(cat, main =""){
  plot(x = cat$index,
       y = cat$mag,
       ylim = c(1,4),
       main = main,
       ylab = "magnitude",
       xlab = "event index",
       bty = 'n',
       pch = 16,
       cex = 0.8)
}

## Plot distinctly the catalogues that:

plot_catalogue(cat_full, main = paste0(nrow(cat_full),' happened'))

plot_catalogue(cat_cons, main = paste0(nrow(cat_cons),' in conservative model'))
lines(x = cat_cons$index, y = cat_cons$v_cons, col = "purple", lwd = 2)

plot_catalogue(cat_step, main = paste0(nrow(cat_step),' in stepped model'))
lines(x = cat_step$index, y = cat_step$v_step, col = "darkorange", lwd = 2)

## Plot these but overlaid (USED IN PAPER)
file_name <- "catalogue_plot_motivating_example.pdf"
path <- here(OUTPUT_PATH, file_name)
pdf(file = path, width = 8, height = 4)
par(mar = c(4.1,4.1,2.1,2.1))
color <- if_else(cat_full$mag >= v_step, true = 'grey60', false = "grey80")
plot(x = cat_full$index,
     y = cat_full$mag,
     pch = 16,
     cex = 0.8,
     col = color,
     bty = 'n',
     ylim = c(1,4),
     xlab = 'event index',
     ylab = 'magnitude')
lines(x = c(0,1000), y = c(1.65,1.65), lwd = 2, col = "darkorange")
lines(x = c(0,500.5,500.5,1000), y = c(1.65, 1.65,1.05,1.05), lty = 2, lwd = 2, col = "black")
dev.off()


# 3: Obtain mles for many conservative, stepped and extended catalogues  -------

#NB: extended catalogues have n_extended = n_stepped but all above v_conserative

## 3.1: Set up true values and storage -----------------------------------------

sig_v_cons <- sig_u + (v_cons - u) * xi
sig_v_step <- sig_u + (v_step - u) * xi

mles_cons <- data.frame(sig_u = rep(NA, N_SIMS), xi = rep(NA, N_SIMS))
mles_step <- data.frame(sig_u = rep(NA, N_SIMS), xi = rep(NA, N_SIMS))
mles_xtra <- data.frame(sig_u = rep(NA, N_SIMS), xi = rep(NA, N_SIMS))

# 3.2: Find conservative, stepped and extended data set mles -------------------

for (i in 1:N_SIMS) {

  # Generate magnitudes for all earthquakes that happened
  mags_full <- rgpd_rd(
    n = cat_size,
    shift_latent = u,
    scale = sig_u,
    shape =  xi,
    to_nearest = to_nearest)

  # Make full, conservative and stepped catalogues
  cat_full <- data.frame(index = seq_along(mags_full),
                         mag = mags_full,
                         v_cons = v_cons,
                         v_step = v_step)
  cat_cons <- filter(cat_full, mag > v_cons)
  cat_step <- filter(cat_full, mag > v_step)

  # Simulate some extra earthquakes for the extended catalogue
  n_extra <- nrow(cat_step) - nrow(cat_cons)
  mags_extra <- rgpd_rd(
    n = n_extra,
    scale = sig_v_cons[1],
    shape = xi,
    shift_latent = v_cons[1],
    to_nearest = to_nearest)
  v_extra <- rep(v_cons[1], n_extra)

  # Record the estimated parameters based on each threshold
  mles_cons[i,] <- mle_gpd_rd(
    x = cat_cons$mag,
    v = cat_cons$v_cons,
    sigxi = c(1,0),
    u = u,
    to_nearest = to_nearest,
    llh_val = FALSE,
    hessian = FALSE
  )

  mles_step[i,] <- mle_gpd_rd(
    x = cat_step$mag,
    v = cat_step$v_step,
    sigxi = c(1,0),
    u = u,
    to_nearest = to_nearest,
    llh_val = FALSE,
    hessian = FALSE
  )

  mles_xtra[i,] <- mle_gpd_rd(
    x = c(cat_cons$mag, mags_extra),
    v = c(cat_cons$v_cons, v_extra),
    sigxi = c(1,0),
    u = u,
    to_nearest = to_nearest,
    llh_val = FALSE,
    hessian = FALSE
  )
  if (i %% 10 == 0) print(i)
}

# 4: Plots relating to estimated parameters   ----------------------------------

## 4.1: parameter recovery plot ------------------------------------------------

file_name <- "bootstrap_mles_motivating_example.pdf"
path <- here(OUTPUT_PATH, file_name)

pdf(path, width = 5, height = 5)
par(mar = c(4.5,4.5,1.0,1.0))
plot(
  rbind(mles_cons, mles_step, mles_xtra),
  type = 'n',
  xlab = expression(hat(sigma)[u]),
  ylab = expression(hat(xi)), bty = 'n')
  points(mles_cons, pch = 16, cex = 0.8 , col = "gray30") #darkest
  points(mles_xtra, pch = 16, cex = 0.8 , col = "gray60") #lightest
  points(mles_step, pch = 16, cex = 0.8 , col = "gray70")
  points(x = sig_u, y = xi, pch = "+")
  legend('topright',
         c("Conservative", "Extended", "Stepped", "True"),
         col = c("grey30", "grey60", "grey70","black"),
         pch = list(16,16,16,43),
         cex = rep(1.2,4),
         bty = 'n')
dev.off()

## 4.2: MSE decomposition plot -------------------------------------------------

# Calculate the MSE decomposition under each threshold
MSE_cons <- split_MSE(
  mles = mles_cons,
  pars_true = c(sig_u, xi),
  par_names = c('sigma','xi'),
  main = 'MSE decomposition (cons)',
  plot_it = FALSE)

MSE_step <- split_MSE(
  mles = mles_step,
  pars_true = c(sig_u, xi),
  par_names = c('sigma','xi'),
  main = 'MSE decomposition (step)',
  plot_it = FALSE)

MSE_xtra <- split_MSE(
  mles_xtra,
  pars_true = c(sig_u, xi),
  par_names = c('sigma','xi'),
  main = 'MSE decomposition (xtra)',
  plot_it = FALSE)

# Create plot to compare the MSE decompositions under each threshold

file_name <- "mse_decomposition_motivating_example.pdf"
path <- here(OUTPUT_PATH, file_name)

pdf(file = path, width = 5, height = 5)
par(mar = c(4.5,4.5,1.0,1.0))
MSE_df <- compare_MSEs(
  MSEs = list(MSE_cons, MSE_xtra, MSE_step),
  legend_text = c("Conservative", "Extended", "Stepped"),
  legend_cex =  1.2,
  bar_names = c( expression(bias^2*(hat(sigma)[u])),
                 expression(bias^2*(hat(xi))),
                 expression(Var(hat(sigma)[u])),
                 expression(Var(hat(xi)))),
  bar_colors = c("grey30","grey80", "grey60"),
  las = 1
)
dev.off()


## 4.3: Calculate MSE with each method  and save as image ----------------------

# Get reduction factors as compared to conservative

names(MSE_df) <- stringr::str_remove(names(MSE_df),"\ \\*")

MSE_df$MSE <- apply(MSE_df, MARGIN = 1, FUN = sum)

MSE_df <- MSE_df %>%
  mutate(MSE_reduction_factor = MSE / max(MSE)) %>%
  add_rownames(var = "threshold_type")

file_name <- "mse_decomposition_motivating_example.csv"
path <- here(OUTPUT_PATH, file_name)
readr::write_csv(MSE_df, file = path)


# 5: Return level plots --------------------------------------------------------

## 5.1: Calculate point estimates and confidence intervals ----------------------
return_periods <- 1/(0.1^(seq(0,3,length.out = 101)))

return_level(scale = mles_cons$sig_u[1], shape = mles_cons$xi[1], shift = 1.05, v = 1.45, period = return_periods)

rl_cons <- mle_return_levels(mles_cons)
rl_step <- mle_return_levels(mles_step)
rl_xtra <- mle_return_levels(mles_xtra)

rl_CI_cons <- return_level_confidence_interval(rl_cons)
rl_CI_step <- return_level_confidence_interval(rl_step)
rl_CI_xtra <- return_level_confidence_interval(rl_xtra)

rl_true <- qgpd(
  p = 1 - 1/return_periods,
  scale = sig_u + xi * (v_cons[1] - u),
  shape = xi,
  shift = v_cons[1])

## 5.2: Create return level plot ----------------------------------------------

file_name <- "return_levels_motivating_example.pdf"
path <- here(OUTPUT_PATH, file_name)

pdf(path, width = 7, height = 5)
par(mar = c(4.1,4.1,1.1,1.1))
plot(
  y = c(rl_CI_cons$upper_limit, rl_CI_cons$lower_limit),
  x = rep(return_periods, 2),
  log = 'x',
  type = 'n',
  ylab = "conditional return level",
  xlab = "return period",
  bty = "n")
CI_lines(rl_CI_cons, col = "grey30", lwd = 2, lty = 2)
CI_lines(rl_CI_xtra, col = "grey70", lwd = 2, lty = 3)
CI_lines(rl_CI_step, col = "grey50", lwd = 2, lty = 4)
lines(x = return_periods, y = rl_true, col = "black", lwd = 2)
legend('topleft',
       c("True","Conservative", "Extended", "Stepped"),
       col = c("black","grey40", "grey60", "grey80"),
       lty = c(1,2,3,4),
       lwd = 2,
       cex = 1.2,
       bty = 'n')
dev.off()









