##
# Small simulation study to show the benefits of stepped threshold
# Forms section 3 of paper on threshold selection
#
# For conservative, stepped and extended catalouges show:
# parameter recovery, mse decomposition and conditional return level estimation
##

# 1: Source necessary code -----------------------------------------------------
source('../00_src/_gpd.R')
source('../00_src/round_to_nearest.R')
source('../00_src/llh_gpd_rd_varu.R')
source('../00_src/mle_gpd_rd_varu.R')
source('../00_src/variable_threshold_benefits/split_MSE.R')

library(dplyr)

# 2: Simulate example catalogue ------------------------------------------------

## 2.1: Set underlying parameters ----------------------------------------------
cat_size <- 1000
u <- 1.05
sig_u <- 0.3
xi <- 0.1
to_nearest <- 0.1
v_cons <- rep(1.65, cat_size)
v_step <- rep(c(1.65, 1.05), each = cat_size / 2)

## 2.2: Create full catagloue --------------------------------------------------
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

## 2.3: Thin by censoring and model inclusion ----------------------------------
cat_censored <- filter(cat_full, mag >= v_step)
cat_step <- cat_censored
cat_cons <- filter(cat_censored, mag >= v_cons)

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
# Happened
plot_catalogue(cat_full, main = paste0(nrow(cat_full),' happened'))
# were observed
plot_catalogue(cat_censored, main = paste0(nrow(cat_censored),' observed'))
# were above the conservative threshold
plot_catalogue(cat_cons, main = paste0(nrow(cat_cons),' in conservative model'))
lines(x = cat_cons$index, y = cat_cons$v_cons, col = 2, lwd = 1.5)
# were above the stepped theshold
plot_catalogue(cat_step, main = paste0(nrow(cat_step),' in stepped model'))
lines(x = cat_step$index, y = cat_step$v_step, col = 2, lwd = 1.5)

## Plot these but overlaid (USED IN PAPER)
pdf(file = "./output/plots/cat_plot_motivating_example.pdf",width = 8,height = 4)
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
lines(x = c(0,1000), y = c(1.65,1.65), lty = 2, col = 2, lwd = 2)
dev.off()


# 3: Obtain mles for many conservative, stepped and extended catalogues  -------

#NB: extended catalgoues have n_extended = n_stepped but all above v_conserative

## 3.1: Set up true values and storage -----------------------------------------
n_sims <- 1000
sig_v_cons <- sig_u + (v_cons - u) * xi
sig_v_step <- sig_u + (v_step - u) * xi

mles_cons <- data.frame(sig_u = rep(NA, n_sims), xi = rep(NA, n_sims))
mles_step <- data.frame(sig_u = rep(NA, n_sims), xi = rep(NA, n_sims))
mles_xtra <- data.frame(sig_u = rep(NA, n_sims), xi = rep(NA, n_sims))

# 3.2: Find conservative, stepped and extended data set mles -------------------
for(i in 1:n_sims){
  # Generate magnitudes for all earthquakes that happened
  mags_full <- rgpd_rd(
    n = cat_size,
    mu_latent = u,
    sig = sig_u,
    xi =  xi,
    to_nearest = to_nearest)

  # Make full, conservative and stepped catalouges
  cat_full <- data.frame(index = seq_along(mags_full),
                         mag = mags_full,
                         v_cons = v_cons,
                         v_step = v_step)
  cat_cons <- filter(cat_full, mag >= v_cons)
  cat_step <- filter(cat_full, mag >= v_step)

  # Simulate some extra earthquakes for the extended catalogue
  n_extra <- nrow(cat_step) - nrow(cat_cons)
  mags_extra <- rgpd_rd(
    n = n_extra,
    sig = sig_v_cons[1],
    xi = xi,
    mu_latent = v_cons[1],
    to_nearest = to_nearest)
  v_extra <- rep(v_cons[1], n_extra)

  # Record the estimated parameters based on each threshold
  mles_cons[i,] <- mle_gpd_rd_varu(
    x = cat_cons$mag,
    v = cat_cons$v_cons,
    sigxi = c(1,0),
    u = u,
    to_nearest = to_nearest,
    llh_val = FALSE,
    hessian = FALSE
  )

  mles_step[i,] <- mle_gpd_rd_varu(
    x = cat_step$mag,
    v = cat_step$v_step,
    sigxi = c(1,0),
    u = u,
    to_nearest = to_nearest,
    llh_val = FALSE,
    hessian = FALSE
  )

  mles_xtra[i,] <- mle_gpd_rd_varu(
    x = c(cat_cons$mag, mags_extra),
    v = c(cat_cons$v_cons, v_extra),
    sigxi = c(1,0),
    u = u,
    to_nearest = to_nearest,
    llh_val = FALSE,
    hessian = FALSE
  )
  if(i %%10 == 0) print(i)
}

# 4: Plots relating to estimated parameters   ----------------------------------

## 4.1: parameter recovery plot ------------------------------------------------
pdf("./output/plots/mle_values_motivating_example.pdf", width = 5, height = 5)
par(mar = c(4.5,4.5,1.0,1.0))
plot(
  rbind(mles_cons, mles_step, mles_xtra),
  type = 'n',
  xlab = expression(hat(sigma)[u]),
  ylab = expression(hat(xi)), bty = 'n')
points(mles_cons, pch = 16, cex = 0.8 , col = "gray30") #darkest
points(mles_xtra, pch = 16, cex = 0.8 , col = "gray70") #lightest
points(mles_step, pch = 16, cex = 0.8 , col = "gray60")
points(x = sig_u, y = xi, pch = "+")
legend('topright',
       c("Conservative", "Extended", "Stepped", "True"),
       col = c("grey30", "grey70", "grey60","black"),
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

compare_MSEs <- function(MSEs, bar_colors = NULL, bar_names = NULL, legend_position = 'topleft', legend_text = letters[seq_along(MSEs)], legend_cex = 1,...){
  # Function to plot many MSE decompositions side-by-side
  heights <- NULL
  if(is.null(bar_names)){ bar_names <- c('bais^2(sigma)', 'bais^2(xi)', 'var(sigma)','var(xi)')}
  if(is.null(bar_colors)){ bar_colors <- grey(ppoints(length(MSEs))) }

  for(i in seq_along(MSEs)){
    heights_i <- c(MSEs[[i]]$sq_bias, MSEs[[i]]$variance)
    heights <- rbind(heights,heights_i)
  }

  barplot(height = heights, beside = TRUE, names.arg = bar_names, col = bar_colors, ...)
  legend(legend_position, legend = legend_text, col = bar_colors, pch = 15, bty = 'n', cex = legend_cex)

  heights_df <- as.data.frame(heights, row.names = legend_text)
  colnames(heights_df) <- bar_names
  invisible(heights_df)
}


pdf(file = "output/plots/mse_decomposition_motivating_example.pdf",
    width = 5,
    height = 5)
par(mar = c(4.5,4.5,1.0,1.0))
MSE_df <- compare_MSEs(
  MSEs = list(MSE_cons, MSE_xtra, MSE_step),
  #ylim = c(0,0.003),
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
MSEs <- apply(MSE_df, MARGIN = 1, FUN = sum)
1/(MSEs["Stepped"] / MSEs["Conservative"])
1/(MSEs["Extended"] / MSEs["Conservative"])

pdf("output/plots/MSE_reduction_factors_motivating_example.pdf",
    width = 5,
    height = 5)
plot.new()
text(x = c(0,0,0),
     y = c(1,0.5,0),
     cex = 1.2,
     pos = 4,
     paste0(
       "MSE reduction factor:",
       c(" conservative = ", "stepped = ", " extended = "),
       c(1,
         round(MSEs["Conservative"] / MSEs["Stepped"],3),
         round(MSEs["Conservative"] / MSEs["Extended"],3))
       )
    )
dev.off()

# 5: Return level plots --------------------------------------------------------

## 5.1: Calculate point estimaes and confidence intervals ----------------------
return_periods <- 1/(0.1^(seq(0,3,length.out = 101)))

return_levels <- function(sigxi_u, u, v, return_periods_v){
  sig_u <- sigxi_u[1]
  xi <- sigxi_u[2]
  sig_v <- sig_u + (v - u) * xi

  qgpd(p = (1-1/return_periods_v), shape = xi, scale = sig_v, mu = v)
}

return_levels(
  sigxi = as.numeric(mles_cons[1,]),
  u = 1.05,
  v = 1.45,
  return_periods_v = return_periods)

rl_cons <- apply(X = unname(as.matrix(mles_cons)), MARGIN = 1, FUN = return_levels, u = u, v = v_cons[1], return_periods_v = return_periods)
rl_step <- apply(X = unname(as.matrix(mles_step)), MARGIN = 1, FUN = return_levels, u = u, v = v_cons[1], return_periods_v = return_periods)
rl_xtra <- apply(X = unname(as.matrix(mles_xtra)), MARGIN = 1, FUN = return_levels, u = u, v = v_cons[1], return_periods_v = return_periods)

alpha <- 0.05
rl_CI_cons <- apply(X = rl_cons, MARGIN = 1, FUN = quantile, probs = c(alpha/2, 0.5, 1-alpha/2))
rl_CI_step <- apply(X = rl_step, MARGIN = 1, FUN = quantile, probs = c(alpha/2, 0.5, 1-alpha/2))
rl_CI_xtra <- apply(X = rl_xtra, MARGIN = 1, FUN = quantile, probs = c(alpha/2, 0.5, 1-alpha/2))

rl_true <- qgpd(p = 1 - 1/return_periods, scale = sig_u + xi * (v_cons[1] - u), shape = xi, mu = v_cons[1])

lines_CI <- function(ymatrix, x, col, type,...){
  for(i in 1:nrow(ymatrix)){
  lines(x = x, y = ymatrix[i,], lty = type[i], col = col[i],...)
  }
}

## 5.2: CCreate return level plot ----------------------------------------------

pdf("output/plots/return_levels_motivating_example.pdf", width = 7, height = 5)
par(mar = c(4.1,4.1,1.1,1.1))
plot(
  y = rbind(rl_CI_cons, rl_CI_step, rl_CI_xtra),
  x = rep(return_periods, each = 9),
  col = rep(1:3, each = 3),
  log = 'x',
  type = 'n',
  ylab = "conditional return level",
  xlab = "return period",
  bty = "n")
lines(x = return_periods, y = rl_true, col = 1, lwd = 2)
lines_CI(x = return_periods, ymatrix = rl_CI_cons, col = rep("grey30",3), type = c(2,2,2), lwd = 2)
lines_CI(x = return_periods, ymatrix = rl_CI_xtra, col = rep("grey70",3), type = c(3,3,3), lwd = 2)
lines_CI(x = return_periods, ymatrix = rl_CI_step, col = rep("grey50",3), type = c(4,4,4), lwd = 2)
legend('topleft',
       c("True","Conservative", "Extended", "Stepped"),
       col = c("black","grey30", "grey70", "grey50"),
       lty = c(1,2,3,4),
       lwd = 2,
       cex = 1.2,
       bty = 'n')
dev.off()









