##  types_of_censoring_plots_LOCAL/main.r
###
## Creates plots of example catalogues with hard and phased censoring
##
## TO BE RUN LOCALLY
##
###
## Author: Zak Varty
##

jobid  <- 1


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

source("./src/get_mle_and_bootstrap.R")
source("./src/get_dq1.R")

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
## Set seeds for each catalogue (will use 1-500) ----
###
set.seed(1234)
seeds <- sample(1:1e7, size = 1e3, replace = FALSE)
print("seeds 1:10 are:")
seeds[1:10]

#_______________________________________________________________________________
###
## Set true parameters and simulate rounded GPD data ----
###

# True parameters
sig_u <- 0.55
xi  <- - 0.1
u   <- 0
v <- c(0.83,0.42)
sig_0 <- sig_u + (0 - u) * xi
n   <- 4000
to_nearest <- 0.1

###
## Simulate rounded GPD data ----
###
cat_id <- 1
# set the seed for this catalogue
  set.seed(seeds[cat_id])
    # Simulate catalogue (reproducibly)
    m_raw <- rgpd(n = n, scale = sig_u, shape = xi, mu = u)
    cat_full <- data.frame(
      m_raw = m_raw,
      u = rep(u, n),
      v = rep(v, each = n/2),
      index = seq_along(m_raw))

    # reduce to rounded, censored catalogue
    cat_hard <- cat_full %>%
      filter(m_raw >= v) %>%
      mutate(m_rounded = round_to_nearest(m_raw, to_nearest),
             index = seq_along(m_rounded)) %>%
      select(m = m_rounded, u, v, index)

    cat_phased <- cat_full %>%
      mutate(distance_below = pmax(v - m_raw, 0),
             detection_prob = pexp(q = distance_below ,rate = 7, lower.tail = FALSE),
             u = runif(n),
             keep = u < detection_prob) %>%
      filter(keep) %>%
      mutate(m_rounded = round_to_nearest(m_raw, to_nearest),
             index = seq_along(m_rounded)) %>%
      select(m = m_rounded, u, v, index)

    pdf_path <- "./Output/plots/example_cat_stepped_hard.pdf"
    pdf(pdf_path, width = 7, height = 5)
    par(mar = c(5.1,5.1,2.1,1.1))
    plot(cat_hard$index,
         cat_hard$m,
         pch = 16,
         cex.axis = 1.5,
         cex.lab = 1.5,
         cex = 0.8,
         bty = 'n',
         xlab = 'event index',
         ylab = 'magnitude', ylim = c(0,3.4))
    lines(x = cat_hard$index, y = cat_hard$v, col = 2, lwd = 3)
    dev.off()

    pdf_path <- "./Output/plots/example_cat_stepped_phased.pdf"
    pdf(pdf_path, width = 7, height = 5)
    par(mar = c(5.1,5.1,2.1,1.1))
    plot(cat_phased$index,
         cat_phased$m,
         pch = 16,
         cex.axis = 1.5,
         cex.lab = 1.5,
         cex = 0.8,
         bty = 'n',
         xlab = 'event index',
         ylab = 'magnitude', ylim = c(0,3.4))
    lines(x = cat_phased$index, y = cat_phased$v, col = 2, lwd = 3)
    dev.off()


    ###
    ## Save catalouge size and  true step lengths ----
    ###
    n_tot <- NROW(cat)
    n_1 <- sum(cat$v == v[1])
    n_2 <- sum(cat$v == v[2])


    ###
    ## Set up Bayesian optimization ----
    ###
    # Bounds on each threshold parameter
    bounds <- list(
      threshold_1 = c(0,1.5),
      threshold_2 = c(0,1.5),
      change_location = c(100,1000))

    # Initial evaluation points (chosen randomly)
    n_init <- 20
    initGrid <- data.frame(
      threshold_1 = runif(n = n_init, min = bounds[[1]][1], max = bounds[[1]][2]),
      threshold_2 = runif(n = n_init, min = bounds[[2]][1], max = bounds[[2]][2]),
      change_location = runif(n = n_init, min = bounds[[3]][1], max = bounds[[3]][2]))

    # Output container
    opt_object <- bayesOpt(
      FUN = FUN_3
      , bounds = bounds
      , initGrid = initGrid
      , iters.n = 1
      , plotProgress = FALSE
      , verbose = 2
      , acq = "ei"
    )

    ###
    ## Do Bayesian optimization ----
    ###
    iterations <- 79
    opt_object <- addIterations(opt_object, iters.n = iterations, verbose = 2, plotProgress = FALSE)

    # save output
    rds_path <- paste0("./Output/data/opt_objects/opt_object_",cat_id,".RDS")
    saveRDS(object = opt_object, file = rds_path)

    # save best parameters found
    rds_path <- paste0("./Output/data/selected/selected_threshold_parameters_",cat_id,".RDS")
    saveRDS(object = as.data.frame(getBestPars(opt_object)), file = rds_path)

    print(paste0("DONE cat_id = ", cat_id))
  })
}

print(paste0("DONE jobid = ", jobid))
