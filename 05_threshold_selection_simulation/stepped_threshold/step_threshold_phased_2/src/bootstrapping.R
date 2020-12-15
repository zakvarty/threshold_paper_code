#=====================================================================
# Functions for simulating bootstrap catalogues and mles for rounded GPD data
# Author: Zak Varty
# Date: 2020-08-11
#=====================================================================

#' Generate a list of bootstrapped catalogues
#'
#' @param gpd_mle parameter vector (sigma_u, xi) for undelying GPD exceedances of u
#' @param u  Constant threshold value for which MLE is given
#' @param n_v_mle Expected number of exceedances of v(\tau) in the original catalgoue
#' @param v_step_locations vector in index time of start, step and end times. (\tau_0, \tau_1,\tau_{max})
#' @param v_step_threshold modelling threshold values on diff(v_step_locations)
#' @param to_nearest Level of rounding applied to underlying GPD values
#' @param n_sim Number of catalogues to simulate
#' @param verbose logical, print progress to console?
#'
#' @return list of data frames giving the magnitude, index time and threshold of each bootstrapped event
#'
#' @example
#' boot_cat <- bootstrap_catalogues_step(
#'  gpd_mle = c(0.55,0.1),
#'  u = 0.95,
#'  w = sample(x = c(rep(1,63), rep(0,41), rep(0.267,14)), size = 63+41+14,replace = FALSE),
#'  v_step_points = c(0, 50.5, 100),
#'  v_step_values = c(1.44,1.07),
#'  to_nearest = 0.1,
#'  n_cats = 2,
#'  verbose = TRUE)
#'
#'
bootstrap_catalogues_step <- function(gpd_mle, u, w, v_step_points, v_step_values, to_nearest, n_cats, verbose = TRUE) {

  # Extract mles
  sig_u_mle <- gpd_mle[1]
  xi_mle  <- gpd_mle[2]

  # Get Bootstrap catalogue sizes:
  # Sample observed count
  #on_Av_samples <- rbinom(n = length(w), size = rep(1, length(w)), prob = w)
  on_Av_samples <- matrix(
    data = rbinom(n = length(w) * n_cats, size = rep(1, length(w)), prob = w),
    nrow = n_cats,
    byrow = TRUE)
  n_v_samples <- rowSums(on_Av_samples)
  # Bootstrap sample MLE of Poission rate \tilde\Lambda(A_v)
  cat_size_mles <- rpois(n = n_cats, lambda = n_v_samples)
  # Sample number of exceedances in each bootstrap catalogue
  cat_sizes <- rpois(n = n_cats, lambda = cat_size_mles)

  # Get probability of sample falling within each step
  step_lower <- v_step_points[-length(v_step_points)]
  step_upper <- v_step_points[-1]
  step_width <- step_upper - step_lower
  step_count <- length(step_width)
  step_prob <- step_width * (1 - pgpd(q = v_step_values, scale = sig_u_mle, shape = xi_mle, mu =  u))

  # Set up storage for bootstrap catalogues
  boot_cats <- vector(mode = 'list', length = n_cats)

  for(i in seq_along(boot_cats)){
    cat <- data.frame(m = rep(NA, cat_sizes[i]),
                      tau = rep(NA, cat_sizes[i]),
                      v = rep(NA, cat_sizes[i]))
    cat_size <- cat_sizes[i]
    # sample the step to which each event belongs
    sampled_steps <- sample(
      x = 1:step_count,
      size =  cat_size,
      replace = TRUE,
      prob = step_prob)

    # sample times uniformly within those steps
    U <- runif(n =  cat_size)
    sampled_times <- step_lower[sampled_steps] + U * step_width[sampled_steps]

    # look up threshold value corresponding to sampled step
    sampled_thresholds <- v_step_values[sampled_steps]

    # sample magnitude conditional on time and threshold
    sampled_magnitudes <- rgpd(
      n = cat_size,
      shape = xi_mle,
      scale = sig_u_mle + xi_mle * (sampled_thresholds - u),
      mu = sampled_thresholds)

    sampled_magnitudes_rounded <- round_to_nearest(
      x = sampled_magnitudes,
      to_nearest = to_nearest)

    # collate in data.frame
    cat <- data.frame(
      m = sampled_magnitudes_rounded,
      v = sampled_thresholds,
      tau = sampled_times,
      m_raw = sampled_magnitudes)

    cat <- dplyr::arrange(cat, sampled_times)

    boot_cats[[i]] <- cat

    if(verbose) print(paste("simulated", i, "of", n_cats))
  }

  return(boot_cats)
}

#' Get expected number of exceedances of modelling threshold from rounded catalgoue
#'
#' @param x Vector of rounded magnitude values
#' @param v Vector of modelling threshold at each rounded magnitude
#' @param to_nearest Level at which rounding is applied to x
#' @param gpd_mle Maximum likelihood estimates of undelying GPD parameters (sigma_u, xi)
#' @param u Constant threshold for which MLE values are given
#'
#' @return scalar numeric, expected number of exceedances of v in catalogue.
#'
#' @examples
#' # simulate catalogue and set example modelling thresholds
#'   y <- rgpd(n = 1000,shape = 0.1, scale = 0.55, mu = 0.95)
#'   x <- round_to_nearest(y, to_nearest = 0.1)
#'   v_1 <- rep(1.05, 1000)
#'   v_2 <- rep(1.02, 1000)
#' # returns sample count for bin edges
#'   sum(y > 1.05)
#'   sum(x > 1.05)
#'   get_n_v_mle(x = x, v = v_1, to_nearest = 0.1, gpd_mle = c(0.55,0.95), u = 0.95)
#' # Interpolates according to fitted gpd when threshold is between bin edges
#'   sum(y > 1.02)
#'   sum(x > 1.02)
#'   get_n_v_mle(x = x, v = v_2, to_nearest = 0.1, gpd_mle = c(0.55,0.95), u = 0.95)
#'
get_n_v_mle <- function(x, v, to_nearest, gpd_mle, u){
  # extract and format parameters
  sig_u_mle <- rep(gpd_mle[1], length(x))
  xi_mle <- rep(gpd_mle[2], length(x))
  u <- rep(u, length(x))

  # bin edges
  delta <- to_nearest / 2
  x_high <- x + delta
  x_low <-  pmax(v, x - delta)
  x_bottom <- x - delta

  # CDF values
  p_high <- pgpd(q = x_high, shape = xi_mle, scale = sig_u_mle, mu = u)
  p_low <-  pgpd(q = x_low, shape = xi_mle, scale = sig_u_mle, mu = u)
  p_bottom <- pgpd(q = x_bottom, shape = xi_mle, scale = sig_u_mle, mu = u)

  # Probability each observation is above threshold
  p_above_v <- (p_high - p_low) / (p_high - p_bottom)

  return(sum(p_above_v))
}

#' Get expected number of exceedances of modelling threshold from rounded catalgoue
#'
#' @param x Vector of rounded magnitude values
#' @param v Vector of modelling threshold at each rounded magnitude
#' @param to_nearest Level at which rounding is applied to x
#' @param gpd_mle Maximum likelihood estimates of undelying GPD parameters (sigma_u, xi)
#' @param u Constant threshold for which MLE values are given
#'
#' @return scalar numeric, expected number of exceedances of v in catalogue.
#'
#' @examples
#' # simulate catalogue and set example modelling thresholds
#'   y <- rgpd(n = 1000,shape = 0.1, scale = 0.55, mu = 0.95)
#'   x <- round_to_nearest(y, to_nearest = 0.1)
#'   v_1 <- rep(1.05, 1000)
#'   v_2 <- rep(1.02, 1000)
#' # returns sample count for bin edges
#'   sum(y > 1.05)
#'   sum(x > 1.05)
#'   w <- get_w_vector(x = x, v = v_1, to_nearest = 0.1, gpd_mle = c(0.55,0.95), u = 0.95)
#'   w
#'   sum(w)
#' # Interpolates according to fitted gpd when threshold is between bin edges
#'   sum(y > 1.02)
#'   sum(x > 1.02)
#'   w <- get_w_vector(x = x, v = v_2, to_nearest = 0.1, gpd_mle = c(0.55,0.95), u = 0.95)
#'   w
#'   sum(w)
get_w_vector <- function(x, v, to_nearest, gpd_mle, u){
  # extract and format parameters
  sig_u_mle <- rep(gpd_mle[1], length(x))
  xi_mle <- rep(gpd_mle[2], length(x))
  u <- rep(u, length(x))

  # bin edges
  delta <- to_nearest / 2
  x_high <- x + delta
  x_low <-  pmax(v, x - delta)
  x_bottom <- x - delta

  # CDF values
  p_high <- pgpd(q = x_high, shape = xi_mle, scale = sig_u_mle, mu = u)
  p_low <-  pgpd(q = x_low, shape = xi_mle, scale = sig_u_mle, mu = u)
  p_bottom <- pgpd(q = x_bottom, shape = xi_mle, scale = sig_u_mle, mu = u)

  # Probability each observation is above threshold
  w <- (p_high - p_low) / (p_high - p_bottom)

  return(w)
}


#' loglikelihood for rounded GPD data above v, where all border points are above v
#'
#' @param sigxi GPD parameter values at which to evaluate (negative) loglikelihood
#' @param u threshold for given GPD parameters
#' @param v modelling threshold at each observation
#' @param x  rounded GPD observations
#' @param to_nearest level of rounding applied to GPD observations
#' @param negative logical. Return negative loglikelihood?
#'
#' @return scalar loglikelihood value
#'
bootstrap_llh_gpd_rd_varu <- function(sigxi, u, v, x, to_nearest, negative){

  # Check all scale params positive
  if(sigxi[1] <= 0){ llh <- -10e6 ; return((-1)^negative * llh)}

  # check x all be above v (adjusted for rounding)
  lep_fail <- any(x <= v - to_nearest/2)
  if(lep_fail){stop('lower endpoint failure:  !all(x > v - to_nearest / 2).')}

  # check all x below UEP (if it exists, adjusting for rounding)
  latent_uep <- NULL
  uep_fail <- FALSE
  if(sigxi[2] < 0){
    latent_uep <- u - sigxi[1] / sigxi[2]
    uep_fail <- max(x) >= (latent_uep + to_nearest / 2)
  }
  if(!is.null(latent_uep) & uep_fail){
    llh <- -10e6
    return((-1)^negative * llh)
  }

  #format inputs
  u <- rep(u, length(x))
  xi <- rep(sigxi[2], length(x))
  sig_u <- rep(sigxi[1], length(x))
  sig_v <- sig_u + xi * (v - u)

  # calculate llh
  x_high <- x + to_nearest / 2
  x_low <- pmax(x - to_nearest / 2, v)

  p_high <- pgpd(q = x_high, shape = xi, scale = sig_v, mu = v)
  p_low  <- pgpd(q = x_low , shape = xi, scale = sig_v, mu = v)

  p_in <- p_high - p_low
  llh <- sum(log(p_in))

  return(((-1)^negative) * llh)
}

# bootstrap_llh_gpd_rd_varu(
#   sigxi = c(1,0),
#   u =  0.9,
#   v = boot_cat[[1]]$v,
#   x = boot_cat[[1]]$m,
#   to_nearest = 0.1,
#   negative = TRUE)
# bootstrap_mle_gpd_rd_varu(bootstrap_x = boot_cat[[1]]$m, bootstrap_v = boot_cat[[1]]$v, sigxi_init = c(1,0),u = 0.95, to_nearest = 0.1, llh_val = TRUE, hessian = TRUE)

bootstrap_mle_gpd_rd_varu <- function(bootstrap_x,bootstrap_v, sigxi_init, u,  to_nearest, llh_val = FALSE, hessian = FALSE) {

  # Numerically minimise negative log-likelihood
  temp <- optim(par = sigxi_init,
                fn = bootstrap_llh_gpd_rd_varu,
                u = u,
                v = bootstrap_v,
                x = bootstrap_x,
                to_nearest = to_nearest,
                negative = TRUE,
                hessian = hessian)

  #format output
  if(!llh_val & !hessian){out <-  temp$par} else {out <- list(params = temp$par)}
  if(llh_val) out$loglik <- -temp$value
  if(hessian) out$hessian <- -temp$hessian

  return(out)
}


