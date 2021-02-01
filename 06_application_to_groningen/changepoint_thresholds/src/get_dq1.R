
#' Sample unrounded threshold exceedances and calculate d(q,1)
#'
#' @param cat catalogue with all observed magnitudes in column "m"
#' @param threshold vector of threshold values at each event in cat
#' @param to_nearest level of rounding applied to cat$m
#' @param u value at which to calculate scale mle (\hat \sigma_u, \hat\xi)
#' @param mle_obj mle-object obtained using get_mle_and_bootstrap
#' @param n_eval_pts number of equally spaced probabilities at which to evaluate quantiles
#' @param y_samples_per_mle number of times to sample unrounded values per bootstrap mle
#'
#' @return numeric value for d(q,1) for the given catalogue, threshold and set of bootstrap mles.
#'
get_dq1 <- function(cat, threshold, to_nearest, u, mle_obj, n_eval_pts = 501, y_samples_per_mle = 10){

  # 1) filter catalogue to those with positive prob of being above threshold
  cat_filtered <- cat %>%
    mutate(threshold = threshold) %>%
    filter(m  > (threshold - to_nearest/2 + 1e-8))

  # 2) Calculate probability each filtered event is above threshold
  w_vector <- get_w_vector(
    x = cat_filtered$m,
    v = cat_filtered$threshold,
    to_nearest = to_nearest,
    gpd_mle = mle_obj$mle,
    u = u)

  # 3)
  eval_probs <- ppoints(n_eval_pts) # evaluation points equally spaced on probability scale
  model_quants <- qexp(p = eval_probs, rate = 1)
  ## Calculates discrepancy from Exp(1) qq-plot
  qq_dist <- function(sample_q,  model_q){ifelse(test = is.finite(sample_q), yes = model_q - sample_q, no =  NA)}

  # 4) sample unrounded magnitudes
  y_samples <- pmap(
    .f = sample_latent_magnitudes,
    .l = list(sig_u = rep(mle_obj$boot_mles[,1], y_samples_per_mle),
              xi = rep(mle_obj$boot_mles[,2], y_samples_per_mle)),
    x = cat_filtered$m,
    u = u,
    to_nearest = to_nearest,
    v = cat_filtered$threshold)

  # 5) Transform unrounded magnitudes to Exp(1) margins
  z_samples <- pmap(
    .f = standardise_gpd_sample,
    .l = list(
      y = y_samples,
      sig_u = rep(mle_obj$boot_mles[,1], y_samples_per_mle),
      xi = rep(mle_obj$boot_mles[,2], y_samples_per_mle)),
    u = u,
    v = cat_filtered$threshold)

  # 6) For each set of Exp(1) obs, calculate d_i(q,1).
  sample_quants <- map(.f = quantile, .x = z_samples, probs = eval_probs, na.rm = TRUE)
  qq_deviation_vectors <- map(.x = sample_quants, .f = qq_dist, model_q = model_quants)
  qq_MAE <- map_dbl(.x = qq_deviation_vectors, .f = function(d){mean(abs(d), na.rm = TRUE)})

  # 7) Calculate the expected value over simulations of underlying dataset
  qq_EMAE <- mean(qq_MAE)
  return(qq_EMAE)
}
