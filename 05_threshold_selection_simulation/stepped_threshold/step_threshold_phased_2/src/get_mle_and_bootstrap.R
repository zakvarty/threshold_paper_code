#' Get MLE of GPD parameters for rounded observations and stepped threshold
#'
#' @param cat catalogue with all observed magnitudes in column "m"
#' @param threshold vector of threshold values at each event in cat
#' @param to_nearest level of rounding applied to cat$m
#' @param u value at which to calculate scale mle (\hat \sigma_u, \hat\xi)
#' @param n_bootstrap_mles number of bootstrap MLE values to simulate
#' @param verbose Print map timings when bootstrapping?
#'
#' @return list of mle, llh and optionally boot_mles
#'
get_mle_and_bootstrap <- function(cat, threshold, to_nearest, u, n_bootstrap_mles = 0, verbose = TRUE ){

  # 1) filter catalogue to those with positive prob of being above threshold
  cat_filtered <- cat %>%
    mutate(threshold = threshold) %>%
    filter(m  > (threshold - to_nearest/2 + 1e-8))

  # 2) calculate mle (\hat sigma_u, \hat xi) and log-likelihood at mle
  MLE <- mle_gpd_rd_varu(
    sigxi = c(1,0),
    u = u,
    x = cat_filtered$m,
    v = cat_filtered$threshold,
    to_nearest = to_nearest,
    llh_val = TRUE)

  out <- list(mle = MLE$params, llh = MLE$loglik, boot_mles = NULL)

  # 3) If needed, get bootstrap sample of mle values
  if(n_bootstrap_mles > 0){

    # 3.1) Calculate prob each filtered event above threshold
    w_vector <- get_w_vector(
      x = cat_filtered$m,
      v = cat_filtered$threshold,
      to_nearest = to_nearest,
      gpd_mle = MLE$params,
      u = u
      )

    # 3.2) Get the number of exceedances of each unique threhsold value
    threshold_table <- table(cat_filtered$threshold)
    step_lengths <- as.double(threshold_table)
    step_values <- as.double(dimnames(threshold_table)[[1]])

    # 3.3) Simulate parameterc bootstrap MLE values
    mle_samples <- sample_mles_gpd_rd_stepped(
      mle = MLE$params,
      u = u,
      step_lengths = step_lengths,
      step_values = step_values,
      w = w_vector,
      to_nearest = to_nearest,
      n_sim = n_bootstrap_mles,
      verbose = verbose)

    out$boot_mles <- t(mle_samples)
  }

return(out)
}
