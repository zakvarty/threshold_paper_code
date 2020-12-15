modified_pp_plot <- function(m, threshold, mle_obj, u, to_nearest, n_rep, n_eval_pts = 501, CI_alpha = 0.05, ...){

  ## Sample true magnitudes above that threshold
  latent_samples <- purrr::pmap(
    .l = list(sig_u = rep(mle_obj$boot_mles[,1], n_rep),
              xi = rep(mle_obj$boot_mles[,2], n_rep)),
    .f = sample_latent_magnitudes,
    x = m,
    v = threshold,
    u = u,
    to_nearest = to_nearest
  )

  ## Transform to have common Exp(1) distribution
  latent_samples_exp <- purrr::pmap(
    .f = standardise_gpd_sample,
    .l = list(
      y = latent_samples,
      sig_u = rep(mle_obj$boot_mles[,1], n_rep),
      xi = rep(mle_obj$boot_mles[,2], n_rep)),
    v = threshold,
    u = u
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
  ## Construct PP plot on Gaussian margins
  ###

  # vector of x-coords for pp plot
  p <- ppoints(n_eval_pts)
  #std_norm_quants <- qnorm(p = p)
  std_exp_quants <- qexp(p = p)

  get_sample_ecdf <- function(x, q){ecdf(x)(q)}
  CI_probs <- c(CI_alpha/2, 1 - CI_alpha / 2)

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
    asp = 1,
    ...)
  polygon(
    x = c(p, rev(p)),
    y = c(theoretical_ecdf_CI[1,], rev(theoretical_ecdf_CI[2,])),
    col = "grey60",
    border = NA
  )
  for(j in seq_along(p)){
    lines(x = rep(p[j],2), y = observed_ecdf_CI[,j], col = point_colours[j])
  }
}

