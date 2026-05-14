pp_gpd <- function(x, u, seed = NULL, n_boot = 500, verbose = FALSE, alpha = 0.05, jitter = 0, q = NULL, plotting_colour = rgb(1, 0.55, 0)){

  if ( is.null(q) ) { q <- seq(from = u, to = u + 1.2 * (max(x) - u), length.out = 201) }

  # allow local seed
  old <- .Random.seed
  on.exit( { .Random.seed <<- old } )
  if (!is.null(seed)) { set.seed(seed) }

  y <- x[x > u]
  n_y <- length(y)

  # prepare plotting data ------------------------------------------------------
  mle <- optim(par = c(1,0), fn = llh_gpd, method = "L-BFGS-B", u = u, x = y, negative = TRUE)

  bootstrap <- list(
    id =  1:n_boot,
    n_excesses = rpois(n = n_boot, lambda = n_y),
    data = vector("list", length = n_boot),
    mle = vector("list", length = n_boot))

  for (i in 1:n_boot) {
    bootstrap$data[[i]] <- sample(y, size = bootstrap$n_excesses[i], replace = TRUE)
    bootstrap$mle[[i]] <- optim(
      par = c(1,0),
      fn = llh_gpd,
      method = "L-BFGS-B",
      u = u,
      x = bootstrap$data[[i]],
      negative = TRUE)
    if (verbose) { cat(i, "\n") }
  }

  boot_sigmas <- sapply(X = bootstrap$mle, FUN = function(a){a$par[1]})
  boot_xis <- sapply(X = bootstrap$mle, FUN = function(a){a$par[2]})

  mcdf <- ecdf(y + runif(n = n_y, min = -jitter/2, max = jitter/2))

  sample_prob <- mcdf(q)
  model_prob <- pgpd(q = q, scale = mle$par[1], shape = mle$par[2], mu = u)

  boot_probs <- purrr::map2(boot_xis, boot_sigmas, .f = pgpd, q = q, mu = u)
  boot_probs_mat <- do.call(rbind, boot_probs)
  boot_probs_CI <- apply(boot_probs_mat, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2))

  # construct plot -------------------------------------------------------------

  interval_colour <- colorspace::adjust_transparency(plotting_colour, alpha = 0.3)

  plot(
    x = model_prob,
    y = sample_prob,
    xlab = "model probability",
    ylab = "sample probability",
    pch = 20,
    col = rgb(0, 0, 0, 0.3),
    asp = 1,
    xlim = c(0, 1),
    ylim = c(0, 1),
    bty = "n")
  polygon(
    x = c(model_prob, rev(model_prob)),
    y = c(boot_probs_CI[1, ], rev(boot_probs_CI[2, ])),
    col = interval_colour,
    border = NA,
    lty = NULL)
  abline(a = 0, b = 1, col = plotting_colour)

  # quietly return data from plot ----------------------------------------------
  plotting_df <- data.frame(
    q = q,
    model_probability = model_prob,
    sample_probability = sample_prob,
    lower = boot_probs_CI[1, ],
    upper = rev(boot_probs_CI[2, ]))

  invisible(
    list(
      plotting_df = plotting_df,
      mle_sigma = mle$par[1],
      mle_xi = mle$par[2],
      boot_sigmas = boot_sigmas,
      boot_xis = boot_xis))
}
