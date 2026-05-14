qq_gpd <- function(x, u, seed = NULL, n_boot = 500, verbose = FALSE, alpha = 0.05, jitter = 0, lim = NULL, plotting_colour = rgb(1, 0.55, 0), ...){

  # Local seeding
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

  boot_sigmas <- sapply(X = bootstrap$mle, FUN = function(a){a$par[[1]]})
  boot_xis <- sapply(X = bootstrap$mle, FUN = function(a){a$par[2]})

  p_points <- ppoints(n_y)
  boot_quantiles <- purrr::pmap(
    .l = list(shape = boot_xis, scale = boot_sigmas),
    .f = qgpd,
    p = p_points,
    mu = u)
  boot_quantiles_mat <- do.call(rbind, boot_quantiles)
  boot_quantilles_CI <- apply(boot_quantiles_mat, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2))

  sample_jittered <- y + runif(n = n_y, min = -jitter / 2, max = jitter / 2)
  model_quantile <- qgpd(p = p_points, scale = mle$par[1], shape = mle$par[2], mu = u)
  sample_quantile <- sort(sample_jittered)

  # construct plot -------------------------------------------------------------

  if (is.null(lim)) lim <-  c(u, u + 1.2 * (max(x) - u))
  interval_colour <- colorspace::adjust_transparency(plotting_colour, alpha = 0.3)

  plot(
    x = model_quantile,
    y = sample_quantile,
    xlab = "model quantile",
    ylab = "sample quantile",
    pch = 20,
    col = rgb(0, 0, 0, 0.3),
    asp = 1,
    xlim = lim,
    ylim = lim,
    bty = "n",
    ...)
  polygon(
    x = c(model_quantile, rev(model_quantile)),
    y = c(boot_quantilles_CI[1, ], rev(boot_quantilles_CI[2, ])),
    col = interval_colour,
    border = NA,
    lty = NULL)
  abline(a = 0, b = 1, col = plotting_colour)

  # quietly return data from plot ----------------------------------------------
  plotting_df <- data.frame(
    p = p_points,
    model_quantile = model_quantile,
    sample_quantile = sample_quantile,
    lower = boot_quantilles_CI[1, ],
    upper = rev(boot_quantilles_CI[2, ]))

  invisible(
    list(
      plotting_df = plotting_df,
      mle_sigma = mle$par[1],
      mle_xi = mle$par[2],
      boot_sigmas = boot_sigmas,
      boot_xis = boot_xis))
}
