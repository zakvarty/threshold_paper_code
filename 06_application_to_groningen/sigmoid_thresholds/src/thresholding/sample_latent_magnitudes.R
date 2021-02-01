#' Sample unrounded GPD magnitudes over variable threshold
#'
#' @param n number of samples to draw
#' @param x vector of observed, rounded magnitudes
#' @param sig_u latent scale parameter for exceedances of u
#' @param xi latent shape parameter
#' @param u threshold for which latent parameters are given
#' @param to_nearest level of rounding applied to observations
#' @param v optional vector of latent thresholds at each observation. Y values below this replaced by NA
#'
#' @return matrix of sampled values, one column per entry of X
#' @export
sample_latent_magnitudes <- function(x, sig_u, xi, u, to_nearest, v = NULL){

  x_low <-  x - (to_nearest / 2)
  x_high <- x + (to_nearest / 2)

  q_low <- pgpd(q = x_low, scale = sig_u, shape = xi, mu = u)
  q_high <- pgpd(q = x_high, scale = sig_u, shape = xi, mu = u)

  p <- runif(length(x))

  #matrix of sampled quantiles of Y one column per value of X
  sampled_latent_quantiles <- p * (q_high - q_low) + q_low
  sampled_latent_values <- qgpd(p = sampled_latent_quantiles,
                                shape = xi,
                                scale = sig_u,
                                mu = u)
  if(!is.null(v)){
    sampled_latent_values[sampled_latent_values < v] <- NA
  }
  return(sampled_latent_values)
}



