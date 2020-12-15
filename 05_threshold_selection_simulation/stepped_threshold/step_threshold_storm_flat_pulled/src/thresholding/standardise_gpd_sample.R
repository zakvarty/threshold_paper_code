#' Transform a GPD sample to have standard exponential distribtuion
#'
#' @param y Vector of observed values
#' @param v Modelling threshold for y values
#' @param u threshold for which parameters are given
#' @param sig_u  (scalar) scale parameter above u
#' @param xi (scalar) shape parameter
#'
#' @return vector of standard exponential observations
#'
standardise_gpd_sample <- function (y, v, u, sig_u, xi){
  y[y < v] <- NA
  if(xi < 0) y[y > (u - sig_u / xi)] <- NA

  sig_v <- sig_u + (v - u) * xi
  # values standardised to Exp(1)
  std_exp <- (1 / xi) * log( 1 + xi * ((y - v) / sig_v))
  return(std_exp)
}
