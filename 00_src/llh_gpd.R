#' log-likelihood of generalised Pareto distribtuion
#'
#' @param sigxi scale and shape parameters (scale must be positive)
#' @param u threshold
#' @param x data values
#' @param negative logical. Should negative log-likelihood be returned instead?
#'
#' @return (negative) log-probability of `x` values exceeding `u` at given parameter values `sigxi`.
#' @export
#'
#' @examples
#' y <- rnorm(1000)
#' llh_gpd(sigxi = c(1,0), u = 2, x = y)
llh_gpd <- function(sigxi, u, x, negative = FALSE){

  sig_u <- sigxi[1]
  xi  <- sigxi[2]
  stopifnot(all(x > u))

  # check parameters are valid
  if (sig_u <= 0) {llh <- -10e6 ; return((-1)^negative * llh)}

  if ((xi < 0) & any(x >= u - sig_u / xi)) {
    llh <- -10e6
    return((-1)^negative * llh)
  }

  # calculate log-likelihood
  llh_terms <- dgpd(x = x, mu = u, scale = sig_u, shape = xi, log = TRUE)
  llh <- sum(llh_terms)
  return((-1)^negative * llh)
}
