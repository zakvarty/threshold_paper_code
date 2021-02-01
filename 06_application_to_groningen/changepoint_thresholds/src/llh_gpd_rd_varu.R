# dependencies
#source('src/_gpd.R')
#source('src/round_to_nearest.R')

#' Evaluate the (negative) log-likelihood of rounded GPD data with variable threshold
#'
#' @param sigxi GPD parameters at threshold u
#' @param u (Low) latent threshold of continuous GPD for which sigxi is given
#' @param v vector of latent threshold values for each observation
#' @param x Vector of observed (rounded) values
#' @param to_nearest level to which observed values are rounded
#' @param negative logical. Return negative log-likelihood?
#'
#' @return (negative) log-likelihood of rounded GPD data with variable threshold
#' @export
#'
#' @examples
#' # simulate with variable threshold
#' v <- rep(c(1.5,1.1), each = 50)
#' x <- rgpd_rd(n = 100,sig = 1.1,xi = 0.1,mu = v,to_nearest = 0.1)
#' plot(x) ; points(x = seq_along(v), y = v, pch = "-")
#'
#' # evaluate log-likelihood
#' llh_gpd_rd_varu(sigxi = c(1,-0.06),u = 1.0,v = v, x = x)
#' # obervation above UEP
#' llh_gpd_rd_varu(sigxi = c(1,-0.1),u = 1.0,v = v, x = x)
#'
llh_gpd_rd_varu <- function(sigxi, u, v, x, to_nearest = 0.1, negative = FALSE){
  # latent GPD parameters
  sig_u <- sigxi[1]
  xi  <- sigxi[2]
  delta <- to_nearest / 2
  sig_v <- sig_u + xi * (v - u)

  # check x all be above v (adjusted for rounding)
  #lep_fail <- any(x <= v - to_nearest/2)
  lep_fail <- any(x <  v - to_nearest/2)
  if(lep_fail){stop('lower endpoint failure:  !all(x > v - to_nearest / 2).')}

  # check all x below UEP (if it exists, adjusting for rounding)
  latent_uep <- NULL
  uep_fail <- FALSE
  if(xi < 0){
    latent_uep <- u - sig_u / xi
    uep_fail <- max(x) >= (latent_uep + delta)
  }
  if(!is.null(latent_uep) & uep_fail){
    llh <- -10e6
    return((-1)^negative * llh)
  }

  # Make df of parameters for each observation
  x_high <- x + delta
  x_low <- pmax(x - delta, v)
  x_bottom <- x - delta

  # Check all scale params positive
  if(sig_u <= 0){ llh <- -10e6 ; return((-1)^negative * llh)}

  # Calculate log-likelhood
   p_high <- pgpd(q = x_high, shape = xi, scale = sig_v, mu = v)
   p_low  <- pgpd(q = x_low,  shape = xi, scale = sig_v, mu = v)
   p_bottom <- pgpd(q = x_bottom,  shape = xi, scale = sig_v, mu = v)

   weight <- (pgpd(q = x_high,  shape = xi, scale = sig_u, mu = u) -
              pgpd(q = x_low,  shape = xi, scale = sig_u, mu = u)) /
              (pgpd(q = x_high,  shape = xi, scale = sig_u, mu = u) -
              pgpd(q = x_bottom,  shape = xi, scale = sig_u, mu = u))

  llh <- sum(weight * log(p_high - p_low))
  return((-1)^negative * llh)
}
