# dependencies
#source("src/llh_gpd_rd_varu.R")

#' Maximum likelihood estimation for rounded GPD data with variable threshold values
#'
#' @param sigxi inital parameter values for optimisation (latent gpd parameters at threshold u)
#' @param u threshold of latent gpd at which to return optimal parameter values
#' @param x vector of rounded gpd observations
#' @param v vector of latent gpd thresholds
#' @param to_nearest level to which observed values are rounded
#' @param llh_val logical. Return the log-likelihood value at mle?
#' @param hessian logical. Return numerically estimated hessian at mle?
#' @param ... addtional parameters to be passed to optim
#'
#' @return vector of parameter mles or list containing mles, llh_val and hessian.
#' @export
#'
#' @examples
#' v <- rep(c(1.5,1.1), each = 50)
#' x <- rgpd_rd(n = 100,sig = 1.1,xi = 0.1,mu = v,to_nearest = 0.1)
#' plot(x) ; points(x = seq_along(v), y = v, pch = "-")
#'
#' mle_gpd_rd_varu(sigxi = c(1,-0.05),u = 1,x,v,to_nearest = 0.1)
#'
mle_gpd_rd_varu <- function(sigxi, u, x, v, to_nearest = 0.1, llh_val = TRUE, hessian = FALSE,...){

  sig_init <- sigxi[1]
  xi_init <- sigxi[2]
  sigs_init <-  sig_init + xi_init * (v - u)

  # check valid starting point
  stopifnot(all(sigs_init > 0))
  stopifnot(all(x >= round_to_nearest(v, to_nearest)))
  if(any(x < round_to_nearest(v, to_nearest))){
    stop(
      paste(
        "data x[i] below bin containing v[i] for elements. Occurred at i = ",
        which(x < round_to_nearest(v, to_nearest))
        )
      )
  }
  if((xi_init < 0) & ((u - sig_init/xi_init) <= (max(x) - 0.5*to_nearest))){
    stop("Invalid starting point. Data above upper end point of distribution.")
  }

  # Numerically minimise negative log-likelihood
  temp <- optim(fn = llh_gpd_rd_varu,
                par = sigxi,
                u = u,
                v = v,
                x = x,
                to_nearest = to_nearest,
                negative = TRUE,
                hessian = hessian,
                ...)

  #format output
  if(!llh_val & !hessian){out <-  temp$par} else {out <- list(params = temp$par)}
  if(llh_val) out$loglik <- -temp$value
  if(hessian) out$hessian <- -temp$hessian

  return(out)
}
