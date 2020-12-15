#' Calculate mean squared error between Exp(1) cdf and ecdf(x)
#'
#' @param x Sample thought to have Exp(1) margin
#' @param n_eval_pts Number of points at which to evaluate (F_n(x) - F(x))^2
#'
#' @return Approximate mean squared error between Exp(1) cdf and ecdf(x)
#'
exp_pp_mean_squared_distance <- function(x, n_eval_pts){
  eval_quants <- ppoints(n_eval_pts)            # set quantiles of Exp(1) at which to evaluate
  eval_pts <- qexp(p = eval_quants, rate = 1)   # convert to points at which to evaluate
  emp_cdf <- empirical_cdf(x, eval_pts, na.rm = TRUE)  # calculate empirical cdf at these points

  sq_distance <- mean((eval_pts - emp_cdf)^2)
  return(sq_distance)
}
