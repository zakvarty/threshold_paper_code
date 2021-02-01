return_levels <- function(sigxi, u, v, return_periods_v){
  sig_u <- sigxi[1]
  xi <- sigxi[2]
  sig_v <- sig_u + (v - u) * xi
  qgpd(p = (1-1/return_periods_v), shape = xi, scale = sig_v, mu = v)
}


plot_return_levels <- function(mle_and_bootstraps, u, v, v_return_periods = NULL, alpha = 0.05, plot = TRUE, add = FALSE,...){

  if(is.null(v_return_periods)){
   v_return_periods =  1/(0.1^(seq(0,3,length.out = 101)))
  }

  point_rls <- return_levels(
    sigxi = mle_and_bootstraps$mle,
    u = u,
    v = v,
    return_periods_v = v_return_periods)

  boot_rls <-  apply(
    mle_and_bootstraps$boot_mles,
    MARGIN = 1,
    FUN = return_levels,
    u = u,
    v = v,
    return_periods_v = v_return_periods)

  CI_rls <- apply(boot_rls, MARGIN = 1, quantile, probs = c(alpha / 2, 1 - alpha /2))

  out <- list(
    return_levels = data.frame(
      return_period = v_return_periods,
      return_level_mle = point_rls),
      CI_low = CI_rls[1,],
      CI_high = CI_rls[2,])

    if(plot & !add){
      plot(x = rep(v_return_periods,3), y = c(point_rls, CI_rls[1,], CI_rls[2,]), ylab = "return level", xlab = 'return period', ..., type = 'n')
      lines(x = v_return_periods,y =  point_rls, lty = 1, ...)
      lines(x = v_return_periods,y =  CI_rls[1,], lty = 2, ...)
      lines(v_return_periods, CI_rls[2,], lty = 2, ...)
    } else if (plot & add) {
      lines(x = v_return_periods,y =  point_rls, lty = 1, ...)
      lines(x = v_return_periods,y =  CI_rls[1,], lty = 2, ...)
      lines(v_return_periods, CI_rls[2,], lty = 2, ...)
    }
  invisible(out)
}

#plot_return_levels(mle_and_bootstraps = conservative_mle, u = -1, v = 1.45, col = 3, log = 'x')
#plot_return_levels(selected_MLE, u = -1, v = 1.45, col = 2, add = TRUE)
