
bootstrap_stability_plot_xi <- function(mles, mle_samples, threshold_vector, alpha = 0.05, true_value = NULL, title = NULL, colour = 1,  ...){
  # Set confidence interval bounds
  CI_probs <- c(alpha / 2, 1 - (alpha / 2))

  # Point esitmates
  point_estimates <- map(mles, function(x){x[2]})

  # Calculate CIs
  CIs <- matrix(NA, nrow =  length(threshold_vector), ncol = 2)
  for(i in seq_along(threshold_vector)){
    CIs[i,] <- quantile(mle_samples[[i]][2,], probs = CI_probs)
  }

  # Construct plot
  plot(
    x = rep(thresholds_vec,2),
    y = c(CIs[,1], CIs[,2]),
    xlab = 'threshold value',
    ylab = 'shape parameter',
    bty = 'n',
    type = 'n',
    main = ifelse(
      test = is.null(title),
      yes = paste('point estimate and ',
                  (1 - alpha) * 100,
                  '% confidence intervals'),
      no = title),
    ...
    )
  abline(h = true_value, lty = 2)
  for( i in seq_along(threshold_vector)){
    lines(x = rep(threshold_vector[[i]],2), y = CIs[i,], lwd = 2)
  }
  points(x = thresholds_vec, y = point_estimates, pch = 16, col = colour)
  return(NULL)
}

bootstrap_stability_plot_sig <- function(mles, mle_samples, threshold_vector, alpha = 0.05, true_value = NULL, title = NULL, colour = 1,  ...){
  # Set confidence interval bounds
  CI_probs <- c(alpha / 2, 1 - (alpha / 2))

  # Point esitmates
  point_estimates <- map(mles, function(x){x[1]})

  # Calculate CIs
  CIs <- matrix(NA, nrow =  length(threshold_vector), ncol = 2)
  for(i in seq_along(threshold_vector)){
    CIs[i,] <- quantile(mle_samples[[i]][1,], probs = CI_probs)
  }

  # Construct plot
  plot(
    x = rep(thresholds_vec,2),
    y = c(CIs[,1], CIs[,2]),
    xlab = 'threshold value',
    ylab = 'scale parameter',
    bty = 'n',
    type = 'n',
    main = ifelse(
      test = is.null(title),
      yes = paste('point estimate and ',
                  (1 - alpha) * 100,
                  '% confidence intervals'),
      no = title),
    ...
  )
  abline(h = true_value, lty = 2)
  for( i in seq_along(threshold_vector)){
    lines(x = rep(threshold_vector[[i]],2), y = CIs[i,], lwd = 2)
  }
  points(x = thresholds_vec, y = point_estimates, pch = 16, col = colour)
  return(NULL)
}


bootstrap_stability_plot_RL <- function(mles, mle_samples, threshold_vector, u, v, CRL = 0.9 ,alpha = 0.05, true_value = NULL, title = NULL, colour = 1,  ...){
  # Set confidence interval bounds
  CI_probs <- c(alpha / 2, 1 - (alpha / 2))

  # Point esitmates
  point_estimates <- rep(NA, length(threshold_vector))
  for (i in seq_along(mles)){
    point_estimates[i] <- qgpd(p = CRL, shape = mles[[i]][2], scale = mles[[i]][1] + (v - u)*mles[[i]][2], mu = v)
  }


  # Calculate CIs
  CIs <- matrix(NA, nrow = length(threshold_vector), ncol = 2)
  for (i in seq_along(threshold_vector)){
    return_levels <- qgpd(p = CRL, shape = mle_samples[[i]][2,], scale = mle_samples[[i]][1,] + (v - u)*mle_samples[[i]][2,], mu = v)
    CIs[i,] <- quantile(return_levels, CI_probs)
  }


  # Construct plot
  plot(
    x = rep(thresholds_vec,2),
    y = c(CIs[,1], CIs[,2]),
    xlab = 'threshold value',
    ylab = paste('conditional return level  (v = ',v, ', p = ',CRL, ')'),
    bty = 'n',
    type = 'n',
    main = ifelse(
      test = is.null(title),
      yes = paste('point estimate and ',
                  (1 - alpha) * 100,
                  '% confidence intervals'),
      no = title),
    ...
  )
  abline(h = true_value, lty = 2)
  for( i in seq_along(threshold_vector)){
    lines(x = rep(threshold_vector[[i]],2), y = CIs[i,], lwd = 2)
  }
  points(x = thresholds_vec, y = point_estimates, pch = 16, col = colour)
  invisible(NULL)
}
