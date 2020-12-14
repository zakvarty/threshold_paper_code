#' Calculate squared bias and variance components for each parameter from a vector of maximum likelihood estimates
#'
#' @param mles data.frame where each row represents a maximum likelihood estimate of a parameter vector
#' @param pars_true numeric vector giving true parameter values
#' @param par_names character vector giving parameter names
#' @param plot_it logical should error compontents be plotted in a bar chart?
#' @param ... additional arguments to bar chart
#'
#' @return silently returns data.frame of error components
#' @export
#'
split_MSE<- function(mles, pars_true, par_names, plot_it = TRUE, ...){

  n_par <- length(pars_true)
  errors <- mles
  biases <- rep(NA_real_, length(pars_true))
  variances <- rep(NA_real_, length(pars_true))

  for(par in 1: ncol(mles)){
    errors[,par] <- mles[,par] - pars_true[par]
    biases[par] <- mean(errors[,par])
    variances[par] <- var(mles[,par])
  }

  sq_biases <- biases^2

  out <- data.frame(parameter = par_names, sq_bias = sq_biases, variance = variances)

  if(plot_it){
  names_bias <- paste0('bias^2(', par_names, ')')
  names_var <- paste0('var(', par_names, ')')

  barplot(height = c(out$sq_bias, out$variance),
          names.arg = c(names_bias, names_var),...)
  }

  invisible(out)
}
