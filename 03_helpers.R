#helper functions for motivating example

compare_MSEs <- function(MSEs, bar_colors = NULL, bar_names = NULL, legend_position = 'topleft', legend_text = letters[seq_along(MSEs)], legend_cex = 1, ...){

  # Function to plot many MSE decompositions side-by-side
  heights <- NULL
  if (is.null(bar_names)) {
    bar_names <- c('bais^2(sigma)', 'bais^2(xi)', 'var(sigma)','var(xi)')
    }
  if (is.null(bar_colors)) { bar_colors <- grey(ppoints(length(MSEs))) }

  for (i in seq_along(MSEs)) {
    heights_i <- c(MSEs[[i]]$sq_bias, MSEs[[i]]$variance)
    heights <- rbind(heights,heights_i)
  }

  barplot(height = heights, beside = TRUE, names.arg = bar_names, col = bar_colors, ...)
  legend(legend_position, legend = legend_text, col = bar_colors, pch = 15, bty = 'n', cex = legend_cex)

  heights_df <- as.data.frame(heights, row.names = legend_text)
  colnames(heights_df) <- bar_names
  invisible(heights_df)
}

# Calculate conditional return levels for excesses of V for GPD models fitted
# to excesses of u.
mle_return_levels <- function(mle_obj, shift, v){
  out <- purrr::pmap(
    .l = list(scale = mle_obj$sig_u, shape = mle_obj$xi),
    .f = return_level,
    period = return_periods,
    shift = shift,
    v = v)

  # reformatting output
  out <- do.call(rbind, out)
  colnames(out) <- stringr::str_c(return_periods)

  as_numeric_factor <- function(x) {as.numeric(levels(x))[x]}

  out %>% as.data.frame() %>%
    tibble::rownames_to_column("boot_id") %>%
    tidyr::pivot_longer(-boot_id,
                        names_to = "return_period",
                        values_to = "return_level") %>%
    mutate(
      boot_id = as.numeric(boot_id),
      return_period = as.numeric(return_period))
}



return_level_confidence_interval <- function(rl_obj, alpha = 0.05){
  as_numeric_factor <- function(x) {as.numeric(levels(x))[x]}

  rl_obj %>%
    mutate(return_period = as.factor(return_period)) %>%
    summarise(
      #return_period = return_period,
      lower_limit = unname(quantile(return_level, alpha / 2)),
      median = unname(quantile(return_level, 0.5)),
      upper_limit = unname(quantile(return_level, 1 - alpha / 2)),
      .by = return_period) %>%
    mutate(return_period = as_numeric_factor(return_period))
}

CI_lines <- function(CI_obj, ...){
  lines(x = CI_obj$return_period, y = CI_obj$median, ...)
  lines(x = CI_obj$return_period, y = CI_obj$lower_limit, ...)
  lines(x = CI_obj$return_period, y = CI_obj$upper_limit, ...)
}
