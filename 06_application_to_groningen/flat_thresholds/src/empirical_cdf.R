empirical_cdf <- function(x, eval_pts, na.rm = FALSE){
  # function to calculate ecdf at single evaluation point
  ecdf_value <- function(data,eval_pt){
    mean(data <= eval_pt, na.rm =  na.rm)
  }
  # apply funciton at each evaluation point
  purrr::map_dbl(.x = eval_pts, .f = ecdf_value, data = x)
}



empirical_pmf <- function(x, eval_pts, na.rm = FALSE){
  #outer(X = 1:10, Y = 3:5, FUN = dplyr::near)
  near <- outer(X = x, Y = eval_pts, FUN = dplyr::near)
  prob <- apply(near, 2, mean, na.rm = na.rm)
  return(prob)
}

# # Examples
# empirical_cdf(x = c(1,1,2,3,4,4,5,6,7,8), eval_pts = c(1,2,3))
# empirical_pmf(x = c(1,1,2,3,4,4,5,6,7,8), eval_pts = c(1,1.5,2))
