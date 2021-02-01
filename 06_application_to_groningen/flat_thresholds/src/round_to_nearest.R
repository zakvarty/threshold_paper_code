round_to_nearest<- function(x, to_nearest = 1, floor = FALSE, ceiling = FALSE ){
  y <- x/to_nearest
  if(floor == FALSE & ceiling == FALSE){
    y <- round(y)
    # Issue: round() rounds xx.5 to nearest even number.
    #decimal_part <- y %% 1
    #ifelse(test = dplyr::near(decimal_part, 0.5) | decimal_part > 0.5,
    #       yes  = y <- ceiling(y),
    #       no   = y <- floor(y))
    # better, butfloating point issues remain very close to boundary.
    # e.g. round_to_nearest(1.44999999, 0.1)  is 1.4
    #      round_to_nearest(1.449999999, 0.1) is 1.5
  }else if(floor == FALSE & ceiling == TRUE){
    y <- ceiling(y)
  }else if(floor == TRUE & ceiling == FALSE){
    y <- floor(y)
  }else if(floor == TRUE & ceiling == TRUE){
    stop('At most one of floor and ceiling may be TRUE')
  }
  y <- y * to_nearest
  return(y)
}
