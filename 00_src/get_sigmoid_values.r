get_sigmoid_values <- function(tau, mu, sigma, v_1, v_2){
  v_1 + pnorm((tau - mu)/ sigma) * (v_2 - v_1)
}
