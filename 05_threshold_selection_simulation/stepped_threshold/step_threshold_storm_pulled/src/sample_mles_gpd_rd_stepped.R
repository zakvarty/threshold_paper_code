

sample_mles_gpd_rd_stepped <- function(mle, u, step_lengths, step_values, expected_exceedances, to_nearest, n_sim, verbose) {

  sig_u_mle <- mle[1]
  xi_mle  <- mle[2]
  delta <- to_nearest / 2

  # Calculate probability of being allocated to each step
  step_survivor_value <- (1 - pgpd(q = step_values,
                           shape = xi_mle,
                           scale = sig_u_mle,
                           mu = u))
  step_probs <- step_lengths * step_survivor_value

  # Simulate n_sim alternative catalogues
  sim_threshs <- vector(mode = "list", length = n_sim)
  sim_mags <- vector(mode = "list", length = n_sim)

  for(i in 1:n_sim){
    # sample number of exceedances of each unique threshold
    ## lambda given x
    sim_Lambda <- rgamma(n = 1, shape = 1 + expected_exceedances, rate = 1)
    sim_count <- rpois(n = 1, lambda = sim_Lambda)
    sim_v <- sample(step_values, size = sim_count, replace = TRUE, prob = step_probs)
    sim_threshs[[i]] <- sim_v

    # calculate scale parameters and simulate continuous data
    sim_sig_vs <- sig_u_mle + (sim_v - u) * xi_mle
    sim_cts <- rgpd(n = length(sim_v), mu = sim_v, scale = sim_sig_vs, shape = xi_mle)

    # round continuous values
    sim_rnd <- round_to_nearest(sim_cts, to_nearest)
    sim_mags[[i]] <- sim_rnd
  }

  # Find mle of each bootstrap data set

  if(verbose){
    message('starting map...')
    t <- Sys.time()
  }

   sim_mles <- purrr::pmap(
    .l = list(bootstrap_x = sim_mags, bootstrap_v = sim_threshs),
    .f = bootstrap_mle_gpd_rd_varu,
    sigxi = c(1,0),
    u = u,
    to_nearest = to_nearest,
    llh_val = FALSE
  )

  if(verbose){
    t <- Sys.time() - t
    message('ended map', round(t,2), ".")
  }

  sim_mles <- do.call(rbind, sim_mles)
  t(sim_mles)
}




