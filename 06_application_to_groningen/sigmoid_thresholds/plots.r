###
### This is sigmmoid_Groningen/plots.r
###
source('src/get_sigmoid_values.r')

library(dplyr)

cat <- readRDS("./Output/data/cat.RDS")
gron_cat <- readRDS("./Output/data/gron_cat.RDS")


################################################################################
################################################################################
###########         Part 1: Joint parameter selection          #################
################################################################################
################################################################################
## Get the top n_ranks threshold parmameters from each BO initialisation ----
n_inits = 5
n_ranks = 5

top_thresholds <- data.frame(
  init = rep(1:n_inits, each = n_ranks),
  rank = rep(1:n_ranks, n_inits),
  v_1 = rep(NA, n_inits * n_ranks),
  v_2 = rep(NA, n_inits * n_ranks),
  mu = rep(NA, n_inits * n_ranks),
  sigma = rep(NA, n_inits * n_ranks),
  dq1 = rep(NA, n_inits * n_ranks)
)
for( i in 1:5){
  opt_object_loaded <- readRDS(file = paste0('Output/data/opt_object_',i,'.RDS'))
  best_pars_loaded <- opt_object_loaded$scoreSummary %>%
    arrange(-Score) %>%
    mutate(dq1 = 1/ Score) %>%
    select(v_1, v_2,mu, sigma, dq1) %>%
    head(n_ranks)
  top_thresholds[(i-1)*n_ranks + 1:5, 3:7] <- best_pars_loaded
}
rm(best_pars_loaded, opt_object_loaded, i)
head(top_thresholds)



## Plot all the selected thresholds on index and natural timescales ---------
chosen_thresholds <- filter(top_thresholds, rank == 1)

pdf("./Output/plots/BO_selected_thresholds.pdf", width = 7, height = 5)
  par(mar = c(5.1,5.1,1.1,1.1))
  ## Index time-scale
  plot(cat, main = "", pch = 16, col = 'grey70', ylab = 'magnitude', cex.axis = 1.6, cex.lab = 1.6)
  for( i in 1:NROW(chosen_thresholds)){
    chosen_threshold <- get_sigmoid_values(
      tau = cat$index,
      mu = chosen_thresholds$mu[i],
      sigma = chosen_thresholds$sigma[i],
      v_1 = chosen_thresholds$v_1[i],
      v_2 = chosen_thresholds$v_2[i])
    lines(cat$index, chosen_threshold, lwd = 2.5, col = i+1)
  }
  abline(col = 9, h = 1.45, lty = 2, lwd = 2)
  #legend("topleft",
  #       legend = 1:NROW(chosen_thresholds),
  #       col = (1:NROW(chosen_thresholds)) + 1,
  #       pch = 16,
  #       title = "init")

  ## Natural time-scale
  plot(x = gron_cat$date,
       y = gron_cat$mag,
       main = "",
       pch = 16,
       col = 'grey70',
       xlab = 'date',
       ylab = 'magnitude',
       cex.axis = 1.6,
       cex.lab = 1.6
       )
  for( i in 1:NROW(chosen_thresholds)){
    chosen_threshold <- get_sigmoid_values(
      tau = cat$index,
      mu = chosen_thresholds$mu[i],
      sigma = chosen_thresholds$sigma[i],
      v_1 = chosen_thresholds$v_1[i],
      v_2 = chosen_thresholds$v_2[i])
    lines(gron_cat$date, chosen_threshold, lwd = 2.5, col = i+1)
  }
  abline(col = 9, h = 1.45, lty = 2, lwd = 2)
  #legend("topleft",
  #       legend = 1:NROW(chosen_thresholds),
  #       col = (1:NROW(chosen_thresholds)) + 1,
  #       pch = 16,
  #      title = "init")
dev.off()


## Plot the chosen and runner-up thresholds found from each starting point -------
pdf("./Output/plots/BO_good_thresholds_per_initialisation.pdf", width = 7, height = 5)
for(i in 1:n_inits){
  best_thresholds <- filter(top_thresholds, init == i)
  plot(cat, main = paste("Initial set", i), pch = 16, col = 'grey70', xlab = 'date', ylab = 'magnitude')
  for(j in 1:n_ranks){
    chosen_threshold <- get_sigmoid_values(
      tau = cat$index,
      mu = best_thresholds$mu[j],
      sigma = best_thresholds$sigma[j],
      v_1 = best_thresholds$v_1[j],
      v_2 = best_thresholds$v_2[j]
    )
    lines(x = cat$index, y = chosen_threshold, lwd = 3, col = j, lty = 1 + (j>1))
  }
}
dev.off()

## Plot the variablity in metric values for each selected threshold
rerun_values_selected <- readRDS("./Output/data/rerun_values_selected_thresholds.RDS")
pdf("./Output/plots/BO_dq1_reruns_at_selected_thresholds.pdf", width = 7, height = 5)
with(rerun_values_selected,
  { plot(
      x = init,
      y = dq1,
      col = init + 1,
      main = "Groningen BO selection score reruns")
    points(
      x = 0:4,
      colMeans(1/rerun_values),
      col = 1,
      pch = "-",
      cex = 2)
  }
)
dev.off()

## Plot the variability in metric values for the top 5 thresholds of each initialisation
rerun_values_good <- readRDS("./Output/data/rerun_values_good_thresholds.RDS")
pdf("./Output/plots/rerun_values_good_thresholds.pdf", width = 7, height = 5)
with(rerun_values_good,
     plot(x = 10 *init + rank, y = dq1, col = init +1, xlab = "10 init + rank", ylab = "d(q,1)")
)
with(rerun_values_good %>%
       group_by(init, rank) %>%
       summarise(result = mean(dq1)),
     points(x = init * 10 + rank, y = result, pch = "-", cex = 2, )
)
dev.off()


################################################################################
################################################################################
###########         Part 2: fixing end levels           ########################
################################################################################
################################################################################

###
## Get the top n_ranks threshold parmameters from each BO initialisation ----
###
n_inits = 5
n_ranks = 5

top_thresholds <- data.frame(
  init = rep(1:n_inits, each = n_ranks),
  rank = rep(1:n_ranks, n_inits),
  v_1 = rep(1.15, n_inits * n_ranks),
  v_2 = rep(0.76, n_inits * n_ranks),
  mu = rep(NA, n_inits * n_ranks),
  sigma = rep(NA, n_inits * n_ranks),
  dq1 = rep(NA, n_inits * n_ranks)
)
for( i in 1:5){
  opt_object_loaded <- readRDS(file = paste0('Output/data/restricted_selection/opt_object_',i,'.RDS'))
  best_pars_loaded <- opt_object_loaded$scoreSummary %>%
    arrange(-Score) %>%
    mutate(dq1 = 1/ Score) %>%
    select(mu, sigma, dq1) %>%
    head(n_ranks)
  top_thresholds[(i-1)*n_ranks + 1:n_ranks, 5:7] <- best_pars_loaded
}
rm(best_pars_loaded, opt_object_loaded, i)
head(top_thresholds)


###
## Calculate the expected number of exceedances of each top threshold -----
###
source("./src/_gpd.R")
source("./src/llh_gpd_rd_varu.R")
source("./src/mle_gpd_rd_varu.R")
source("./src/bootstrapping.R")
expected_exccedances <- rep(NA, 25)
for(i in 1:NROW(top_thresholds)){
  temp_v <- get_sigmoid_values(
    tau = cat$index,
    mu = top_thresholds$mu[i],
    sigma = top_thresholds$sigma[i],
    v_1 = top_thresholds$v_1[i],
    v_2 = top_thresholds$v_2[i]
  )
  temp_cat <- cat %>% mutate(v = temp_v) %>%  filter(m > (v - 0.1/2 + 1e-08))
  temp_mle <- mle_gpd_rd_varu(
    sigxi = c(1,0),
    u = 0,
    x = temp_cat$m,
    v = temp_cat$v,
    to_nearest = 0.1,
    llh_val = FALSE,
    hessian = FALSE)
  temp_expected_exceedances <- get_n_v_mle(
    x = temp_cat$m,
    v = temp_cat$v,
    to_nearest = 0.1,
    gpd_mle = temp_mle,
    u = 0
  )
  expected_exccedances[i] <- temp_expected_exceedances
}
round(expected_exccedances)


###
## Plot the selected thresholds on Index and natural time scales
###
pdf_path <- "./Output/plots/restricted_selection/BO_selected_thresholds_restricted.pdf"
pdf(pdf_path, width = 7, height = 5)
    par(mar = c(4.1,4.1,1.1,1.1))
    # Plot the selected thresholds- Index time-scale
    plot(cat, main = "", pch = 16, col = 'grey70', ylab = 'magnitude', cex.axis =1.5, cex.lab = 1.5)
    for( i in 1:NROW(top_thresholds)){
      if(top_thresholds$rank[i] == 1){
        temp_threshold <- get_sigmoid_values(
          tau = cat$index,
          mu = top_thresholds$mu[i],
          sigma = top_thresholds$sigma[i],
          v_1 = top_thresholds$v_1[i],
          v_2 = top_thresholds$v_2[i])
        lines(cat$index, temp_threshold, lwd = 3, col = top_thresholds$init[i] + 1)
      }
    }
    abline(col = 9, h = 1.45, lty = 2, lwd = 2)
    #legend("topleft",
    #       legend = paste0(unique(top_thresholds$init), " (",round(expected_exccedances[seq(1,21,by = 5)]),")"),
    #       col = seq_along(unique(top_thresholds$init)) + 1,
    #       pch = 16,
    #       title = "init  & N(A) mle")

    # Plot the selected thresholds - natural timescale
    plot(x = gron_cat$date, gron_cat$mag, xlab = "date", pch = 16, col = 'grey70', ylab = 'magnitude', cex.axis = 1.5, cex.lab = 1.5)
    for( i in 1:NROW(top_thresholds)){
      if(top_thresholds$rank[i] == 1){
        temp_threshold <- get_sigmoid_values(
          tau = cat$index,
          mu = top_thresholds$mu[i],
          sigma = top_thresholds$sigma[i],
          v_1 = top_thresholds$v_1[i],
          v_2 = top_thresholds$v_2[i])
        lines(gron_cat$date, temp_threshold, lwd = 3, col = top_thresholds$init[i] + 1)
      }
    }
    abline(col = 9, h = 1.45, lty = 2, lwd = 2)
    abline(v = as.Date("2014-01-01"), lwd = 2) # G-netork installation starts
    abline(v = as.Date("2014-11-05"), lwd = 2) # First G-network sensors activated
    abline(v = as.Date("2018-01-01"), lwd = 2) # Last G-network sensors activated
    text(x = as.Date("2012-06-01"), y = 3.5, labels = "A", cex = 2) # "Installation \n start"
    text(x = as.Date("2016-01-01"), y = 3.5, labels = "B", cex = 2) # Sensors \n active")
    text(x = as.Date("2019-01-01"), y = 3.5, labels = "C", cex = 2) #All \n active",)

    #abline(v = as.Date("2015-01-11"), lwd = 2) # (approx) change location for turquoise (between events 746 and 747)
    #abline(v = as.Date("2016-11-20"), lwd = 2) # (approx) change loction for reg (event 950.5)

    #legend("topleft",
    #       legend = paste0(unique(top_thresholds$init), " (",round(expected_exccedances[seq(1,21,by = 5)]),")"),
    #       col = seq_along(unique(top_thresholds$init)) + 1,
    #       pch = 16,
    #       title = "init  & N(A) mle")
dev.off()

###
## Plot the rerun metric values for selected and runner-up thresholds -----
###

rerun_values_good <- readRDS('./Output/data/restricted_selection/rerun_values_good_thresholds_1500-10.RDS')

pdf_path <- "./Output/plots/restricted_selection/BO_rerun_values_restricted_1500-10.pdf"
pdf(pdf_path, width = 7, height = 5)
    # plot re-run values
    with(rerun_values_good,
         plot(x = 10 *init + rank, y = dq1, col = init +1, xlab = "10*init + rank", ylab = "d(q,1)", main = "reruns of good thresholds. (sigmoid restricted)")
    )
    # add the mean re-run value of each threshold
    with(rerun_values_good %>%  group_by(init,rank) %>% summarise(dq1_mean = mean(dq1)),
         points(x = init * 10 + rank, y = dq1_mean, pch = '-', cex = 2))

    # add the original run value for each threshold
    with(top_thresholds,
         points(x = init * 10 + rank, y = dq1, col = init + 1, pch = 16))
dev.off()

