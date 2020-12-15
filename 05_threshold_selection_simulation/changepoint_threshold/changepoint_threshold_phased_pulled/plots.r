
# Get the indices of sucessful runs

downloaded_selected <- as.numeric(stringr::str_sub(
  string = list.files("./Output/data/selected/"),start = 31, end = -5))

downloaded_catalogue_sizes <- as.numeric(stringr::str_sub(
  string = list.files("./Output/data/catalogue_sizes/"),start = 16, end = -5))

downloaded_opt_objects <- as.numeric(stringr::str_sub(
  string = list.files("./Output/data/opt_objects/"),start = 12, end = -5))

downloaded_runs <- intersect(downloaded_selected, downloaded_catalogue_sizes)

n_downloaded <- length(downloaded_runs)

# Load true parameters
load(file = "./Output/data/true_parameter_values.RData")

# Storage for catalogue sizes and true change locations for each catalogue
catalogue_sizes <- data.frame(
  jobid = rep(NA_real_, n_downloaded),
  n_total = rep(NA_real_, n_downloaded),
  n_1 = rep(NA_real_, n_downloaded),
  n_2 = rep(NA_real_, n_downloaded))

# Storage for selected thresholds and change location for each catalogue
selected <- data.frame(
  threshold_1 = rep(NA_real_, n_downloaded),
  threshold_2 = rep(NA_real_, n_downloaded),
  change_location = rep(NA_real_, n_downloaded)
)


# Load true change locations and selected parameters into memory
for(cat_id in downloaded_runs){
  print(cat_id)

  # Load the metrics for that cat_id
  rds_path <- paste0("./Output/data/catalogue_sizes/catalogue_size_",cat_id,".RDS")
  temp <- readRDS(rds_path)
  catalogue_sizes[cat_id,] <- temp

  # load the selected threshold for that cat_id
  rds_path <- paste0("./Output/data/selected/selected_threshold_parameters_",cat_id,".RDS")
  temp <- readRDS(rds_path)
  selected[cat_id,] <- temp
}

head(selected)
plot(selected)

plot(catalogue_sizes$n_total)

errors <- data.frame(
  threshold_1 = selected$threshold_1 - v[1],
  threshold_2 = selected$threshold_2 - v[2],
  change_location = selected$change_location - catalogue_sizes$n_1
)

# Make some plots of selection errors
pairs(errors, pch = 16,  col = rgb(0,0,0,0.5),lower.panel = NULL, cex.axis = 1.5)

pdf("./Output/plots/selection_errors_changepoint_phased.pdf", width = 5, height = 5)
GGally::ggpairs(errors,upper = list(continuous = "blank")) +
  ggplot2::ggtitle("changepoint threshold selection errors, hard censoring.") +
  ggplot2::theme_bw()
dev.off()

pdf("./Output/plots/location_error_v_step_size_changepoint_phased.pdf", width = 7, height =5)
plot( x = errors$change_location ,
      y = selected$threshold_1 - selected$threshold_2,
      xlab = 'error in change location',
      ylab = 'chosen step size',
      main = 'changepoint threshold, phased censoring')
abline(h = v[1] - v[2], col = 2, lty = 2)
abline(v = 0, col = 2, lty = 2)
dev.off()

pdf("./Output/plots/selection_errors_II_changepoint_phased.pdf", width = 8, height = 6)
par(mfrow=c(2,3))
plot(errors$threshold_1, xlab = 'catalogue', ylab = 'threshold_1 error')
abline(h = 0, col = 2, lty = 2, lwd = 2)
plot(errors$threshold_2, xlab = 'catalogue', ylab = 'threshold_2 error')
abline(h = 0, col = 2, lty = 2, lwd = 2)
plot(errors$change_location, xlab = 'catalogue', ylab = 'change_location error')
abline(h = 0, col = 2, lty = 2, lwd = 2)

plot(density(errors$threshold_1, na.rm = TRUE), xlab = 'threshold_1 error', main = '')
abline(v = 0, col = 2, lty = 2, lwd = 2)
plot(density(errors$threshold_2, na.rm = TRUE), xlab = 'threshold_2 error', main = '')
abline(v = 0, col = 2, lty = 2, lwd = 2)
plot(density(errors$change_location, na.rm = TRUE), xlab = 'change_location error', main = '')
abline(v = 0, col = 2, lty = 2, lwd = 2)
dev.off()
#plot(density(errors$threshold_2, na.rm = TRUE, bw = 0.04), xlab = 'threshold_1 error', main = '')
#abline(v = 0, col = 2, lty = 2, lwd = 2)
#plot(density(errors$threshold_2, na.rm = TRUE, bw = 0.04), xlab = 'threshold_1 error', main = '')
#abline(v = 0, col = 2, lty = 2, lwd = 2)

par(mfrow = c(1,1))
plot(x = errors$threshold_2,y =  errors$change_location )

errors_trimmed <- errors[downloaded_runs, ]
change_location_errors_1 <- errors_trimmed%>%
  filter(between(threshold_2,-0.1,0)) %>%
  pull(threshold_2)

change_location_errors_2 <- errors_trimmed%>%
  filter(between(threshold_2,0,0.1)) %>%
  pull(threshold_2)

ttest <- t.test(x = change_location_errors_1, y = change_location_errors_2)
CIs <- c(mean(change_location_errors_1),
mean(change_location_errors_2),
quantile(change_location_errors_1, prob = 0.025),
quantile(change_location_errors_2, prob = 0.025),
quantile(change_location_errors_1, prob = 0.975),
quantile(change_location_errors_2, prob = 0.975))
plot(x = rep(c(1,2), 3),
     y = CIs ,
     pch = "-",
     cex = 2,
     ylab = "change_location_error",
     xlab = "threshold 2 error mode",
     main = paste0("CI on difference = (", round(ttest$conf.int[1],4),",",round(ttest$conf.int[2],4),")" ))

###
# Integrated absolute errors
###

changepoint_IAE <- function(selected_v1, selected_v2, selected_tau, true_v1, true_v2, true_tau, tau_max){

  if(any(is.na(c(selected_v1, selected_v2, selected_tau, true_v1, true_v2, true_tau, tau_max)))){
    return(NA)
  }

  bin_width <- tau_max / 5001
  eval_pts <- seq(0,tau_max,length.out = 5001)
  selected_threshold_values <- ifelse(eval_pts < selected_tau, selected_v1, selected_v2)
  true_threshold_values <- ifelse(eval_pts < true_tau, true_v1, true_v2)
  IAE <- sum(bin_width * abs(selected_threshold_values - true_threshold_values))
  IAE
}

IAEs <- purrr::pmap(
  .l = list(
    selected_v1 = selected$threshold_1,
    selected_v2 = selected$threshold_2,
    selected_tau = selected$change_location,
    true_tau = catalogue_sizes$n_1,
    tau_max = catalogue_sizes$n_total
  ),
  .f = changepoint_IAE,
  true_v1 = v[1],
  true_v2 = v[2]
)

pdf_path <- "./Output/plots/threshold_IAE.pdf"
pdf(pdf_path, width = 7, height = 5)

plot(y = selected$threshold_1 - selected$threshold_2,
     x = errors$change_location,
     xlab = "change location error",
     ylab = "selected step size")

plot(y = abs(selected$threshold_1 - selected$threshold_2) * abs(errors$change_location),
     x = errors$change_location,
     xlab = "change location error",
     ylab = "middle area")
plot(y = unlist(IAEs),
     x = errors$change_location,
     xlab = "change_location error",
     ylab = "threshold IAE",
     main = "Integrated absolute error in threshold selection",
     pch = 16,
     ylim = c(0,500))
dev.off()

plot(errors$change_location)
plot(density(errors$threshold_2, na.rm = TRUE, bw = 0.05))
plot(density(errors$threshold_1,na.rm = TRUE, bw = 0.04))
plot(density(errors$change_location, na.rm = TRUE, bw = 30))
colMeans(errors, na.rm = TRUE)


#++++++

