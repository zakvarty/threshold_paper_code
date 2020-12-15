
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

plot(errors$threshold_2, xlab = 'catalogue', ylab = 'threshold_2 error')
abline(h = 0, col = 2, lty = 2, lwd = 2)
plot(density(errors$threshold_2, na.rm = TRUE), xlab = 'threshold_1 error', main = '')


abline(v = 0, col = 2, lty = 2, lwd = 2)

plot(errors$threshold_2, )
plot(errors$change_location)
plot(density(errors$threshold_2, na.rm = TRUE, bw = 0.05))
plot(density(errors$threshold_1,na.rm = TRUE, bw = 0.04))
plot(density(errors$change_location, na.rm = TRUE, bw = 30))
colMeans(errors, na.rm = TRUE)


#++++++

