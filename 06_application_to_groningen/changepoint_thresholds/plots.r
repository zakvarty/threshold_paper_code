
# Get the indices of sucessful runs
downloaded_runs <- as.numeric(stringr::str_sub(
  string = list.files("./Output/data/selected/"),start = 31, end = -5))
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
  threhsold_2 = selected$threshold_2 - v[2],
  change_location = selected$change_location - catalogue_sizes$n_1
)

pairs(errors, pch = 16,  col = rgb(0,0,0,0.5),lower.panel = NULL)

pdf("./Output/plots/selection_errors.pdf", width = 5, height = 5)
GGally::ggpairs(errors,upper = list(continuous = "blank"))+ ggplot2::theme_bw()
dev.off()
GGally::ggally_text()
#++++++

