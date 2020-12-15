
n_jobs <- 100
v_values <- c(0.42,0.42)

selected <- data.frame(
  jobid = rep(NA_integer_, n_jobs),
  v_1_pp_EWMSE = rep(NA_real_, n_jobs),
  v_2_pp_EWMSE = rep(NA_real_, n_jobs),
  v_1_pp_EWMAE = rep(NA_real_, n_jobs),
  v_2_pp_EWMAE = rep(NA_real_, n_jobs),
  v_1_qq_EMSE  = rep(NA_real_, n_jobs),
  v_2_qq_EMSE  = rep(NA_real_, n_jobs),
  v_1_qq_EMAE  = rep(NA_real_, n_jobs),
  v_2_qq_EMAE  = rep(NA_real_, n_jobs)
)

for(jobid in seq(1,n_jobs)){
# Load the metrics for that jobid
rds_path <- paste0("./Desktop/storm_pulls/step_threshold_storm_flat/Output/data/metrics/metrics_",jobid,".RDS")
metrics <- readRDS(rds_path)
print(jobid)
# extract the threshold pair selected by each metric
selected_this_run <- data.frame(
  jobid = jobid,
  v_1_pp_EWMSE = metrics$v_1[which.min(metrics$pp_EWMSE)],
  v_2_pp_EWMSE = metrics$v_2[which.min(metrics$pp_EWMSE)],
  v_1_pp_EWMAE = metrics$v_1[which.min(metrics$pp_EWMAE)],
  v_2_pp_EWMAE = metrics$v_2[which.min(metrics$pp_EWMAE)],
  v_1_qq_EMSE  = metrics$v_1[which.min(metrics$qq_EMSE)],
  v_2_qq_EMSE  = metrics$v_2[which.min(metrics$qq_EMSE)],
  v_1_qq_EMAE  = metrics$v_1[which.min(metrics$qq_EMAE)],
  v_2_qq_EMAE  = metrics$v_2[which.min(metrics$qq_EMAE)]
)

# add the selected thresholds to a central data.frame
selected[jobid,] <- selected_this_run[1,]
}
rm(selected_this_run)

head(selected,20)
dplyr::summarise_all(.tbl = selected, .funs = min)
dplyr::summarise_all(.tbl = selected, .funs = max)

epsilon <- 0.003
pdf_path <- "./Desktop/storm_pulls/step_threshold_storm_flat/Output/plots/selected.pdf"
pdf(file = pdf_path, width = 6, height = 6)


par(mfrow =c(2,2))
plot(
  x = selected$v_1_pp_EWMAE + runif(n = n_jobs, -epsilon, epsilon),
  y = selected$v_2_pp_EWMAE + runif(n = n_jobs, -epsilon, epsilon),
  pch = 16,
  col = rgb(0,0,0,0.2),
  asp = 1,
  xlab = 'v_1',
  ylab = 'v_2',
  main = 'd(p,1)',
  xlim = range(metrics$v_1),
  ylim = range(metrics$v_2))
points(x = v_values[1],, y = v_values[2], pch = 16, col= 2, cex = 0.6)

plot(
  x = selected$v_1_pp_EWMSE + runif(n = n_jobs, -epsilon, epsilon),
  y = selected$v_2_pp_EWMSE + runif(n = n_jobs, -epsilon, epsilon),
  pch = 16,
  col = rgb(0,0,0,0.2),
  asp = 1,
  xlab = 'v_1',
  ylab = 'v_2',
  main = 'd(p,2)',
  xlim = range(metrics$v_1),
  ylim = range(metrics$v_2))
points(x = v_values[1],, y = v_values[2], pch = 16, col= 2, cex = 0.6)

plot(
  x = selected$v_1_qq_EMAE + runif(n = n_jobs, -epsilon, epsilon),
  y = selected$v_2_qq_EMAE + runif(n = n_jobs, -epsilon, epsilon),
  pch = 16,
  col = rgb(0,0,0,0.2),
  xlab = 'v_1',
  ylab = 'v_2',
  main = 'd(q,1)',
  xlim = range(metrics$v_1),
  ylim = range(metrics$v_2))
points(x = v_values[1],, y = v_values[2], pch = 16, col= 2, cex = 0.6)

plot(
  x = selected$v_1_qq_EMSE + runif(n = n_jobs, -epsilon, epsilon),
  y = selected$v_2_qq_EMSE + runif(n = n_jobs, -epsilon, epsilon),
  pch = 16,
  col = rgb(0,0,0,0.2),
  xlab = 'v_1',
  ylab = 'v_2',
  main = 'd(q,2)',
  xlim = range(metrics$v_1),
  ylim = range(metrics$v_2))
points(x = v_values[1],, y = v_values[2], pch = 16, col= 2, cex = 0.6)
dev.off()
