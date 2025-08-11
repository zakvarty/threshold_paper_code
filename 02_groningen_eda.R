# Groningen catalogue EDA:

# - Filter full catalogue to date range [1995-01-01 to 2025-01-01)
# - Plot on natural time scale with loess mean
# - Plot on index time scale with loess mean
# - Plot QQ and PP plots for magnitudes exceeding 0.0, 0.76 and 1.45 M_l
# - Plot QQ and PP plots using uniformly jittered data


## SOURCE ----------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(here)
library(readr)
library(threshold)
library(lubridate)

SEED <- 1234
OUTPUT_PATH <- here("00_outputs", "02_motivation_and_model")

set.seed(SEED)

## LOAD DATA -------------------------------------------------------------------

# earthquake catalogue
cat <- threshold::groningen_earthquakes
cat <- cat %>%
  mutate(index = order(decimal_date)) %>%
  mutate(date_time = lubridate::as_datetime(stringr::str_c(date, time))) %>%
  filter(date >= as_date("1995-01-01")) %>%
  filter(date < as_date("2025-01-01")) %>%
  select(mag, index, date, date_time)

# network development dates
network_dates <- threshold::network_dates

date_to_index <- function(date, catalog = cat){sum(catalog$date <= date)}
network_dates$index <- vapply(network_dates$date, date_to_index, FUN.VALUE = 1)
network_dates$label_index <- network_dates$index + c(-40, 40, 40)

## PLOT CATALOGUES -------------------------------------------------------------

p1_natural_timescale <- ggplot(cat, aes(x = date_time, y = mag)) +
  geom_point(colour = rgb(0.7,0.7,0.7,0.7), size = 1.2) +
  geom_smooth(method = 'loess', col = 'black', fill = rgb(0,0,0,0.5)) +
  ylab('magnitude') +
  xlab('event time') +
  theme_minimal(base_size = 22) +
  geom_vline(data = network_dates, mapping = aes(xintercept = date)) +
  geom_text(data = network_dates,
            mapping = aes(x = label_date, label = code),
            y = -0.25,
            size = 7)

p2_index_timescale <- ggplot(cat, aes(x = index, y = mag)) +
  geom_point(colour = rgb(0.7,0.7,0.7,0.7), size = 1.2) +
  geom_smooth(method = 'loess', col = 'black', fill = rgb(0,0,0,0.5)) +
  ylab('magnitude') +
  xlab('event index') +
  theme_minimal(base_size = 22) +
  geom_vline(data = network_dates, mapping = aes(xintercept = index)) +
  geom_text(data = network_dates,
            mapping = aes(x = label_index, label = code),
            y = -0.25,
            size = 7)


## SAVE CATALGOUE PLOTS --------------------------------------------------------

file_name <- "groningen_catalogue_natural.png"
path <- here(OUTPUT_PATH, file_name)
ggsave(p1_natural_timescale, filename = path, width = 8, height = 5)
cat(file_name, "complete \n")

file_name <- "groningen_catalogue_index.png"
path <- here(OUTPUT_PATH, file_name)
ggsave(p2_index_timescale, filename = path, width = 8, height = 5)
cat(file_name, "complete \n")


## QQ, PP and ACF plots --------------------------------------------------------


thresholds_str <- c("0.00", "0.76", "1.45")
thresholds_str_minimal <- stringr::str_remove(thresholds_str, "\\.")
thresholds <- as.numeric(thresholds_str)

qqs <- list()
qqs_jittered <- list()

for (t in seq_along(thresholds)) {

  base_file_name <- stringr::str_c("groningen_", thresholds_str_minimal[t])

  # unjittered QQ plot ------
  file_name <- stringr::str_c(base_file_name, "_qq.pdf")
  path <- here(OUTPUT_PATH, file_name)

  pdf(file = path, width = 4, height = 4)
  par(mar = c(4.1, 4.1, 1.1, 1.1))
  qqs[[t]] <- qq_gpd(x = cat$mag, u = thresholds[t], seed = SEED)
  dev.off()
  cat(file_name, "complete \n")

  # unjittered PP plot -------
  file_name <- stringr::str_c(base_file_name, "_pp.pdf")
  path <- here(OUTPUT_PATH, file_name)

  pdf(file = path, width = 4, height = 4)
  par(mar = c(4.1, 4.1, 1.1, 1.1))
  pp_gpd(x = cat$mag, u = thresholds[t], seed = SEED)
  dev.off()
  cat(file_name, "complete \n")

  # jittered QQ plot ------
  file_name <- stringr::str_c(base_file_name, "_qq_jittered.pdf")
  path <- here(OUTPUT_PATH, file_name)

  pdf(file = path, width = 4, height = 4)
  par(mar = c(4.1, 4.1, 1.1, 1.1))
  qqs_jittered[[t]] <- qq_gpd(x = cat$mag, u = thresholds[t], jitter = 0.1, seed = SEED)
  dev.off()
  cat(file_name, "complete \n")

  # jittered PP plot ------
  file_name <- stringr::str_c(base_file_name, "_pp_jittered.pdf")
  path <- here(OUTPUT_PATH, file_name)

  pdf(file = path, width = 4, height = 4)
  par(mar = c(4.1, 4.1, 1.1, 1.1))
  pp_gpd(x = cat$mag, u = thresholds[t], jitter = 0.1, seed = SEED)
  dev.off()
  cat(file_name, "complete \n")

  # acf plot ------
  file_name <- stringr::str_c(base_file_name, "_acf.pdf")
  path <- here(OUTPUT_PATH, file_name)
  excesses <- cat$mag[cat$mag > thresholds[t]] - thresholds[t]

  pdf(file = path, width = 4, height = 4)
  par(mar = c(4.1, 4.1, 1.1, 1.1))
  acf(x = excesses, bty = "n", lag.max = 30)
  dev.off()
  cat(file_name, "complete \n")

  # acf plot jittered ------
  n_excesses <- length(excesses)
  excesses_jittered <- excesses + runif(n_excesses, min = -0.05, max = 0.05)

  file_name <- stringr::str_c(base_file_name, "_acf_jittered.pdf")
  path <- here(OUTPUT_PATH, file_name)

  pdf(file = path, width = 4, height = 4)
  par(mar = c(4.1, 4.1, 1.1, 1.1))
  acf(x = excesses_jittered, bty = "n", lag.max = 30)
  dev.off()
  cat(file_name, "complete \n")
}

cat("Done.")

names(qqs) <- thresholds_str_minimal
names(qqs_jittered) <- thresholds_str

## Record summaries of model fits above trial thresholds -----------------------


parameter_names <- c("sigma_u", "xi", "sigma_0")

extract_mles <- function(qq_obj, u){
  c(
    "sigma" = qq_obj$mle_sigma,
    "xi" = qq_obj$mle_xi,
    "sigma_0" = qq_obj$mle_sigma - u * qq_obj$mle_xi)
}

extract_CI <- function(qq_obj, u, prob){
  bounds <- c(
    quantile(qq_obj$boot_sigmas, probs = prob),
    quantile(qq_obj$boot_xis, probs = prob),
    quantile(qq_obj$boot_sigmas - u * qq_obj$boot_xis, probs = prob)
  )
  names(bounds) <- parameter_names
  bounds
}

mles <- c(
  extract_mles(qqs[[1]], u = thresholds[[1]]),
  extract_mles(qqs[[2]], u = thresholds[[2]]),
  extract_mles(qqs[[3]], u = thresholds[[3]]))

ci_lower <- c(
  extract_CI(qqs[[1]], u = thresholds[[1]], prob = 0.025),
  extract_CI(qqs[[2]], u = thresholds[[2]], prob = 0.025),
  extract_CI(qqs[[3]], u = thresholds[[3]], prob = 0.025))

ci_upper <- c(
  extract_CI(qqs[[1]], u = thresholds[[1]], prob = 0.975),
  extract_CI(qqs[[2]], u = thresholds[[2]], prob = 0.975),
  extract_CI(qqs[[3]], u = thresholds[[3]], prob = 0.975))


fit_summaries <- tibble(
  u = rep(thresholds, each = 3),
  parameter = rep(parameter_names, times = length(thresholds)),
  mle = mles,
  ci_lower = ci_lower,
  ci_upper = ci_upper
)

readr::write_csv(x = fit_summaries,file = here(OUTPUT_PATH, "fit_summaries.csv"))


mles_jittered <- c(
  extract_mles(qqs_jittered[[1]], u = thresholds[[1]]),
  extract_mles(qqs_jittered[[2]], u = thresholds[[2]]),
  extract_mles(qqs_jittered[[3]], u = thresholds[[3]]))

ci_lower_jittered <- c(
  extract_CI(qqs_jittered[[1]], u = thresholds[[1]], prob = 0.025),
  extract_CI(qqs_jittered[[2]], u = thresholds[[2]], prob = 0.025),
  extract_CI(qqs_jittered[[3]], u = thresholds[[3]], prob = 0.025))

ci_upper_jittered <- c(
  extract_CI(qqs_jittered[[1]], u = thresholds[[1]], prob = 0.975),
  extract_CI(qqs_jittered[[2]], u = thresholds[[2]], prob = 0.975),
  extract_CI(qqs_jittered[[3]], u = thresholds[[3]], prob = 0.975))


fit_summaries_jittered <- tibble(
  u = rep(thresholds, each = 3),
  parameter = rep(parameter_names, times = length(thresholds)),
  mle = mles_jittered,
  ci_lower = ci_lower_jittered,
  ci_upper = ci_upper_jittered
)

readr::write_csv(x = fit_summaries_jittered,file = here(OUTPUT_PATH, "fit_summaries_jittered.csv"))

