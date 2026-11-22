demo_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    dirname(normalizePath(sub("^--file=", "", file_arg)))
  } else {
    getwd()
  }
})

library(ggplot2)
library(ggpubr)
library(zoo)
source(file.path(demo_dir, "R", "src", "main.R"))
source(file.path(demo_dir, "R", "helpers", "aemet_models.R"))
source(file.path(demo_dir, "R", "helpers", "plot_bands.R"))
source(file.path(demo_dir, "R", "helpers", "plot_surface_3d.R"))

## Figures for the summer-temperature application. The generic band and
## panel helpers live in R/helpers/plot_bands.R; the labels and layout below
## are specific to this dataset.
plot_latitude_summer <- function(fit, figure_dir = "output/aemet") {
  titles <- c("Intercept", "Year > 2008", expression("CO"[2] * " Emission"),
              "Wind Speed", "Max Precipitation")
  n_day <- nrow(fit$betaHat[[1]])
  n_latitude <- ncol(fit$betaHat[[1]])
  day_breaks <- seq(1, n_day, length.out = 4)
  day_labels <- c("Jun", "Jul", "Aug", "Sep")
  latitude_breaks <- seq(1, n_latitude, length.out = 4)
  latitude_labels <- c(35, 38, 41, 44)
  shown <- c(1, 2, 3, 5)
  save_panels(
    bivariate_panels(fit, titles, x = "Day", y = "Latitude",
                     x_breaks = day_breaks, x_labels = day_labels,
                     y_breaks = latitude_breaks, y_labels = latitude_labels),
    file.path(figure_dir, "temp_altitude_ols.jpeg"),
    ncol = 3, nrow = 2, width = 14, height = 6.66, legend = "right"
  )
  save_panels(
    marginal_panels(fit, titles, side = "t", x = "Day",
                    x_breaks = day_breaks, x_labels = day_labels)[shown],
    file.path(figure_dir, "temp_altitude_ols_day.jpeg"),
    ncol = 4, nrow = 1, width = 14, height = 3.33, legend = "top", common_legend = TRUE
  )
  save_panels(
    marginal_panels(fit, titles, side = "s", x = "Latitude",
                    x_breaks = latitude_breaks, x_labels = latitude_labels)[shown],
    file.path(figure_dir, "temp_altitude_ols_altitude.jpeg"),
    ncol = 4, nrow = 1, width = 14, height = 3.33, legend = "top", common_legend = TRUE
  )
  plot_latitude_summer_3d(fit, figure_dir)
}

plot_latitude_summer_3d <- function(fit, figure_dir = "output/aemet") {
  titles <- list(
    "Intercept", "Year > 2008", expression("CO"[2] * " Emission"),
    "Wind Speed", "Max Precipitation"
  )
  files <- tempfile(paste0("surface", seq_along(fit$betaHat)), fileext = ".jpeg")
  on.exit(unlink(files), add = TRUE)
  for (p in seq_along(fit$betaHat)) {
    save_surface_3d(
      fit$betaHat[[p]],
      files[p],
      x_title = "Latitude",
      y_title = "Day",
      z_title = "Coef"
    )
  }
  save_panels(
    Map(jpeg_panel, files, titles),
    file.path(figure_dir, "temp_altitude_ols_3d.jpeg"),
    ncol = 3, nrow = 2, width = 15, height = 9, legend = "none"
  )
}

bootstrap.opt <- list(B = 50, M = 10000, method = "nonparametric")

data_file <- file.path(demo_dir, "data", "aemet_dat.RData")
if (!file.exists(data_file)) {
  stop("data/aemet_dat.RData is missing.\n",
       "Build it first with:  Rscript data/prepare_aemet.R\n",
       "See data/README.md for the AEMET API key setup.",
       call. = FALSE)
}
load(data_file)
out_dir <- file.path(demo_dir, "output", "aemet")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
summer_days <- 153:273

latitude <- prepare_latitude_summer(dat, summer_days)
fit_prec <- fit_temperature(
  temperature ~ whether_late + co2_ems + wind + max_prec,
  latitude$data,
  S = latitude$S,
  knots = c(10, 8),
  bootstrap.opt = bootstrap.opt
)
save(fit_prec, file = file.path(out_dir, "temp_altitude_ols.RData"))
gc()
plot_latitude_summer(fit_prec, out_dir)
