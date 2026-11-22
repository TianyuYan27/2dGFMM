helpers_dir <- local({
  ofile <- NULL
  for (i in seq_len(sys.nframe())) {
    candidate <- sys.frame(i)$ofile
    if (!is.null(candidate)) ofile <- candidate
  }
  dirname(normalizePath(ofile))
})
source(file.path(helpers_dir, "plot_bands.R"))
source(file.path(helpers_dir, "plot_surface_3d.R"))

save_surface_panels <- function(betas, titles, path, width, height) {
  files <- tempfile(paste0("surface", seq_along(betas)), fileext = ".jpeg")
  on.exit(unlink(files), add = TRUE)
  for (p in seq_along(betas)) {
    save_surface_3d(
      betas[[p]], files[p],
      x_title = "S", y_title = "T", z_title = "Coef"
    )
  }
  save_panels(
    Map(jpeg_panel, files, titles),
    path,
    ncol = length(betas), nrow = 1, width = width, height = height, legend = "none"
  )
}

synthetic_axes <- function(fit) {
  n_t <- nrow(fit$betaHat[[1]])
  n_s <- ncol(fit$betaHat[[1]])
  x_breaks <- unique(c(1, round(n_t / 3), round(2 * n_t / 3), n_t))
  y_breaks <- unique(c(1, round(n_s / 2), n_s))
  list(
    x_breaks = x_breaks,
    x_labels = sprintf("%.2f", (x_breaks - 1) / (n_t - 1)),
    y_breaks = y_breaks,
    y_labels = sprintf("%.2f", (y_breaks - 1) / (n_s - 1))
  )
}

save_synthetic_bands <- function(fits, out_dir, prefix) {
  ax <- synthetic_axes(fits[[1]]$fit)

  bivariate <- do.call(c, lapply(fits, function(item) {
    bivariate_panels(
      item$fit, c("Intercept", "Slope"), x = "T", y = "S",
      x_breaks = ax$x_breaks, x_labels = ax$x_labels,
      y_breaks = ax$y_breaks, y_labels = ax$y_labels
    )
  }))
  save_panels(
    bivariate, file.path(out_dir, paste0(prefix, "_bands_bivariate.jpeg")),
    ncol = 2, nrow = 2, width = 10, height = 7.2, legend = "right"
  )

  marginal <- do.call(c, lapply(fits, function(item) {
    marginal_panels(
      item$fit, c("Intercept", "Slope"), side = "t", x = "T",
      x_breaks = ax$x_breaks, x_labels = ax$x_labels
    )
  }))
  save_panels(
    marginal, file.path(out_dir, paste0(prefix, "_bands_marginal_t.jpeg")),
    ncol = 2, nrow = 2, width = 10, height = 7.2, legend = "top",
    common_legend = TRUE
  )
}

plot_synthetic <- function(beta_true, fits, out_dir, prefix) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  save_surface_panels(
    beta_true,
    list("True intercept", "True slope"),
    file.path(out_dir, paste0(prefix, "_true.jpeg")),
    width = 10, height = 4.6
  )
  save_surface_panels(
    lapply(fits, function(item) item$fit$betaHat[[2]]),
    lapply(fits, function(item) item$name),
    file.path(out_dir, paste0(prefix, "_slope.jpeg")),
    width = 10, height = 4.6
  )
  save_synthetic_bands(fits, out_dir, prefix)
}
