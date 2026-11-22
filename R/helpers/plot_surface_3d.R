## 3D coefficient surface, the plotly surface used in demo_simu.R.
## x follows the columns of beta (S), y follows the rows (T).

surface_3d_dir <- local({
  ofile <- NULL
  for (i in seq_len(sys.nframe())) {
    candidate <- sys.frame(i)$ofile
    if (!is.null(candidate)) ofile <- candidate
  }
  dirname(normalizePath(ofile))
})

save_surface_3d <- function(beta, path, x_title = "S", y_title = "T", z_title = "Coef") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  matrix_path <- tempfile(fileext = ".csv")
  on.exit(unlink(matrix_path), add = TRUE)
  utils::write.table(beta, matrix_path, sep = ",", row.names = FALSE, col.names = FALSE)
  status <- system2(
    Sys.which("python3"),
    c(file.path(surface_3d_dir, "plot_surface_3d.py"),
      matrix_path, normalizePath(path, mustWork = FALSE), x_title, y_title, z_title)
  )
  if (status != 0) stop("plotly surface export failed for ", path)
  invisible(path)
}
