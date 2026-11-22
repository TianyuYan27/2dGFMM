## Synthetic demos for one Gaussian sample and one binomial sample.
## Each demo fits the non-robust and robust estimators, then tests H0: beta(s,t) = 0.
## From the Github directory:
##   Rscript demo_synthetic.R

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
source(file.path(demo_dir, "R", "src", "main.R"))
source(file.path(demo_dir, "R", "helpers", "sim_generate.R"))
source(file.path(demo_dir, "R", "helpers", "evaluate.R"))
source(file.path(demo_dir, "R", "helpers", "plot_synthetic.R"))

## H0: beta(s,t) = 0 for the intercept and the slope.
test_coefficients <- function(fits, const = 0) {
  rows <- lapply(fits, function(item) {
    p <- pvalues.gfmm2d(item$fit, const = const, tol = 0.001, delta = 0)
    data.frame(
      method = item$name,
      intercept = unname(p[1]),
      slope = unname(p[2]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

run_synthetic_demo <- function(family, n, S, Tlen, bootstrap.opt, out_dir) {
  set.seed(123)
  data <- sim_generate(
    family = family, bivariate = TRUE, n = n, S = S, T = Tlen,
    type = 2, rho = 1.5,
    SNR_B = 1, SNR_sigma = 1
  )
  beta_true <- true_beta(S, Tlen, type = 2, family = family)
  knots <- c(max(round(S / 4), 4), min(round(Tlen / 4), 35))
  method_names <- if (family == "gaussian") c("TSE-OLS", "TSE-Huber") else c("TSE-ML", "TSE-RQL")

  if (family == "binomial") {
    eta_var <- estimate_eta_var_local(
      data = data, formula = Y ~ X, S = S, Tlen = Tlen, silence = TRUE
    )$eta.var
  }

  fit_one <- function(robust) {
    if (family == "gaussian") {
      gfmm2d_gaussian(
        formula = Y ~ X, data = data, S = S,
        robust = robust, smoother = "sandwich", knots = knots,
        bootstrap.opt = bootstrap.opt, pcb = TRUE, scb = TRUE, parallel = FALSE
      )
    } else {
      gfmm2d_binomial(
        formula = Y ~ X, data = data, S = S,
        robust = robust, smoother = "sandwich", knots = knots,
        eta.var = eta_var,
        bootstrap.opt = bootstrap.opt, pcb = TRUE, scb = TRUE, parallel = FALSE
      )
    }
  }
  fits <- list(
    list(name = method_names[1], fit = fit_one(FALSE)),
    list(name = method_names[2], fit = fit_one(TRUE))
  )
  metrics <- do.call(rbind, lapply(fits, function(item) {
    evaluate_fit(item$fit, beta_true, item$name)
  }))
  pvalues <- test_coefficients(fits)

  cat("\n", family, " estimation\n", sep = "")
  print(metrics)
  cat("\n", family, " hypothesis tests, H0: beta(s,t) = 0\n", sep = "")
  print(pvalues)

  save(fits, beta_true, metrics, pvalues, file = file.path(out_dir, paste0(family, "_fit.RData")))
  plot_synthetic(beta_true, fits, out_dir, family)
  invisible(list(metrics = metrics, pvalues = pvalues))
}

n <- 50
S <- 10
Tlen <- 100
bootstrap.opt <- list(B = 50, M = 10000, method = "nonparametric")

out_dir <- file.path(demo_dir, "output", "synthetic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
gaussian <- run_synthetic_demo("gaussian", n, S, Tlen, bootstrap.opt, out_dir)
binomial <- run_synthetic_demo("binomial", n, S, Tlen, bootstrap.opt, out_dir)
cat("synthetic demos finished\n")
