## ISE is the grid average of squared error.
## PCB is the pointwise band estimate +/- 2 standard errors.
## SCB is the simultaneous band estimate +/- qn standard errors.
## Coverage is the proportion of grid points whose band contains the truth.
## IAW is the average band width.

integrated_squared_error <- function(beta_hat, beta_true) {
  ise <- vapply(seq_along(beta_hat), function(p) {
    mean((beta_hat[[p]] - beta_true[[p]])^2)
  }, numeric(1))
  names(ise) <- c("intercept", "slope")[seq_along(ise)]
  ise
}

calculate_cover <- function(fitted_model, beta_true) {
  if (is.null(fitted_model$betaHat.cov)) {
    stop("calculate_cover requires betaHat.cov. Fit with pcb = TRUE and scb = TRUE.")
  }
  P <- length(fitted_model$betaHat)
  cover_pcb <- iaw_pcb <- cover_scb <- iaw_scb <- rep(NA_real_, P)
  names(cover_pcb) <- names(iaw_pcb) <- names(cover_scb) <- names(iaw_scb) <-
    c("intercept", "slope")[seq_len(P)]

  for (p in seq_len(P)) {
    estimate <- as.vector(fitted_model$betaHat[[p]])
    truth <- as.vector(beta_true[[p]])
    se <- sqrt(diag(fitted_model$betaHat.cov[, , p]))
    qn <- fitted_model$qn[p]

    pcb_up <- estimate + 2 * se
    pcb_low <- estimate - 2 * se
    cover_pcb[p] <- mean(truth <= pcb_up & truth >= pcb_low)
    iaw_pcb[p] <- mean(pcb_up - pcb_low)

    scb_up <- estimate + qn * se
    scb_low <- estimate - qn * se
    cover_scb[p] <- mean(truth <= scb_up & truth >= scb_low)
    iaw_scb[p] <- mean(scb_up - scb_low)
  }

  list(cover_pcb = cover_pcb, iaw_pcb = iaw_pcb,
       cover_scb = cover_scb, iaw_scb = iaw_scb)
}

evaluate_fit <- function(fit, beta_true, method) {
  ise <- integrated_squared_error(fit$betaHat, beta_true)
  cover <- calculate_cover(fit, beta_true)
  data.frame(
    method = method,
    ise_intercept = unname(ise["intercept"]),
    ise_slope = unname(ise["slope"]),
    cover_pcb_intercept = unname(cover$cover_pcb["intercept"]),
    cover_pcb_slope = unname(cover$cover_pcb["slope"]),
    iaw_pcb_intercept = unname(cover$iaw_pcb["intercept"]),
    iaw_pcb_slope = unname(cover$iaw_pcb["slope"]),
    cover_scb_intercept = unname(cover$cover_scb["intercept"]),
    cover_scb_slope = unname(cover$cover_scb["slope"]),
    iaw_scb_intercept = unname(cover$iaw_scb["intercept"]),
    iaw_scb_slope = unname(cover$iaw_scb["slope"]),
    stringsAsFactors = FALSE
  )
}
