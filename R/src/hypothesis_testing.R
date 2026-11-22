#' P-values for testing H0: beta(s,t) = const
#'
#' @param object gfmm2d model
#' @param const the constant value for the null hypothesis, default is 0 to test no effects
#' @param tol tolerance of difference for binary search algorithm, default is 0.01
#' @param sig_level significance level for the hypothesis test, default is 0.05
#' @param delta difference between alpha from upper bound and lower bound, only used for algorithm correction
#' @return p-values of all covariates in gfmm2d model

pvalues.gfmm2d <- function(object, const = 0, tol = 0.01, sig_level = 0.05, delta = 0.6) {
  vm <- object$vm
  p_values <- rep(0, nrow(vm))
  for (p in 1:length(p_values)) {
    betaHat <- object$betaHat[[p]]
    betaHat.var <- object$betaHat.cov[, , p]

    alpha_left <- 0
    alpha_right <- 1
    while (alpha_right - alpha_left > tol) {
      alpha <- (alpha_left + alpha_right) / 2
      qn <- quantile(vm[p, ], 1 - alpha)
      betaHat.up <- as.vector(betaHat) + qn * sqrt(diag(betaHat.var))
      # Widen the band when any upper bound falls below const.
      if (any(betaHat.up < const)) {
        alpha_right <- alpha
      } else {
        # Narrow the band when every upper bound stays above const.
        alpha_left <- alpha
      }
    }
    alpha_high <- (alpha_left + alpha_right) / 2

    alpha_left <- 0
    alpha_right <- 1
    while (alpha_right - alpha_left > tol) {
      alpha <- (alpha_left + alpha_right) / 2
      qn <- quantile(vm[p, ], 1 - alpha)
      betaHat.low <- as.vector(betaHat) - qn * sqrt(diag(betaHat.var))
      # Widen the band when any lower bound falls above const.
      if (any(betaHat.low > const)) {
        alpha_right <- alpha
      } else {
        # Narrow the band when every lower bound stays below const.
        alpha_left <- alpha
      }
    }
    alpha_low <- (alpha_left + alpha_right) / 2

    print(paste0("alpha_low =", alpha_low, " alpha_high =", alpha_high))

    # When the bounds are close and at least one alpha is below the significance
    # level, report the larger alpha if only one side is significant, and the
    # smaller alpha if both are. Otherwise report the smaller alpha.
    if ((alpha_low < sig_level | alpha_high < sig_level) & abs(alpha_high - alpha_low) < delta) {
      if (alpha_low < sig_level & alpha_high < sig_level) {
        p_values[p] <- min(alpha_high, alpha_low)
      } else {
        p_values[p] <- max(alpha_high, alpha_low)
      }
    } else {
      p_values[p] <- min(alpha_high, alpha_low)
    }
  }
  p_values
}
