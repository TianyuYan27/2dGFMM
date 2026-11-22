## Synthetic longitudinal functional data: sim_generate(), intercept(), f1(), f2().
#' Simulate longitudinal functional response Y (matrix) and a scalar predictor X (vector)
#' @param family distribution of longitudinal functional data, including "gaussian", "binomial".
#' @param bivariate whether the beta is 2d.
#' @param n number of subjects.
#' @param S number of S grid points on the functional domain.
#' @param T number of T grid points on the functional domain.
#' @param type 1 uses the localized slope f1(); 2 uses the smooth slope f2().
#' @param SNR_B relative importance of random effects
#' @param SNR_sigma signal-to-noise ratio in Gaussian data generation
#' @return a data frame containing generated predictors and longitudinal functional outcomes.

sim_generate <- function(family = "gaussian", bivariate = TRUE, n = 100, S = 10, T = 100,
                         type, rho = 1.5,
                         SNR_B = 1, SNR_sigma = 1) {
  library(mvtnorm)

  subj <- as.factor(rep(1:n, each = S))
  S.grid <- seq(0, 1, length.out = S)
  T.grid <- seq(0, 1, length.out = T)
  beta0 <- intercept(S = S.grid, T = T.grid)
  if (type == 1) {
    beta1 <- f1(S = S.grid, T = T.grid)
  } else {
    beta1 <- f2(S = S.grid, T = T.grid)
  }
  if (family == "binomial") {
    ## Binomial coefficient surfaces are scaled so the linear predictor stays
    ## away from complete separation. Gaussian surfaces are unchanged.
    binomial_beta_scale <- 0.25
    beta0 <- binomial_beta_scale * beta0
    beta1 <- binomial_beta_scale * beta1
  }

  phi <- array(0, dim = c(2, T))
  phi[1, ] <- (1.5 - sin(2 * T.grid * pi) - cos(2 * T.grid * pi))
  phi[1, ] <- phi[1, ] / sqrt(sum(phi[1, ]^2))
  phi[2, ] <- sin(4 * T.grid * pi)
  phi[2, ] <- phi[2, ] / sqrt(sum(phi[2, ]^2))

  phi1 <- array(0, dim = c(2, T))
  phi1[1, ] <- cos(2 * pi * T.grid)
  phi1[1, ] <- phi1[1, ] / sqrt(sum(phi1[1, ]^2))
  phi1[2, ] <- sin(2 * pi * T.grid)
  phi1[2, ] <- phi1[2, ] / sqrt(sum(phi1[2, ]^2))

  X_des <- lapply(1:n, function(i) rnorm(S, 0, 2))
  if (bivariate) {
    Xbeta <- lapply(1:n, function(i) X_des[[i]] * beta1)
    fixef <- lapply(1:n, function(i) beta0 + Xbeta[[i]])
  } else {
    beta0 <- apply(beta0, 2, mean)
    beta1 <- apply(beta1, 2, mean)
    fixef <- lapply(1:n, function(i) matrix(rep(beta0, S), nrow = S, byrow = TRUE) + X_des[[i]] %*% t(beta1))
  }

  a <- rmvnorm(n, mean = rep(0, 2), sigma = diag(c(3, 1.5)))
  u0 <- a %*% phi
  a1 <- rmvnorm(n, mean = rep(0, 2), sigma = diag(c(3, 1.5)))
  u1 <- a1 %*% phi1

  Z_des <- lapply(1:n, function(i) 6 * (S.grid - 0.5) + rnorm(S, 0, rho))
  ranef <- lapply(1:n, function(i) t(replicate(S, u0[i, ])) + t(replicate(S, u1[i, ])) * Z_des[[i]])

  ranef <- lapply(1:n, function(i) sd(fixef[[i]]) / sd(ranef[[i]]) / SNR_B * ranef[[i]])
  Y_true <- lapply(1:n, function(i) fixef[[i]] + ranef[[i]])
  p_true <- lapply(1:n, function(i) plogis(Y_true[[i]]))

  Y_obs <- lapply(1:n, function(i) array(0, dim = c(S, T)))
  for (i in 1:n) {
    sd_signal <- sd(Y_true[[i]])
    for (s in 1:S) {
      for (t in 1:T) {
        if (family == "gaussian") {
          Y_obs[[i]][s, t] <- rnorm(1, mean = Y_true[[i]][s, t], sd = sd_signal / SNR_sigma)
        } else if (family == "binomial") {
          Y_obs[[i]][s, t] <- rbinom(1, 1, p_true[[i]][s, t])
        }
      }
    }
  }
  X <- unlist(X_des)
  Y <- do.call("rbind", Y_obs)
  data.frame(ID = subj, X = X, Longit = rep(S.grid, n), Y = I(Y))
}

## Coefficient surfaces on the same orientation as gfmm2d_gaussian()$betaHat: T rows, S columns.
true_beta <- function(S, T, type, family = "gaussian") {
  S.grid <- seq(0, 1, length.out = S)
  T.grid <- seq(0, 1, length.out = T)
  intercept_mat <- intercept(S = S.grid, T = T.grid)
  slope <- if (type == 1) f1(S = S.grid, T = T.grid) else f2(S = S.grid, T = T.grid)
  if (family == "binomial") {
    binomial_beta_scale <- 0.25
    intercept_mat <- binomial_beta_scale * intercept_mat
    slope <- binomial_beta_scale * slope
  }
  list(t(intercept_mat), t(slope))
}

intercept <- function(S, T) {
  intercept <- array(0, dim = c(length(S), length(T)))
  for (s in 1:length(S)) {
    for (t in 1:length(T)) {
      intercept[s, t] <- 3 * sin(pi * (S[s] + 0.5)^2) * cos(1 * pi * (T[t]) + 0.5) + 1
    }
  }
  intercept
}

f1 <- function(S, T) {
  y <- array(0, dim = c(length(S), length(T)))
  for (s in 1:length(S)) {
    for (t in 1:length(T)) {
      y[s, t] <- sin(0.8 * pi * (S[s] + 0.5)^2) * cos(4 * pi * (T[t]))
    }
  }
  block <- y[(2 * length(S) / 10):(5 * length(S) / 10), (14 * length(T) / 100):(38 * length(T) / 100)]
  f <- array(0, dim = c(length(S), length(T)))
  f[(1 * length(S) / 10):(4 * length(S) / 10), (14 * length(T) / 100):(38 * length(T) / 100)] <- 5 * block
  f[(7 * length(S) / 10):(10 * length(S) / 10), (14 * length(T) / 100):(38 * length(T) / 100)] <- -5 * block
  f[(1 * length(S) / 10):(4 * length(S) / 10), (62 * length(T) / 100):(86 * length(T) / 100)] <- -5 * block
  f[(7 * length(S) / 10):(10 * length(S) / 10), (62 * length(T) / 100):(86 * length(T) / 100)] <- 5 * block
  f
}

f2 <- function(S, T) {
  y <- array(0, dim = c(length(S), length(T)))
  for (s in 1:length(S)) {
    for (t in 1:length(T)) {
      y[s, t] <- 5 * sin(0.5 * pi * (S[s] + 0.5)^2) * cos(2 * pi * (T[t]) + 0.5)
    }
  }
  y
}
