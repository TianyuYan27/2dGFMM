## ============================================================================
## Binomial two-stage estimator with a logistic-normal integral link.
##
## Main entry point: gfmm2d_binomial().
##
## At every grid point (s, t) the marginal mean is the logit probability
## integrated over a latent random effect eta ~ N(0, v(s, t)). The link is
## m_v(q) = E[expit(q + sqrt(v) Z)], whose inverse is used as the GLM link so
## that marginal responses can estimate the *conditional* coefficient surface.
##
## The latent variance v(s, t) is estimated once, before the two-stage fit, by
## a local random-intercept model on a configurable window. The default is
## 1 x 3 along the time axis (the target visit and one time step on either
## side), set by variance.control$window = c(1, 3):
##   logit P(Y_i(s, t_j) = 1 | b_i) = x_i(s, t_j)' theta + b_i,
##   b_i ~ N(0, v),  j in {b - 1, b, b + 1}.
## The variance component from glmer() is v(s, t). If the caller supplies
## eta.var, that estimate is used instead.
##
## A common variance at (s, t) does not represent a random slope whose latent
## variance depends on the subject. Do not pass a shared variance for that
## model.
##
## Smoothing is the original sandwich or tensor-product step. Pointwise bands
## resample whole subjects and refit both stages. Simultaneous bands use the
## original marginal FPCA and B-spline Monte Carlo. eta.var stays fixed across
## bootstrap replicates. The binomial factor gamma = 1.5 is the original
## bootstrap inflation; it is not a correction from the marginal coefficient to
## the conditional coefficient.
## ============================================================================


## ============================================================================
## Main function
## ============================================================================

#' Fit a generalized two-dimensional functional mixed model for binary data
#' with a logistic-normal integral link
#'
#' @param formula two-sided formula object in lm() format; the response is a
#'   matrix with one row per (subject, visit) and one column per time point
#' @param data data frame containing the variables in formula and an ID column;
#'   covariates may be vectors or matrices whose columns index the time axis
#' @param S number of visits per subject
#' @param robust whether to use M-estimation (Cantoni-Ronchetti Mqle) at the
#'   pointwise stage
#' @param smoother sandwich smoother ("sandwich") or tensor product smoother
#'   ("te")
#' @param knots c(Sknots, Tknots), the number of knots for the S and T domains
#' @param eta.var optional latent variance on the logit scale. Either a scalar
#'   or a strict S x T matrix. If NULL (default), v(s, t) is estimated from the
#'   data by a local random-intercept model; see variance.control
#' @param variance.control list controlling local variance estimation:
#'   window = c(S_window, T_window), the number of grid points in the local
#'   window along the visit and time axes (default c(1, 3), i.e. the target
#'   visit and one time step on either side; both entries must be odd),
#'   hs = optional half-width on the visit axis that overrides window[1],
#'   ht = optional half-width on the time axis that overrides window[2],
#'   nAGQ = adaptive Gauss-Hermite nodes used by glmer (default 7)
#' @param bootstrap.opt bootstrapping options: B = number of subject
#'   replicates, M = Monte Carlo draws for simultaneous bands, method =
#'   "nonparametric"
#' @param bandwidth.control correlation correction and adjustment for
#'   confidence band bandwidth control: d = tail trimming, adjust = variance
#'   adjustment
#' @param gh.n number of Gauss-Hermite quadrature nodes for the integral link
#' @param robust.c tuning constant for the Huber score in Mqle
#' @param robust.wx covariate weights for glmrob: "none" or "covMcd"
#' @param pcb whether to obtain covariance estimates and pointwise confidence
#'   bands, default is TRUE
#' @param scb whether to obtain simultaneous confidence bands, default is FALSE
#' @param parallel whether to run bootstrap replications and local variance
#'   fits in parallel (Linux/macOS)
#' @param silence whether to suppress progress messages, default is FALSE
#' @param keep_diagnostics whether to return pointwise fit diagnostics
#' @export
#' @return a list containing betaHat (smoothed coefficient surfaces),
#'   betaTilde (raw pointwise surfaces), betaHat.cov (bootstrap covariance
#'   surfaces), qn (band critical values), vm (bootstrap draws used for
#'   p-values), eta.var (the variance surface actually used),
#'   variance_diagnostics, variance_summary, fit_diagnostics,
#'   estimation_settings, and timing information

gfmm2d_binomial <- function(formula,
                            data,
                            S,
                            robust = FALSE,
                            smoother = "sandwich",
                            knots = NULL,
                            eta.var = NULL,
                            variance.control = list(window = c(1, 3), hs = NULL, ht = NULL, nAGQ = 7L),
                            bootstrap.opt = list(B = 50, M = 10000, method = "nonparametric"),
                            bandwidth.control = list(d = 2, adjust = NULL),
                            gh.n = 40L,
                            robust.c = 1.345,
                            robust.wx = "none",
                            pcb = TRUE,
                            scb = FALSE,
                            parallel = FALSE,
                            silence = FALSE,
                            keep_diagnostics = TRUE)
{
  if (isTRUE(scb) && !isTRUE(pcb)) {
    stop("Simultaneous bands reuse the pointwise bootstrap fit; set pcb = TRUE.")
  }
  if (!robust.wx %in% c("none", "covMcd")) {
    stop('robust.wx must be either "none" or "covMcd".')
  }
  if (length(robust.c) != 1 || !is.finite(robust.c) || robust.c <= 0) {
    stop("robust.c must be a positive finite number.")
  }
  if (!smoother %in% c("sandwich", "te")) {
    stop('smoother must be either "sandwich" or "te".')
  }

  model_formula <- as.character(formula)
  stopifnot(model_formula[1] == "~" & length(model_formula) == 3)
  vars <- unlist(strsplit(model_formula[3], split = ".\\+."))
  response <- data[, model_formula[2]]
  n <- nrow(response) / S
  if (n != as.integer(n)) stop("The number of response rows is not divisible by S.")
  n <- as.integer(n)
  Tlen <- ncol(response)
  if (is.null(Tlen)) stop("The response must be a matrix with one column per time point.")

  eta_var_supplied <- !is.null(eta.var)
  gfmm2d_binomial_dependencies(need_lme4 = !eta_var_supplied)

  ## --------------------------------------------------------------------------
  ## Step 0: latent variance v(s, t)
  ## --------------------------------------------------------------------------
  t_start <- proc.time()
  variance_diagnostics <- NULL
  variance_summary <- NULL
  if (!eta_var_supplied) {
    if (!silence) cat("Step 0: local random-intercept variance\n")
    vc <- variance.control
    if (is.null(vc$window)) vc$window <- c(1L, 3L)
    if (is.null(vc$nAGQ)) vc$nAGQ <- 7L
    vfit <- estimate_eta_var_local(
      data = data, formula = formula, S = S, Tlen = Tlen,
      window = vc$window, hs = vc$hs, ht = vc$ht, nAGQ = vc$nAGQ,
      parallel = parallel, silence = silence
    )
    eta.var <- vfit$eta.var
    variance_diagnostics <- vfit$grid
    variance_summary <- vfit$summary
    if (!silence) {
      cat(sprintf(
        "  local var mean %.3f  range [%.3f, %.3f]  singular %d  failed %d\n",
        variance_summary$mean_v, variance_summary$min_v, variance_summary$max_v,
        variance_summary$n_singular, variance_summary$n_fail
      ))
    }
  }
  eta.var <- validate_eta_var(eta.var, S, Tlen)
  t_var <- proc.time()
  quad <- logistic_normal_quad(gh.n)

  t_boot <- NULL
  Y <- lapply(seq_len(n), function(i) response[((i - 1) * S + 1):(i * S), ])
  var.mat <- NULL
  for (v in vars) {
    if (is.array(data[, v]) || is.matrix(data[, v])) var.mat <- c(var.mat, v)
  }

  ## --------------------------------------------------------------------------
  ## Step 1: logistic-normal pointwise estimates
  ## --------------------------------------------------------------------------
  if (!silence) cat("Step 1: logistic-normal pointwise estimate\n")
  link_cache <- new.env(parent = emptyenv())
  link_for <- function(v) {
    key <- sprintf("%.12g", v)
    hit <- link_cache[[key]]
    if (!is.null(hit)) return(hit)
    link_cache[[key]] <- logistic_normal_link(v, quad)
    link_cache[[key]]
  }

  ## Fit one pointwise model at (s, t).
  bim <- function(s, t) {
    Yst <- sapply(Y, function(x) x[s, t])
    Xst <- data[seq(s, (n - 1) * S + s, by = S), vars, drop = FALSE]
    if (length(vars) == 1) {
      if (!is.null(var.mat)) {
        Xst <- array(Xst[, t], dim = c(n, 1))
      } else {
        Xst <- array(unlist(Xst), dim = c(n, 1))
      }
      colnames(Xst) <- vars
    } else if (!is.null(var.mat)) {
      Xst.sub <- NULL
      for (p.mat in seq_along(var.mat)) {
        Xst.sub <- cbind(Xst.sub, Xst[, var.mat[p.mat]][, t])
      }
      colnames(Xst.sub) <- var.mat
      var.others <- setdiff(names(Xst), var.mat)
      Xst <- cbind(Xst[, var.others, drop = FALSE], Xst.sub)
    }
    dat <- data.frame(Xst, Yst, check.names = FALSE)
    term_names <- colnames(stats::model.matrix(
      as.formula(paste(model_formula[1], model_formula[3])), data = dat
    ))[-1]
    form_bi <- as.formula(paste0("Yst ~ ", model_formula[3]))
    v_st <- eta.var[s, t]
    empty <- stats::setNames(rep(NA_real_, length(term_names) + 1), c("intercept", term_names))
    diagnostic <- data.frame(
      s = s, t = t, n_used = sum(is.finite(Yst)), method = if (robust) "glmrob" else "glm",
      link = if (v_st == 0) "logit" else "logistic-normal", eta_var = v_st,
      converged = FALSE, design_rank = NA_integer_, near_boundary = FALSE,
      failure_reason = NA_character_, warning_messages = "", iterations = NA_integer_,
      stringsAsFactors = FALSE
    )

    fail <- function(reason, rank = NA_integer_, warnings = "") {
      diagnostic$failure_reason <<- reason
      diagnostic$design_rank <<- rank
      diagnostic$warning_messages <<- paste(warnings, collapse = " | ")
      list(beta = empty, diagnostic = diagnostic)
    }

    if (all(stats::na.omit(dat$Yst) == 0) || all(stats::na.omit(dat$Yst) == 1)) {
      return(fail("constant_response"))
    }
    Xmm <- try(stats::model.matrix(form_bi, dat), silent = TRUE)
    if (inherits(Xmm, "try-error")) return(fail("design_matrix"))
    rank <- qr(Xmm)$rank
    diagnostic$design_rank <- rank
    if (rank < ncol(Xmm)) return(fail("rank_deficient", rank))

    warns <- character()
    near_boundary <- FALSE
    note_warning <- function(w) {
      msg <- conditionMessage(w)
      warns <<- c(warns, msg)
      if (grepl("fitted probabilities numerically 0 or 1", msg, fixed = TRUE)) {
        near_boundary <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
    family_st <- stats::binomial(link = link_for(v_st))
    if (robust) {
      fit_bi <- try(withCallingHandlers(
        robustbase::glmrob(
          formula = form_bi, data = dat, family = family_st, method = "Mqle",
          weights.on.x = robust.wx,
          control = robustbase::glmrobMqle.control(tcc = robust.c, maxit = 100, acc = 1e-6)
        ),
        warning = note_warning
      ), silent = TRUE)
    } else {
      fit_bi <- try(withCallingHandlers(
        stats::glm(
          formula = form_bi, data = dat, family = family_st,
          control = stats::glm.control(maxit = 100, epsilon = 1e-8)
        ),
        warning = note_warning
      ), silent = TRUE)
    }
    diagnostic$near_boundary <- near_boundary
    diagnostic$warning_messages <- paste(warns, collapse = " | ")
    if (inherits(fit_bi, "try-error")) {
      return(fail(paste("fit_error:", as.character(fit_bi)), rank, warns))
    }
    beta_raw <- tryCatch(stats::coef(fit_bi), error = function(e) NULL)
    diagnostic$iterations <- if (!is.null(fit_bi$iter)) fit_bi$iter else NA_integer_
    converged <- isTRUE(fit_bi$converged)
    diagnostic$converged <- converged
    ## A point with fitted probabilities numerically 0 or 1 is (quasi-)
    ## separated: neither the GLM nor the Mqle equation identifies the slope.
    ## glmrob can report convergence at an arbitrarily large root (|beta| in
    ## the thousands), and a single such point is spread over a whole region
    ## by the least-squares smoother. Drop the point and let the existing
    ## interpolation from its neighbours fill it in.
    if (near_boundary) return(fail("near_boundary", rank, warns))
    if (!converged) return(fail("not_converged", rank, warns))
    if (is.null(beta_raw) || length(beta_raw) != length(empty) || any(!is.finite(beta_raw))) {
      return(fail("non_finite_coefficient", rank, warns))
    }
    names(beta_raw) <- names(empty)
    diagnostic$failure_reason <- NA_character_
    list(beta = beta_raw, diagnostic = diagnostic)
  }

  massm <- vector("list", S)
  for (s in seq_len(S)) {
    massm[[s]] <- lapply(seq_len(Tlen), function(t) bim(s, t))
  }
  beta_mat <- lapply(massm, function(one_s) {
    do.call(rbind, lapply(one_s, function(pt) as.data.frame(as.list(pt$beta), check.names = FALSE)))
  })
  betaTilde <- do.call(rbind, beta_mat)
  diagnostics <- if (keep_diagnostics) {
    do.call(rbind, lapply(massm, function(one_s) do.call(rbind, lapply(one_s, `[[`, "diagnostic"))))
  } else {
    NULL
  }

  ## --------------------------------------------------------------------------
  ## Step 2: bivariate smoothing
  ## --------------------------------------------------------------------------
  if (!silence) cat("Step 2: Bivariate smoothing\n")
  if (is.null(knots)) knots <- c(min(round(S / 1.5), 35), min(round(Tlen / 4), 35))
  Sknots <- knots[1]
  Tknots <- knots[2]
  S.argvals <- seq_len(S)
  T.argvals <- seq_len(Tlen)
  sL <- refund:::pspline.setting(
    x = S.argvals, knots = refund:::select_knots(S.argvals, Sknots),
    p = 3, m = 2, periodicity = FALSE, weight = NULL
  )
  tL <- refund:::pspline.setting(
    x = T.argvals, knots = refund:::select_knots(T.argvals, Tknots),
    p = 3, m = 2, periodicity = FALSE, weight = NULL
  )
  B1 <- tL$B
  Bt1 <- Matrix::Matrix(t(as.matrix(B1)))
  P1 <- tL$P
  B2 <- sL$B
  Bt2 <- Matrix::Matrix(t(as.matrix(B2)))
  P2 <- sL$P

  betaHat <- array(0, dim = c(Tlen * S, ncol(betaTilde)))
  lambda <- array(0, dim = c(2, ncol(betaTilde)))
  S1 <- S2 <- S2kS1 <- vector("list", ncol(betaTilde))
  na_before <- numeric(ncol(betaTilde))
  names(na_before) <- names(betaTilde)
  for (p in seq_len(ncol(betaTilde))) {
    sm <- array(betaTilde[[p]], dim = c(Tlen, S))
    na_before[p] <- mean(!is.finite(sm))
    if (all(!is.finite(sm))) {
      stop("Coefficient ", names(betaTilde)[p], " has no valid estimate at any grid point.")
    }
    if (any(!is.finite(sm))) sm <- zoo::na.approx(sm, na.rm = FALSE, rule = 2)
    if (any(!is.finite(sm))) {
      stop("Coefficient ", names(betaTilde)[p], " still has non-finite values after interpolation.")
    }
    y <- refund::fbps(sm, knots = c(Tknots, Sknots))
    if (smoother == "te") {
      te.sm <- data.frame(y = as.vector(sm), t = rep(seq_len(Tlen), S), s = rep(seq_len(S), each = Tlen))
      betaHat[, p] <- mgcv::gam(y ~ te(t, s, k = c(Tknots, Sknots)), data = te.sm)$fitted.values
    } else {
      betaHat[, p] <- as.vector(y$Yhat)
    }
    lambda[, p] <- y$lambda
    S1[[p]] <- as.matrix(B1 %*% solve(Bt1 %*% B1 + lambda[1, p] * P1) %*% Bt1)
    S2[[p]] <- as.matrix(B2 %*% solve(Bt2 %*% B2 + lambda[2, p] * P2) %*% Bt2)
    S2kS1[[p]] <- kronecker(S2[[p]], S1[[p]])
  }
  beta_tilde <- lapply(seq_len(ncol(betaTilde)), function(p) array(betaTilde[[p]], dim = c(Tlen, S)))
  names(beta_tilde) <- names(betaTilde)
  betaHat <- lapply(seq_len(ncol(betaHat)), function(p) array(betaHat[, p], dim = c(Tlen, S)))
  names(betaHat) <- names(betaTilde)
  t_fit <- proc.time()

  settings <- list(
    family = "binomial", robust = robust, smoother = smoother, knots = c(Sknots, Tknots),
    gh.n = as.integer(gh.n), robust.c = robust.c, robust.wx = robust.wx,
    eta_var_mean = mean(eta.var), eta_var_min = min(eta.var), eta_var_max = max(eta.var),
    eta_var_source = if (eta_var_supplied) "user" else "local_glmer",
    variance_window = if (is.null(variance_summary)) NULL else
      c(variance_summary$window_s, variance_summary$window_t),
    gamma = 1.5,
    gamma_meaning = "original binomial bootstrap inflation, not a conditional-coefficient correction"
  )

  if (!pcb) {
    out <- list(
      betaHat = betaHat, betaTilde = beta_tilde,
      eta.var = eta.var,
      variance_diagnostics = variance_diagnostics,
      variance_summary = variance_summary,
      fit_diagnostics = diagnostics, estimation_settings = settings,
      smooth_na_fraction = na_before,
      time = gfmm2d_binomial_elapsed(t_start, t_var, t_fit, t_fit)
    )
    class(out) <- c("gfmm2d_binomial", "gfmm2d")
    return(out)
  }

  ## --------------------------------------------------------------------------
  ## Step 3: subject bootstrap
  ## --------------------------------------------------------------------------
  B <- bootstrap.opt$B
  para <- bootstrap.opt$method
  if (!identical(para, "nonparametric")) stop("Only nonparametric subject bootstrap is implemented.")
  if (!silence) cat("bootstrapping...\n")
  boot_one <- function(b) {
    fmm2d_binomial_boot(
      b = b, data = data, n = n, formula = formula, S = S, robust = robust,
      smoother = smoother, knots = knots, eta.var = eta.var, gh.n = gh.n,
      robust.c = robust.c, robust.wx = robust.wx, parallel = FALSE
    )
  }
  if (parallel) {
    betaHat.boot <- parallel::mclapply(seq_len(B), boot_one, mc.cores = max(1, parallel::detectCores() - 1))
  } else {
    betaHat.boot <- lapply(seq_len(B), boot_one)
  }
  ok <- !vapply(betaHat.boot, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (sum(ok) < 2) stop("Fewer than 2 valid bootstrap replicates.")
  betaHat.boot <- betaHat.boot[ok]
  if (!silence) cat("bootstrapping finished:", sum(ok), "of", B, "\n")
  t_boot <- proc.time()

  gamma <- 1.5
  P <- length(betaHat)
  cov.beta.ts.hat <- array(0, dim = c(Tlen * S, Tlen * S, P))
  for (p in seq_len(P)) {
    draws <- do.call(rbind, lapply(betaHat.boot, function(one) c(one[[p]])))
    cov_raw <- gamma * stats::cov(draws)
    edcomp <- eigen(cov_raw, symmetric = TRUE)
    keep <- which(edcomp$values > 0)
    if (length(keep) == Tlen * S) {
      ker <- cov_raw
    } else {
      vals <- edcomp$values[keep]
      ker <- edcomp$vectors[, keep, drop = FALSE] %*%
        diag(vals, nrow = length(vals)) %*%
        t(edcomp$vectors[, keep, drop = FALSE])
    }
    cov.beta.ts.hat[, , p] <- as.matrix(S2kS1[[p]] %*% ker %*% t(S2kS1[[p]]))
  }

  qn <- rep(NA_real_, P)
  vm <- NULL
  if (isTRUE(scb)) {
    if (is.null(bootstrap.opt$M) || bootstrap.opt$M < 1) {
      stop("bootstrap.opt$M is the number of Monte Carlo draws for simultaneous bands.")
    }
    if (!silence) cat("Step 4: simultaneous confidence bands\n")
    bands <- gfmm2d_binomial_scb(
      betaHat = betaHat, betaHat.boot = betaHat.boot, cov.beta.ts.hat = cov.beta.ts.hat,
      Bt = as.matrix(B1), S = S, Tlen = Tlen, gamma = gamma, M = bootstrap.opt$M,
      bandwidth.control = bandwidth.control
    )
    qn <- bands$qn
    vm <- bands$vm
  }

  out <- list(
    betaHat = betaHat, betaTilde = beta_tilde, betaHat.cov = cov.beta.ts.hat,
    qn = qn, vm = vm,
    eta.var = eta.var,
    variance_diagnostics = variance_diagnostics,
    variance_summary = variance_summary,
    fit_diagnostics = diagnostics, estimation_settings = settings,
    smooth_na_fraction = na_before, bootstrap_ok = sum(ok), bootstrap_requested = B,
    time = gfmm2d_binomial_elapsed(t_start, t_var, t_fit, t_boot)
  )
  class(out) <- c("gfmm2d_binomial", "gfmm2d")
  out
}


## ============================================================================
## Dependencies
## ============================================================================

gfmm2d_binomial_dependencies <- function(need_lme4 = FALSE) {
  needed <- c("refund", "zoo", "Matrix", "parallel", "robustbase", "statmod", "mgcv")
  if (need_lme4) needed <- c(needed, "lme4")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         ". gfmm2d_binomial() does not install packages automatically.")
  }
}


## ============================================================================
## Integral link: quadrature, link object, robust score
## ============================================================================

logistic_normal_quad <- function(gh.n) {
  gh.n <- as.integer(gh.n)
  if (length(gh.n) != 1 || !is.finite(gh.n) || gh.n < 2) {
    stop("gh.n must be a finite integer no smaller than 2.")
  }
  statmod::gauss.quad.prob(n = gh.n, dist = "normal")
}

## Build a class "link-glm" object for m_v. v == 0 returns the ordinary logit
## link. Quadrature nodes already carry standard normal weights: do not
## multiply them by sqrt(2).
logistic_normal_link <- function(v, quad) {
  if (length(v) != 1 || !is.finite(v) || v < 0) {
    stop("logistic_normal_link() requires a finite nonnegative variance v.")
  }
  if (v == 0) return(stats::make.link("logit"))

  z <- quad$nodes
  w <- quad$weights
  shift <- sqrt(v) * z

  linkinv <- function(eta) {
    eta <- as.numeric(eta)
    p <- stats::plogis(outer(eta, shift, "+"))
    as.numeric(p %*% w)
  }
  mu.eta <- function(eta) {
    eta <- as.numeric(eta)
    p <- stats::plogis(outer(eta, shift, "+"))
    as.numeric((p * (1 - p)) %*% w)
  }
  linkfun <- function(mu) {
    mu <- as.numeric(mu)
    vapply(mu, function(m) {
      if (!is.finite(m)) return(NA_real_)
      lo <- -20
      hi <- 20
      flo <- linkinv(lo) - m
      fhi <- linkinv(hi) - m
      for (k in seq_len(8)) {
        if (is.finite(flo) && is.finite(fhi) && flo <= 0 && fhi >= 0) break
        lo <- lo * 2
        hi <- hi * 2
        flo <- linkinv(lo) - m
        fhi <- linkinv(hi) - m
      }
      if (!(flo <= 0 && fhi >= 0)) return(NA_real_)
      stats::uniroot(function(a) linkinv(a) - m, c(lo, hi), tol = 1e-10)$root
    }, numeric(1))
  }
  structure(
    list(
      linkfun = linkfun,
      linkinv = linkinv,
      mu.eta = mu.eta,
      valideta = function(eta) TRUE,
      name = "logistic-normal"
    ),
    class = "link-glm"
  )
}

## Score of the Cantoni-Ronchetti binomial Mqle equation under the integral
## mean. Returns the coefficient-length score vector. wX = 1 when
## robust.wx = "none".
logistic_normal_mqle_score <- function(X, y, beta, v, quad, c_tun, wX = NULL) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  beta <- as.numeric(beta)
  if (is.null(wX)) wX <- rep(1, nrow(X))
  link <- logistic_normal_link(v, quad)
  eta <- as.numeric(X %*% beta)
  mu <- link$linkinv(eta)
  dmu <- link$mu.eta(eta)
  V <- mu * (1 - mu)
  if (any(!is.finite(V)) || any(V <= 0)) {
    stop("The score check encountered a nonpositive binomial variance.")
  }
  r <- (y - mu) / sqrt(V)
  psi <- pmax(-c_tun, pmin(c_tun, r))
  psi1 <- pmax(-c_tun, pmin(c_tun, (1 - mu) / sqrt(V)))
  psi0 <- pmax(-c_tun, pmin(c_tun, -mu / sqrt(V)))
  a_c <- mu * psi1 + (1 - mu) * psi0
  colSums(wX * X * (dmu / sqrt(V)) * (psi - a_c))
}


## ============================================================================
## Local random-intercept variance estimation
## ============================================================================

validate_eta_var <- function(eta.var, S, Tlen) {
  if (is.null(eta.var)) {
    stop("eta.var is missing after variance estimation. ",
         "Use gfmm2d_gaussian() for Gaussian responses.")
  }
  if (is.null(dim(eta.var))) {
    if (length(eta.var) != 1) {
      stop("eta.var must be a scalar or a strict S x T matrix; ",
           "vectors are not recycled into a matrix.")
    }
    eta.var <- matrix(as.numeric(eta.var), nrow = S, ncol = Tlen)
  } else {
    eta.var <- as.matrix(eta.var)
    if (nrow(eta.var) != S || ncol(eta.var) != Tlen) {
      stop("The eta.var matrix must be S x T = ", S, " x ", Tlen,
           ", but it is ", paste(dim(eta.var), collapse = " x "),
           ". Do not transpose the variance to match a T x S coefficient surface.")
    }
  }
  if (any(!is.finite(eta.var)) || any(eta.var < 0)) {
    stop("eta.var must contain finite nonnegative values.")
  }
  eta.var
}

## Build the long-format data used by the local variance model.
## Rows of `data` are subject-major with S rows per subject, matching the
## layout used by gfmm2d_gaussian() and gfmm2d_binomial(). Matrix covariates are indexed by time.
local_long_data <- function(data, formula, S, Tlen) {
  model_formula <- as.character(formula)
  response <- data[, model_formula[2]]
  if (is.null(dim(response))) response <- matrix(response, ncol = Tlen)
  n_block <- nrow(response)
  n <- n_block / S
  if (n != as.integer(n)) stop("The response rows are not divisible by S.")
  n <- as.integer(n)

  id_all <- data$ID
  if (is.null(id_all)) stop("data must contain an ID column.")
  if (length(id_all) == n) id_all <- rep(id_all, each = S)
  if (length(id_all) != n_block) {
    stop("The ID column must have one entry per row of the response.")
  }

  rhs <- model_formula[3]
  rhs_vars <- all.vars(stats::as.formula(paste("~", rhs)))
  S.grid <- seq(0, 1, length.out = S)
  T.grid <- seq(0, 1, length.out = Tlen)

  pieces <- vector("list", Tlen)
  for (tt in seq_len(Tlen)) {
    piece <- data.frame(
      id = id_all,
      s = rep(S.grid, n),
      t = T.grid[tt],
      y = response[, tt],
      stringsAsFactors = FALSE
    )
    for (v in rhs_vars) {
      xv <- data[[v]]
      piece[[v]] <- if (is.matrix(xv) || is.array(xv)) xv[, tt] else xv
    }
    pieces[[tt]] <- piece
  }
  do.call(rbind, pieces)
}

## Compatibility helper: long data for the Y ~ X design used by the smoke
## scripts and the comparison scripts.
long_binomial_xy <- function(dat, S, Tlen) {
  n <- length(unique(dat$ID))
  S.grid <- seq(0, 1, length.out = S)
  T.grid <- seq(0, 1, length.out = Tlen)
  pieces <- vector("list", Tlen)
  for (tt in seq_len(Tlen)) {
    pieces[[tt]] <- data.frame(
      id = dat$ID,
      s = rep(S.grid, n),
      t = T.grid[tt],
      y = dat$Y[, tt],
      x = dat$X,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, pieces)
}

## Fit the local random-intercept model at one target location.
estimate_local_variance <- function(dat, s0, t0, hs, ht, xvars, rhs, nAGQ = 7L) {
  stopifnot(hs > 0, ht > 0)
  required <- c("id", "s", "t", "y", xvars)
  if (!all(required %in% names(dat))) stop("The local data are missing required columns.")
  d <- dat[stats::complete.cases(dat[, required, drop = FALSE]), , drop = FALSE]
  if (!all(d$y %in% c(0, 1))) stop("y must be coded as numeric 0/1.")
  d <- d[abs(d$s - s0) <= hs & abs(d$t - t0) <= ht, , drop = FALSE]
  if (nrow(d) == 0L) stop("No observations in this neighborhood.")
  d$id <- droplevels(factor(d$id))
  if (nlevels(d$id) < 2L) stop("Too few subjects in this neighborhood.")
  if (!any(table(d$id) >= 2L)) stop("No subject has repeated observations in this neighborhood.")
  if (length(unique(d$y)) < 2L) stop("y is constant in this neighborhood.")
  form <- stats::as.formula(paste("y ~", rhs, "+ (1 | id)"))
  ## The singular-fit notice is kept in the returned diagnostics, so only
  ## silence its console message here.
  fit <- suppressMessages(suppressWarnings(lme4::glmer(
    formula = form,
    data = d,
    family = stats::binomial(link = "logit"),
    nAGQ = nAGQ,
    control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )))
  v_raw <- as.numeric(lme4::VarCorr(fit)[["id"]][1, 1])
  messages <- fit@optinfo$conv$lme4$messages
  messages <- if (is.null(messages)) "" else paste(messages, collapse = "; ")
  data.frame(
    s = s0, t = t0, v_raw = v_raw,
    n_subject = nlevels(d$id),
    n_location = nrow(unique(d[, c("s", "t")])),
    n_obs = nrow(d),
    singular = lme4::isSingular(fit, tol = 1e-4),
    message = messages,
    stringsAsFactors = FALSE
  )
}

## Compatibility wrapper for the original one-point interface.
estimate_local_v <- function(dat, s0, t0, hs, ht, xvars = "x", nAGQ = 7) {
  estimate_local_variance(dat, s0, t0, hs, ht,
                          xvars = xvars, rhs = paste(xvars, collapse = " + "),
                          nAGQ = nAGQ)
}

## Fit the local variance at every (s, t) grid point.
local_variance_surface <- function(dat, S, Tlen, hs, ht, nAGQ = 7L,
                                   xvars = "x", rhs = NULL, parallel = TRUE) {
  if (is.null(rhs)) rhs <- paste(xvars, collapse = " + ")
  S.grid <- seq(0, 1, length.out = S)
  T.grid <- seq(0, 1, length.out = Tlen)
  grid <- expand.grid(s = S.grid, t = T.grid, KEEP.OUT.ATTRS = FALSE)
  one <- function(k) {
    tryCatch(
      estimate_local_variance(dat, grid$s[k], grid$t[k], hs, ht, xvars, rhs, nAGQ),
      error = function(e) {
        data.frame(
          s = grid$s[k], t = grid$t[k], v_raw = NA_real_,
          n_subject = NA_integer_, n_location = NA_integer_, n_obs = NA_integer_,
          singular = NA, message = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
  }
  rows <- if (isTRUE(parallel)) {
    parallel::mclapply(seq_len(nrow(grid)), one, mc.cores = max(1, parallel::detectCores() - 1))
  } else {
    lapply(seq_len(nrow(grid)), one)
  }
  v_grid <- do.call(rbind, rows)
  est <- matrix(NA_real_, S, Tlen)
  # expand.grid keeps s varying fastest.
  si <- rep(seq_len(S), times = Tlen)
  ti <- rep(seq_len(Tlen), each = S)
  est[cbind(si, ti)] <- v_grid$v_raw
  list(grid = v_grid, est.var = est)
}

## Estimate the full latent variance surface from the original data.
## window = c(S_window, T_window) gives the number of grid points included
## along the visit and time axes. The default c(1, 3) is the 1 x 3 window:
## only the target visit and one time step on either side. Explicit hs / ht
## override the window-derived half-widths. Failed points are filled with the
## median of the successful points.
estimate_eta_var_local <- function(data, formula, S, Tlen,
                                   window = c(1, 3), hs = NULL, ht = NULL,
                                   nAGQ = 7L, parallel = FALSE, silence = FALSE) {
  window <- as.integer(window)
  if (length(window) != 2 || any(!is.finite(window)) ||
      any(window < 1) || any(window %% 2 == 0)) {
    stop("variance.control$window must contain two positive odd integers ",
         "(the number of grid points along s and t).")
  }
  S.grid <- seq(0, 1, length.out = S)
  T.grid <- seq(0, 1, length.out = Tlen)
  ds <- if (S > 1) S.grid[2] - S.grid[1] else 1
  dt <- if (Tlen > 1) T.grid[2] - T.grid[1] else 1
  if (is.null(hs)) hs <- ((window[1] - 1) / 2) * ds + 1e-8
  if (is.null(ht)) ht <- ((window[2] - 1) / 2) * dt + 1e-8
  if (!is.finite(hs) || hs <= 0 || !is.finite(ht) || ht <= 0) {
    stop("The local variance window half-widths must be positive.")
  }

  rhs <- as.character(formula)[3]
  xvars <- all.vars(stats::as.formula(paste("~", rhs)))
  long <- local_long_data(data, formula, S, Tlen)
  if (!silence) {
    cat(sprintf(
      "  local variance window %d x %d (hs = %.4g, ht = %.4g), nAGQ = %d\n",
      window[1], window[2], hs, ht, as.integer(nAGQ)
    ))
  }
  surf <- local_variance_surface(
    long, S, Tlen, hs, ht, nAGQ = nAGQ,
    xvars = xvars, rhs = rhs, parallel = parallel
  )
  v <- surf$est.var
  ok <- is.finite(v)
  n_fail <- sum(!ok)
  if (n_fail == length(v)) stop("The local variance fit failed at every grid point.")
  if (n_fail > 0) v[!ok] <- stats::median(v[ok])
  summary <- data.frame(
    n_fail = n_fail,
    n_singular = sum(surf$grid$singular %in% TRUE),
    n_zero = sum(v == 0),
    mean_v = mean(v), min_v = min(v), max_v = max(v),
    window_s = window[1], window_t = window[2],
    hs = hs, ht = ht, nAGQ = as.integer(nAGQ),
    stringsAsFactors = FALSE
  )
  list(eta.var = v, grid = surf$grid, summary = summary)
}


## ============================================================================
## Simultaneous bands
## ============================================================================

## Marginal FPCA on the visit domain. Same inputs as utils.R::s.fpca().
s_fpca <- function(eta, S, fpca.opt) {
  S.argvals <- seq_len(S)
  etaList <- vector("list", nrow(eta))
  sList <- vector("list", nrow(eta))
  for (nt in seq_len(nrow(eta))) {
    ind <- which(!is.na(eta[nt, ]))
    etaList[[nt]] <- eta[nt, ind]
    sList[[nt]] <- S.argvals[ind]
  }
  fpca <- suppressWarnings(fdapace::FPCA(etaList, sList, optns = fpca.opt))
  list(pc.s = ncol(fpca$xiEst), xiEst = fpca$xiEst, psi = fpca$phi)
}

## B-spline scores of the FPCA coefficients. Same calculation as utils.R::b.lm().
b_lm <- function(size, Bt, pc.s, xiEst) {
  Tlen <- nrow(Bt)
  ker <- array(0, dim = c(pc.s, ncol(Bt), ncol(Bt)))
  for (j in seq_len(pc.s)) {
    bspl <- matrix(xiEst[, j], nrow = size, ncol = Tlen, byrow = TRUE)
    for (i in seq_len(size)) {
      b <- unname(stats::coef(stats::lm(bspl[i, ] ~ 0 + Bt)))
      ker[j, , ] <- ker[j, , ] + kronecker(t(b), b)
    }
  }
  list(ker = ker)
}

## Original simultaneous-band Monte Carlo. gamma is removed from the standard
## error here, matching the non-analytic branch of gfmm2d_gaussian().
gfmm2d_binomial_scb <- function(betaHat, betaHat.boot, cov.beta.ts.hat, Bt, S, Tlen,
                                gamma, M, bandwidth.control) {
  needed <- c("fdapace", "mvtnorm", "refund")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Simultaneous bands require packages: ", paste(missing, collapse = ", "),
         ". gfmm2d_binomial() does not install packages automatically.")
  }
  B <- length(betaHat.boot)
  if (B < 2) stop("Fewer than 2 valid bootstrap replicates; cannot build simultaneous bands.")
  d <- bandwidth.control$d
  if (is.null(d) || length(d) != 1 || !is.finite(d) || d < 0 || 2 * d >= ncol(Bt)) {
    stop("bandwidth.control$d leaves no basis functions after trimming both ends.")
  }
  fpca.opt <- list(dataType = "Dense", methodSelectK = "FVE")
  P <- length(betaHat)
  qn <- rep(NA_real_, P)
  vm <- matrix(NA_real_, P, M)
  psi.boot <- vector("list", P)
  ker.boot <- vector("list", P)
  beta_avg <- vector("list", P)
  for (p in seq_len(P)) {
    draws <- lapply(betaHat.boot, function(one) one[[p]])
    beta_avg[[p]] <- Reduce("+", draws) / B
    eta <- do.call(rbind, lapply(draws, function(x) x - beta_avg[[p]]))
    eigenS <- s_fpca(eta, S, fpca.opt)
    psi.boot[[p]] <- eigenS$psi
    ker.boot[[p]] <- b_lm(B, Bt, eigenS$pc.s, eigenS$xiEst)$ker
  }
  keep <- (d + 1):(ncol(Bt) - d)
  Bt.r <- Bt[, keep, drop = FALSE]
  for (p in seq_len(P)) {
    sp <- sqrt(diag(cov.beta.ts.hat[, , p]) / gamma)
    qm <- lapply(seq_len(M), function(m) array(0, dim = c(Tlen, S)))
    for (j in seq_len(dim(ker.boot[[p]])[1])) {
      scale <- if (is.null(bandwidth.control$adjust)) 3.2 else bandwidth.control$adjust
      Sig <- ker.boot[[p]][j, , ] / (B * scale)
      Sig.r.t <- Sig[keep, keep, drop = FALSE]
      Sig.r <- refund::fbps(Sig.r.t, knots = round(ncol(Bt.r) / 1.5))$Yhat
      sigma <- if (isSymmetric(Sig.r)) Sig.r else Sig.r.t
      umat <- mvtnorm::rmvnorm(M, mean = rep(0, ncol(Bt.r)), sigma = sigma)
      Btu <- apply(umat, 1, function(x) Bt.r %*% x)
      for (m in seq_len(M)) {
        qm[[m]] <- qm[[m]] + kronecker(Btu[, m], t(psi.boot[[p]][, j]))
      }
    }
    for (m in seq_len(M)) {
      draw <- qm[[m]] + beta_avg[[p]]
      vm[p, m] <- max(abs(draw - betaHat[[p]]) / array(sp, dim = c(Tlen, S)))
    }
    qn[p] <- unname(stats::quantile(vm[p, ], 0.95))
  }
  names(qn) <- names(betaHat)
  list(qn = qn, vm = vm)
}


## ============================================================================
## Bootstrap
## ============================================================================

fmm2d_binomial_boot <- function(b, data, n, formula, S, robust, smoother, knots,
                                eta.var, gh.n, robust.c, robust.wx, parallel) {
  if (isTRUE(parallel)) {
    stop("Pointwise fits inside a bootstrap replicate stay serial.")
  }
  if (b %% 10 == 0) cat(b, "\n")
  ids <- unique(data$ID)
  if (length(ids) != n) stop("The number of unique IDs does not match n.")
  blocks <- split(seq_len(nrow(data)), data$ID)
  for (try_i in seq_len(50)) {
    drawn <- sample(ids, size = n, replace = TRUE)
    key <- as.character(drawn)
    row.ind <- unlist(blocks[key], use.names = FALSE)
    data.boot <- data[row.ind, , drop = FALSE]
    data.boot$ID <- rep(seq_along(drawn), times = lengths(blocks[key]))
    fit <- try(gfmm2d_binomial(
      formula = formula, data = data.boot, S = S, eta.var = eta.var,
      robust = robust, smoother = smoother, knots = knots, gh.n = gh.n,
      robust.c = robust.c, robust.wx = robust.wx, pcb = FALSE, scb = FALSE,
      parallel = FALSE, silence = TRUE, keep_diagnostics = FALSE
    ), silent = TRUE)
    if (!inherits(fit, "try-error") && all(vapply(fit$betaHat, function(x) all(is.finite(x)), logical(1)))) {
      return(fit$betaHat)
    }
  }
  stop("A bootstrap replicate failed to produce finite smoothed coefficients 50 times.")
}


## ============================================================================
## Timing
## ============================================================================

gfmm2d_binomial_elapsed <- function(t_start, t_var, t_fit, t_boot, t_end = proc.time()) {
  elapsed <- function(a, b) unname((b - a)["elapsed"])
  c(
    variance = elapsed(t_start, t_var),
    fit = elapsed(t_var, t_fit),
    bootstrap = elapsed(t_fit, t_boot),
    bands = elapsed(t_boot, t_end),
    total = elapsed(t_start, t_end)
  )
}

