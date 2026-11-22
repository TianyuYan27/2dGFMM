## Shared plotting helpers: bivariate and marginal confidence bands, and
## panel assembly. Nothing in this file is specific to a data application.

## Heat map of one bivariate coefficient surface. beta_hat[[p]] is indexed
## by t (rows) and s (columns). Tiles are coloured by the sign of the
## estimate where the band excludes zero and are left near-white elsewhere.
plot_bivariate_band <- function(beta_hat, beta_cov, p, qn, scb = TRUE, limits = NULL) {
  se <- sqrt(diag(beta_cov[, , p]))
  multiplier <- if (scb) qn[p] else 2
  beta_hat_up <- as.vector(beta_hat[[p]]) + multiplier * se
  beta_hat_low <- as.vector(beta_hat[[p]]) - multiplier * se
  n_t <- nrow(beta_hat[[p]])
  n_s <- ncol(beta_hat[[p]])
  grid <- data.frame(
    betaHat = as.vector(beta_hat[[p]]),
    betaHat.low = beta_hat_low,
    betaHat.up = beta_hat_up,
    t = rep(seq_len(n_t), times = n_s),
    s = rep(seq_len(n_s), each = n_t)
  )
  grid$sign <- ifelse(grid$betaHat.up < 0, grid$betaHat,
                      ifelse(grid$betaHat.low > 0, grid$betaHat, 0))
  if (is.null(limits)) {
    limits <- range(grid$betaHat)
  }
  ggplot(grid, aes(t, s, fill = sign)) +
    geom_tile() +
    theme_bw() +
    theme(axis.text.x = element_text(size = 13), axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 13), axis.title.y = element_text(size = 16),
          legend.title = element_text(size = 14), legend.key.width = unit(1, "cm"),
          plot.title = element_text(size = 19, hjust = 0.5)) +
    scale_fill_gradient2(low = "Steel Blue 3", mid = "grey98",
                         high = "red", midpoint = 0, oob = scales::squish,
                         limits = limits, space = "Lab",
                         name = expression(paste(hat(I), "(s,t)")))
}

## Marginal band of one coefficient surface, averaged over either axis.
## side = "t" averages over s and returns a function of t; side = "s"
## averages over t and returns a function of s.
plot_marginal_band <- function(beta_hat, beta_cov, p, qn, side,
                               zero.line = TRUE, pcb = TRUE, scb = TRUE) {
  n_t <- nrow(beta_hat[[p]])
  n_s <- ncol(beta_hat[[p]])
  if (side == "t") {
    beta_marginal <- apply(beta_hat[[p]], 1, mean)
    cov_marginal <- matrix(0, n_t, n_t)
    for (t in seq_len(n_t)) {
      for (v in t:n_t) {
        cov_tmp <- 0
        for (s in seq_len(n_s)) {
          for (u in seq_len(n_s)) {
            cov_tmp <- beta_cov[(s - 1) * n_t + t, (u - 1) * n_t + v, p] + cov_tmp
          }
        }
        cov_marginal[v, t] <- cov_marginal[t, v] <- cov_tmp / (n_s^2)
      }
    }
    side_len <- n_t
  } else if (side == "s") {
    beta_marginal <- apply(beta_hat[[p]], 2, mean)
    cov_marginal <- matrix(0, n_s, n_s)
    for (s in seq_len(n_s)) {
      for (u in s:n_s) {
        cov_tmp <- 0
        for (t in seq_len(n_t)) {
          for (v in seq_len(n_t)) {
            cov_tmp <- beta_cov[(s - 1) * n_t + t, (u - 1) * n_t + v, p] + cov_tmp
          }
        }
        cov_marginal[u, s] <- cov_marginal[s, u] <- cov_tmp / (n_t^2)
      }
    }
    side_len <- n_s
  } else {
    stop("side must be \"t\" or \"s\".")
  }
  marginal <- data.frame(
    x = seq_len(side_len),
    beta = beta_marginal,
    pcb.low = beta_marginal - 2 * sqrt(diag(cov_marginal)),
    pcb.up = beta_marginal + 2 * sqrt(diag(cov_marginal)),
    scb.low = beta_marginal - qn[p] * sqrt(diag(cov_marginal)),
    scb.up = beta_marginal + qn[p] * sqrt(diag(cov_marginal))
  )
  ggplot() +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    {if (scb) geom_ribbon(aes(x = x, ymax = scb.up, ymin = scb.low),
                         data = marginal, fill = "gray10", alpha = 0.4)} +
    {if (pcb) geom_ribbon(aes(x = x, ymax = pcb.up, ymin = pcb.low),
                         data = marginal, fill = "gray30", alpha = 0.4)} +
    geom_line(aes(x = x, y = beta, color = "Estimate"),
              data = marginal, alpha = 1, lty = 1, lwd = 0.8) +
    scale_colour_manual(name = "", values = c(Estimate = "blue3")) +
    {if (zero.line) geom_hline(yintercept = 0, color = "black", lwd = 0.8, linetype = "dotted")} +
    theme(axis.text.x = element_text(size = 13), axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 13), axis.title.y = element_text(size = 16),
          legend.title = element_text(size = 14), legend.key.width = unit(1, "cm"),
          plot.title = element_text(size = 18, hjust = 0.5))
}

## Arrange ggplot panels into one JPEG.
save_panels <- function(panels, file, ncol, nrow, width, height, legend, common_legend = FALSE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  jpeg(file, units = "in", width = width, height = height, res = 300)
  on.exit(dev.off(), add = TRUE)
  print(ggarrange(plotlist = panels, ncol = ncol, nrow = nrow,
                  legend = legend, common.legend = common_legend))
}

## One heat-map panel per coefficient surface.
bivariate_panels <- function(fit, titles, x, y, x_breaks, x_labels, y_breaks, y_labels) {
  beta_hat <- fit$betaHat
  beta_cov <- fit$betaHat.cov
  qn <- fit$qn
  lapply(seq_along(beta_hat), function(p) {
    plot_bivariate_band(beta_hat, beta_cov, p, qn = qn, scb = FALSE) +
      labs(y = y, x = x, title = titles[p]) +
      scale_y_continuous(breaks = y_breaks, labels = y_labels) +
      scale_x_continuous(breaks = x_breaks, labels = x_labels)
  })
}

## One marginal-band panel per coefficient surface.
marginal_panels <- function(fit, titles, side, x, x_breaks, x_labels) {
  beta_hat <- fit$betaHat
  beta_cov <- fit$betaHat.cov
  qn <- fit$qn
  lapply(seq_along(beta_hat), function(p) {
    plot_marginal_band(beta_hat, beta_cov, p, qn = qn, side = side, zero.line = FALSE) +
      labs(y = expression(hat(beta)(t)), x = x, title = titles[p]) +
      scale_x_continuous(breaks = x_breaks, labels = x_labels)
  })
}

## Wrap a JPEG file in a ggplot panel, so it can join a ggarrange grid.
jpeg_panel <- function(file, title) {
  img <- jpeg::readJPEG(file)
  ggplot() +
    annotation_raster(img, 0, 1, 0, 1) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(title = title) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
}
