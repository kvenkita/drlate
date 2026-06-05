# Diagnostic plots for drlate fits (ggplot2, in Suggests).

#' Diagnostic plots for drlate fits
#'
#' @param x A fitted [drlate()] object (with `keep_data = TRUE`).
#' @param type One of:
#'   * `"overlap"`: histograms of the estimated instrument propensity score
#'     by instrument arm, with the `pstolerance` bounds marked. Mass piling
#'     up near 0 or 1 signals overlap problems.
#'   * `"balance"`: a love plot of standardized mean differences from
#'     [balance()], unweighted vs IPW-weighted, with the conventional
#'     |SMD| = 0.1 reference lines.
#'   * `"weights"`: distributions of the implied IPW weights by arm;
#'     a long right tail means a few observations dominate the estimate.
#' @param bins Number of histogram bins for `"overlap"` and `"weights"`.
#' @param ... Currently unused.
#' @return A `ggplot` object.
#' @export
plot.drlate <- function(x, type = c("overlap", "balance", "weights"),
                        bins = 30, ...) {
  type <- match.arg(type)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("package 'ggplot2' is required for drlate plots; install it with ",
         "install.packages(\"ggplot2\").", call. = FALSE)
  }
  switch(type,
    overlap = plot_overlap(x, bins),
    balance = plot_balance(x),
    weights = plot_weights(x, bins)
  )
}

#' @noRd
need_ps <- function(x) {
  if (is.null(x$ps)) {
    stop("no instrument propensity score is estimated with method = \"ra\"; ",
         "this plot requires an IPW-type fit.", call. = FALSE)
  }
  x$ps
}

#' @noRd
plot_overlap <- function(x, bins) {
  ctx <- need_ctx(x)
  ps <- need_ps(x)
  df <- data.frame(
    ps = ps$ps,
    arm = factor(ifelse(ctx$z == 1, "Z = 1", "Z = 0"),
                 levels = c("Z = 0", "Z = 1"))
  )
  ggplot2::ggplot(df, ggplot2::aes(x = .data$ps, fill = .data$arm)) +
    ggplot2::geom_histogram(bins = bins, alpha = 0.55,
                            position = "identity", color = NA) +
    ggplot2::geom_vline(xintercept = c(ctx$pstolerance,
                                       1 - ctx$pstolerance),
                        linetype = "dashed", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = c("Z = 0" = "#9e9e9e",
                                          "Z = 1" = "#2c6e8f")) +
    ggplot2::labs(
      x = "Estimated instrument propensity score",
      y = "Count", fill = NULL,
      title = "Overlap of the instrument propensity score",
      subtitle = "Mass near 0 or 1 indicates limited overlap"
    ) +
    ggplot2::theme_minimal()
}

#' @noRd
plot_balance <- function(x) {
  b <- balance(x)
  long <- data.frame(
    variable = rep(b$variable, 2L),
    smd = c(b$smd_unweighted, b$smd_weighted),
    weighting = rep(c("Unweighted", "IPW-weighted"), each = nrow(b))
  )
  long$weighting <- factor(long$weighting,
                           levels = c("Unweighted", "IPW-weighted"))
  # Order variables by absolute unweighted imbalance
  ord <- b$variable[order(abs(b$smd_unweighted))]
  long$variable <- factor(long$variable, levels = ord)
  ggplot2::ggplot(long, ggplot2::aes(x = .data$smd, y = .data$variable,
                                     color = .data$weighting,
                                     shape = .data$weighting)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed",
                        linewidth = 0.3) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_color_manual(values = c("Unweighted" = "#9e9e9e",
                                           "IPW-weighted" = "#2c6e8f")) +
    ggplot2::labs(
      x = "Standardized mean difference (Z = 1 vs Z = 0)",
      y = NULL, color = NULL, shape = NULL,
      title = "Covariate balance across instrument arms",
      subtitle = "Dashed lines mark the conventional |SMD| = 0.1 threshold"
    ) +
    ggplot2::theme_minimal()
}

#' @noRd
plot_weights <- function(x, bins) {
  ctx <- need_ctx(x)
  ps <- need_ps(x)
  if (x$estimand == "latt") {
    wts <- ifelse(ctx$z == 1, ctx$w, ctx$w * ps$ps / (1 - ps$ps))
  } else {
    wts <- ifelse(ctx$z == 1, ctx$w / ps$ps, ctx$w / (1 - ps$ps))
  }
  df <- data.frame(
    weight = wts,
    arm = factor(ifelse(ctx$z == 1, "Z = 1", "Z = 0"),
                 levels = c("Z = 0", "Z = 1"))
  )
  q99 <- stats::quantile(df$weight, 0.99)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$weight, fill = .data$arm)) +
    ggplot2::geom_histogram(bins = bins, alpha = 0.55,
                            position = "identity", color = NA) +
    ggplot2::geom_vline(xintercept = q99, linetype = "dotted",
                        linewidth = 0.3) +
    ggplot2::facet_wrap(~arm, ncol = 1L, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c("Z = 0" = "#9e9e9e",
                                          "Z = 1" = "#2c6e8f"),
                               guide = "none") +
    ggplot2::labs(
      x = "Inverse-propensity weight component", y = "Count",
      title = "Distribution of implied weights",
      subtitle = "Dotted line marks the 99th percentile; long tails mean a few observations dominate"
    ) +
    ggplot2::theme_minimal()
}
