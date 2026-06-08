# Nonparametric bootstrap for the (late, num, denom) point estimates.
# Each draw re-runs the sequential point estimation ONLY (no sandwich),
# via the same compute_point() code path as the main fit.

#' Resample a ctx by row indices, recomputing compliance means and case
#' @noRd
ctx_resample <- function(ctx, idx) {
  b <- ctx
  b$y <- ctx$y[idx]
  b$d <- ctx$d[idx]
  b$z <- ctx$z[idx]
  b$w <- ctx$w[idx]
  b$cluster <- if (is.null(ctx$cluster)) NULL else ctx$cluster[idx]
  b$Xo <- ctx$Xo[idx, , drop = FALSE]
  b$Xt <- ctx$Xt[idx, , drop = FALSE]
  b$Xz <- ctx$Xz[idx, , drop = FALSE]
  b$n <- length(idx)
  b$dmeanz1 <- mean(b$d[b$z == 1])
  b$dmeanz0 <- mean(b$d[b$z == 0])
  b$case <- if (b$dmeanz0 %in% c(0, 1) && b$dmeanz1 %in% c(0, 1)) "bothdeg"
            else if (b$dmeanz0 %in% c(0, 1)) "z0deg"
            else if (b$dmeanz1 %in% c(0, 1)) "z1deg"
            else "interior"
  b$osample <- FALSE   # overlap violations inside a draw raise (and skip)
  b
}

#' One bootstrap draw: resample, re-estimate, return the method's estimates vector
#' or NAs on any failure (degenerate resample, non-convergence, overlap).
#' @noRd
boot_draw <- function(ctx, units, unit_rows, stat_names) {
  take <- sample(length(units), replace = TRUE)
  idx <- unlist(unit_rows[take], use.names = FALSE)
  # Suppress per-draw messages/warnings (e.g. GLM non-convergence notes);
  # genuine failures become NA rows and are counted by the caller.
  tryCatch(
    suppressWarnings(suppressMessages(
      compute_point(ctx_resample(ctx, idx))$est$estimates)),
    error = function(e) stats::setNames(rep(NA_real_, length(stat_names)),
                                        stat_names)
  )
}

#' Run the bootstrap; returns draws matrix, SEs, percentile CIs, counts.
#' @noRd
drlate_boot <- function(ctx, reps, seed = NULL, cores = 1L, level = 0.95,
                        stat_names = c("late", "num", "denom")) {
  # Resampling units: whole clusters when clustering, else rows
  if (is.null(ctx$cluster)) {
    units <- seq_len(ctx$n)
    unit_rows <- as.list(units)
  } else {
    units <- unique(ctx$cluster)
    unit_rows <- split(seq_len(ctx$n), ctx$cluster)[as.character(units)]
    if (length(units) < 30L) {
      warning("cluster bootstrap with only ", length(units), " clusters: ",
              "the resampling distribution is coarse and percentile ",
              "intervals may be unreliable; consider the analytic ",
              "cluster-robust variance instead.", call. = FALSE)
    }
  }

  if (!is.null(seed)) {
    old_kind <- RNGkind("L'Ecuyer-CMRG")
    on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
    set.seed(seed)
  }

  if (cores > 1L) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
    draws <- parallel::parLapply(cl, seq_len(reps), function(i, ctx, units,
                                                            unit_rows,
                                                            stat_names) {
      boot_draw(ctx, units, unit_rows, stat_names)
    }, ctx = ctx, units = units, unit_rows = unit_rows, stat_names = stat_names)
  } else {
    draws <- lapply(seq_len(reps), function(i)
      boot_draw(ctx, units, unit_rows, stat_names))
  }
  draws <- do.call(rbind, draws)
  colnames(draws) <- stat_names

  ok <- stats::complete.cases(draws)
  n_fail <- sum(!ok)
  if (n_fail > 0 && n_fail / reps > 0.01) {
    warning(n_fail, " of ", reps, " bootstrap draws failed (degenerate ",
            "resample, non-convergence, or overlap violation) and were ",
            "dropped.", call. = FALSE)
  }
  good <- draws[ok, , drop = FALSE]
  if (nrow(good) < 2L) {
    stop("the bootstrap failed in nearly all draws; check the sample size ",
         "and overlap.", call. = FALSE)
  }

  a <- (1 - level) / 2
  list(
    draws = good,
    se = apply(good, 2, stats::sd),
    ci = t(apply(good, 2, stats::quantile, probs = c(a, 1 - a))),
    reps = reps,
    reps_ok = nrow(good)
  )
}
