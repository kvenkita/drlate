# Sequential weighted GLM fits producing the analytic point estimates
# (Stata's `from()` starting values, which ARE the reported estimates
# because the stacked GMM runs with iterate(0)).

#' Weighted GLM on a subsample; returns coefficients named by columns of X.
#' quasi- families allow non-integer weights with identical score equations
#' (hence identical coefficients) to binomial/poisson MLE.
#' @noRd
fit_wglm <- function(yvar, X, family, w, subset) {
  ys <- yvar[subset]
  Xs <- X[subset, , drop = FALSE]
  ws <- w[subset]
  if (family == "gaussian") {
    fit <- stats::lm.wfit(Xs, ys, ws)
    b <- fit$coefficients
  } else {
    fam <- if (family == "binomial") stats::quasibinomial()
           else stats::quasipoisson()
    fit <- suppressWarnings(
      stats::glm.fit(Xs, ys, weights = ws, family = fam)
    )
    if (!fit$converged) {
      warning("GLM did not converge; coefficients used as-is ",
              "(mirrors Stata's behavior).", call. = FALSE)
    }
    b <- fit$coefficients
  }
  if (anyNA(b)) {
    stop("collinear columns in a model matrix; remove redundant covariates.",
         call. = FALSE)
  }
  stats::setNames(b, colnames(X))
}

#' Out-of-sample predictions on the response scale (Stata `predict` after a
#' subsample fit predicts for all observations).
#' @noRd
predict_glm <- function(b, X, family) {
  fam_linkinv(family)(drop(X %*% b))
}

#' Fit an instrument-arm regression unless that arm is degenerate.
#' Returns list(coefs = named numeric or NULL, mu = n-vector predictions,
#' degenerate_value = NULL or the constant).
#' @noRd
fit_arm <- function(ctx, yvar, X, family, arm, fitw, dmean) {
  if (dmean %in% c(0, 1)) {
    return(list(coefs = NULL, mu = rep(dmean, ctx$n), degenerate_value = dmean))
  }
  subset <- if (arm == 1) ctx$z == 1 else ctx$z == 0
  b <- fit_wglm(yvar, X, family, fitw, subset)
  list(coefs = b, mu = predict_glm(b, X, family), degenerate_value = NULL)
}
