# Covariate balance diagnostics.

#' Covariate balance across instrument arms
#'
#' Computes standardized mean differences (SMDs) of the model covariates
#' between the two instrument arms, before and after weighting by the
#' inverse of the estimated instrument propensity score. Well-balanced
#' weighted covariates (conventionally, absolute SMD below 0.1) indicate
#' that the propensity score model is doing its job.
#'
#' The covariate set is the union of the columns of the instrument,
#' outcome, and treatment model matrices (the intercept is dropped). The
#' SMD denominator is the unweighted pooled standard deviation
#' \eqn{\sqrt{(s_1^2 + s_0^2)/2}} in both columns, so the two columns are
#' directly comparable. Weighted arm means are Hájek means using the
#' inverse-propensity weights implied by the fit (for
#' `estimand = "latt"`, the Z=0 arm uses the ATT odds weights
#' \eqn{p/(1-p)}, matching the estimator).
#'
#' @param object A fitted [drlate()] object (with `keep_data = TRUE`).
#' @param ... Currently unused.
#' @return A data frame with one row per covariate and columns
#'   `variable`, `smd_unweighted`, and `smd_weighted`.
#' @seealso [plot.drlate()] with `type = "balance"` for the love plot.
#' @export
balance <- function(object, ...) UseMethod("balance")

#' @rdname balance
#' @export
balance.drlate <- function(object, ...) {
  ctx <- need_ctx(object)
  if (is.null(object$ps)) {
    stop("no instrument propensity score is estimated with method = \"ra\"; ",
         "balance() requires an IPW-type fit.", call. = FALSE)
  }
  X <- diag_covariates(ctx)
  z <- ctx$z
  w <- ctx$w
  ps <- object$ps$ps

  # Weights per arm: LATE reweights both arms by inverse propensity;
  # LATT keeps the Z=1 arm as-is and odds-weights the Z=0 arm.
  if (object$estimand == "latt") {
    w1 <- w
    w0 <- w * ps / (1 - ps)
  } else {
    w1 <- w / ps
    w0 <- w / (1 - ps)
  }

  smd <- function(x, wts1, wts0) {
    s1 <- z == 1
    pooled <- sqrt((stats::var(x[s1]) + stats::var(x[!s1])) / 2)
    if (!is.finite(pooled) || pooled == 0) return(NA_real_)
    (wmean(x[s1], wts1[s1]) - wmean(x[!s1], wts0[!s1])) / pooled
  }

  data.frame(
    variable = colnames(X),
    smd_unweighted = apply(X, 2, smd, wts1 = w, wts0 = w),
    smd_weighted   = apply(X, 2, smd, wts1 = w1, wts0 = w0),
    row.names = NULL
  )
}

#' Deduplicated union of the three model matrices, intercept dropped
#' @noRd
diag_covariates <- function(ctx) {
  X <- cbind(ctx$Xz, ctx$Xo, ctx$Xt)
  X <- X[, !duplicated(colnames(X)), drop = FALSE]
  X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
}

#' Fetch the retained ctx or fail with guidance
#' @noRd
need_ctx <- function(object) {
  if (is.null(object$ctx)) {
    stop("this fit was created with keep_data = FALSE; re-fit with ",
         "keep_data = TRUE to use diagnostics.", call. = FALSE)
  }
  object$ctx
}

#' First-stage z-statistic (effect of Z on D over its SE)
#' @noRd
firststage_z <- function(x) {
  unname(x$coefficients[3] / sqrt(x$vcov3[3, 3]))
}
