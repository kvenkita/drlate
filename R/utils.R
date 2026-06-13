# Internal utilities: weighted statistics, validation, covariate standardization.

#' Weighted mean (Stata `summarize [iw=w]` r(mean))
#' @noRd
wmean <- function(x, w) sum(w * x) / sum(w)

#' Weighted standard deviation (Stata `summarize [iw=w]` r(sd):
#' denominator sum(w) - 1)
#' @noRd
wsd <- function(x, w) {
  m <- wmean(x, w)
  sqrt(sum(w * (x - m)^2) / (sum(w) - 1))
}

#' Weighted variance (square of [wsd]); reduces to [stats::var] when w = 1.
#' @noRd
wvar <- function(x, w) wsd(x, w)^2

#' Check that a variable takes exactly the values 0 and 1
#' @noRd
check_binary <- function(x, name, role) {
  ux <- unique(x)
  if (length(ux) != 2L || !all(ux %in% c(0, 1))) {
    stop(sprintf("%s `%s` must be binary 0/1.", role, name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate a response against its model family, mirroring drlate_estimate.ado
#' @noRd
check_family <- function(x, name, family, role) {
  if (family == "binomial") check_binary(x, name, role)
  if (family == "poisson" && any(x < 0)) {
    stop(sprintf("%s `%s` must be non-negative for the poisson model.",
                 name, role), call. = FALSE)
  }
  invisible(TRUE)
}

#' Standardize continuous model-matrix columns in place.
#'
#' Mirrors the numerical-stability standardization in drlate_estimate.ado:
#' columns with more than two distinct values are centered and scaled by
#' their weighted mean and sd. Since the intercept is always included,
#' this never changes the column span, so fitted values (and therefore the
#' reported estimates) are unaffected; it only conditions the optimization
#' and the numeric Jacobian.
#' @noRd
standardize_mm <- function(X, w) {
  for (j in seq_len(ncol(X))) {
    cj <- X[, j]
    if (length(unique(cj)) > 2L) {
      s <- wsd(cj, w)
      if (is.finite(s) && s > 0) X[, j] <- (cj - wmean(cj, w)) / s
    }
  }
  X
}

#' Map drlate model strings to family handles used internally
#' @noRd
fam_linkinv <- function(family) {
  switch(family,
    gaussian = identity,
    binomial = stats::plogis,
    poisson  = exp,
    stop("unknown family: ", family)
  )
}

#' Overlap check, mirroring _drlate_ps.ado
#' @noRd
check_overlap <- function(ps, pstolerance, osample) {
  viol <- ps < pstolerance | ps > (1 - pstolerance)
  nviol <- sum(viol)
  if (nviol > 0 && !osample) {
    stop(sprintf(paste0(
      "%d observation%s violate the overlap assumption (instrument ",
      "propensity score outside [%g, 1-%g]). Re-run with `osample = TRUE` ",
      "to identify the violating observations."),
      nviol, if (nviol > 1) "s" else "", pstolerance, pstolerance),
      call. = FALSE)
  }
  viol
}
