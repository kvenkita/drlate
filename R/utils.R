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

#' User-facing outcome/treatment model names mapped to the internal link
#' token, the required response domain, and a printed label. The estimation
#' machinery is link-driven, so the fractional families (flogit, fprobit)
#' share all fit/predict/score code with their binary counterparts (logit,
#' probit) and differ only in the allowed response domain and the label.
#' @noRd
.drlate_models <- list(
  linear  = list(link = "gaussian", domain = "real",   label = "linear"),
  logit   = list(link = "logit",    domain = "binary", label = "logit"),
  probit  = list(link = "probit",   domain = "binary", label = "probit"),
  poisson = list(link = "poisson",  domain = "nonneg", label = "poisson"),
  flogit  = list(link = "logit",    domain = "unit",
                 label = "fractional logit"),
  fprobit = list(link = "probit",   domain = "unit",
                 label = "fractional probit")
)

#' Internal link token for a user-facing model name
#' @noRd
model_link <- function(model) .drlate_models[[model]]$link

#' Printed label for a user-facing model name
#' @noRd
model_label <- function(model) .drlate_models[[model]]$label

#' Check that a variable takes exactly the values 0 and 1
#' @noRd
check_binary <- function(x, name, role) {
  ux <- unique(x)
  if (length(ux) != 2L || !all(ux %in% c(0, 1))) {
    stop(sprintf("%s `%s` must be binary 0/1.", role, name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Check that a variable lies in the unit interval [0, 1] (fractional models)
#' @noRd
check_fractional <- function(x, name, role) {
  if (any(x < 0 | x > 1)) {
    stop(sprintf("%s `%s` must lie in [0, 1] for a fractional model.",
                 role, name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate a response against its (user-facing) model family, mirroring
#' drlate_estimate.ado
#' @noRd
check_family <- function(x, name, model, role) {
  domain <- .drlate_models[[model]]$domain
  if (domain == "binary") {
    check_binary(x, name, role)
  } else if (domain == "unit") {
    check_fractional(x, name, role)
  } else if (domain == "nonneg" && any(x < 0)) {
    stop(sprintf("%s `%s` must be non-negative for the %s model.",
                 role, name, model), call. = FALSE)
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

#' Inverse-link for an internal link token
#' @noRd
fam_linkinv <- function(link) {
  switch(link,
    gaussian = identity,
    logit    = stats::plogis,
    probit   = stats::pnorm,
    poisson  = exp,
    stop("unknown link: ", link)
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
