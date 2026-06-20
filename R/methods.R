# S3 methods for class "drlate".

#' @export
print.drlate <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  title <- if (x$estimand == "late") "Local average treatment effect"
           else "Local average treatment effect on the treated"
  methodd <- toupper(x$method)
  normd <- if (x$ivmodel == "ipt" || x$statnorm == "nrm") "normalized"
           else "unnormalized"
  kappa_names <- c(kappa = "tau_a", kappa0 = "tau_a,0", kappa10 = "tau_a,10")
  est_line <- if (x$method %in% c("ra", "ipwra")) {
    methodd
  } else if (x$method %in% names(kappa_names)) {
    kind <- if (x$method == "kappa10") "normalized" else "unnormalized"
    paste0(methodd, " (", kappa_names[[x$method]], "; ", kind,
           " Abadie kappa weighting)")
  } else if (x$method == "ipw" && x$ivmodel != "ipt") {
    paste0(methodd, " (", normd, "; kappalate ",
           if (normd == "normalized") "tau_u" else "tau_a,1", ")")
  } else {
    paste0(methodd, " (", normd, ")")
  }
  omodeld <- model_label(x$omodel)
  tmodeld <- model_label(x$tmodel)
  if (x$method == "ipw") {
    omodeld <- "weighted mean"
    tmodeld <- "weighted mean"
  }
  if (x$method %in% names(kappa_names)) {
    omodeld <- "none (kappa weighting)"
    tmodeld <- "none (kappa weighting)"
  }
  zmodeld <- switch(x$ivmodel,
    logit = "logit (MLE)", cbps = "logit (CBPS)", ipt = "logit (IPT)",
    probit = "probit (MLE)")

  cat("\n", title, "\n", sep = "")
  cat("Number of obs    : ", format(x$N, big.mark = ","), "\n", sep = "")
  if (!is.null(x$N_clust)) {
    cat("Number of clusters: ", x$N_clust, "\n", sep = "")
  }
  cat("Estimator        : ", est_line, "\n", sep = "")
  cat("Outcome model    : ", omodeld, "\n", sep = "")
  cat("Treatment model  : ", tmodeld, "\n", sep = "")
  cat("Instrument model : ", zmodeld, "\n", sep = "")
  if (identical(x$vcov_method, "bootstrap")) {
    cat("Std. errors      : nonparametric bootstrap (",
        x$boot$reps_ok, " of ", x$boot$reps, " reps)\n", sep = "")
  }
  cat("\n")

  print(coeftable(x), digits = digits)

  fz <- firststage_z(x)
  # With a single binary instrument, z^2 is the first-stage robust F;
  # flag below the conventional F = 10 (|z| ~ 3.16) and show the
  # weak-instrument-robust Fieller set alongside the Wald interval.
  cat("\nFirst stage (Z on D): z = ", format(fz, digits = 4),
      " (z^2 ~ first-stage F = ", format(fz^2, digits = 3), ")",
      sep = "")
  if (is.finite(fz) && fz^2 < 10) {
    cat("  [weak: Wald inference on the ratio may be unreliable]")
    if (!is.null(x$layout$num) && !is.null(x$layout$denom)) {
      f <- fieller_from_fit(x)
      cat("\nFieller 95% confidence set for the ", x$estimand, ": ",
          format_fieller(f, digits = max(4L, digits)), sep = "")
    }
  }
  cat("\n")
  invisible(x)
}

#' @noRd
coeftable <- function(x) {
  b <- x$coefficients
  if (identical(x$vcov_method, "bootstrap")) {
    se <- unname(x$boot$se)
    ci <- unname(x$boot$ci)
    ci_labels <- c("[2.5% boot", "boot 97.5%]")
  } else {
    se <- sqrt(diag(x$vcov3))
    ci <- cbind(b - stats::qnorm(0.975) * se, b + stats::qnorm(0.975) * se)
    ci_labels <- c("[95% conf.", "interval]")
  }
  zstat <- b / se
  p <- 2 * stats::pnorm(-abs(zstat))
  out <- cbind(b, se, zstat, p, ci)
  dimnames(out) <- list(names(b),
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)", ci_labels))
  out
}

#' @export
summary.drlate <- function(object, ...) {
  structure(list(fit = object, coeftable = coeftable(object)),
            class = "summary.drlate")
}

#' @export
print.summary.drlate <- function(x,
                                 digits = max(3L, getOption("digits") - 3L),
                                 ...) {
  print(x$fit, digits = digits, ...)
  invisible(x)
}

#' @export
coef.drlate <- function(object, ...) object$coefficients

#' @export
vcov.drlate <- function(object, full = FALSE, ...) {
  if (full) object$vcov_full else object$vcov3
}

#' @export
nobs.drlate <- function(object, ...) object$N

#' Confidence intervals for drlate fits
#'
#' @param object A fitted [drlate()] object.
#' @param parm Coefficients to include (names or indices); defaults to all
#'   three reported quantities.
#' @param level Confidence level.
#' @param method `"default"` gives Wald intervals from the joint sandwich
#'   (or bootstrap percentile intervals when the fit used
#'   `vcov = "bootstrap"`). `"fieller"` inverts the test of
#'   `num - t * denom = 0` using the joint covariance of the numerator and
#'   denominator, giving a confidence set for the LATE/LATT ratio that
#'   remains valid when the first stage is weak; the set may be an
#'   interval, the complement of an interval, or the whole line, and is
#'   returned as a `"drlate_fieller"` object with its own print method.
#' @param ... Currently unused.
#' @return
#'   For `method = "default"`, a numeric matrix with one row per requested
#'   coefficient (`parm`) and two columns holding the lower and upper
#'   confidence limits. The columns are labelled with the corresponding
#'   percentiles (for the default 95% level, `"2.5 %"` and `"97.5 %"`). The
#'   limits are Wald intervals from the joint sandwich covariance, or
#'   percentile intervals from the resampling draws when the fit was computed
#'   with `vcov = "bootstrap"`.
#'
#'   For `method = "fieller"`, an object of class `"drlate_fieller"`: a list
#'   describing the weak-instrument-robust confidence set for the LATE/LATT
#'   ratio (its endpoints and shape, the estimand name, and the confidence
#'   level), with its own `print` method. Because a Fieller set need not be a
#'   bounded interval, it is returned in this form rather than as a matrix of
#'   endpoints.
#' @export
confint.drlate <- function(object, parm, level = 0.95,
                           method = c("default", "fieller"), ...) {
  method <- match.arg(method)
  b <- object$coefficients
  a <- (1 - level) / 2

  if (method == "fieller") {
    f <- fieller_from_fit(object, level = level)
    f$estimand <- names(b)[1]
    class(f) <- "drlate_fieller"
    return(f)
  }

  if (missing(parm)) parm <- names(b)
  if (identical(object$vcov_method, "bootstrap")) {
    probs <- c(a, 1 - a)
    out <- t(apply(object$boot$draws, 2, stats::quantile, probs = probs))
    rownames(out) <- names(b)
    out <- out[parm, , drop = FALSE]
  } else {
    se <- sqrt(diag(object$vcov3))
    q <- stats::qnorm(1 - a)
    out <- cbind(b - q * se, b + q * se)[parm, , drop = FALSE]
  }
  colnames(out) <- sprintf("%.1f %%", 100 * c(a, 1 - a))
  out
}

#' @export
print.drlate_fieller <- function(x, digits = 4, ...) {
  cat("Fieller ", format(100 * x$level), "% confidence set for ",
      x$estimand, ":\n  ", format_fieller(x, digits = digits), "\n",
      sep = "")
  invisible(x)
}
