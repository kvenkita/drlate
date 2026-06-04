# S3 methods for class "drlate".

#' @export
print.drlate <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  title <- if (x$estimand == "late") "Local average treatment effect"
           else "Local average treatment effect on the treated"
  methodd <- toupper(x$method)
  normd <- if (x$ivmodel == "ipt" || x$statnorm == "nrm") "normalized"
           else "unnormalized"
  est_line <- if (x$method %in% c("ra", "ipwra")) methodd
              else paste0(methodd, " (", normd, ")")
  omodeld <- x$omodel
  tmodeld <- x$tmodel
  if (x$method == "ipw") {
    omodeld <- "weighted mean"
    tmodeld <- "weighted mean"
  }
  zmodeld <- switch(x$ivmodel,
    logit = "logit (MLE)", cbps = "logit (CBPS)", ipt = "logit (IPT)")

  cat("\n", title, "\n", sep = "")
  cat("Number of obs    : ", format(x$N, big.mark = ","), "\n", sep = "")
  if (!is.null(x$N_clust)) {
    cat("Number of clusters: ", x$N_clust, "\n", sep = "")
  }
  cat("Estimator        : ", est_line, "\n", sep = "")
  cat("Outcome model    : ", omodeld, "\n", sep = "")
  cat("Treatment model  : ", tmodeld, "\n", sep = "")
  cat("Instrument model : ", zmodeld, "\n\n", sep = "")

  print(coeftable(x), digits = digits)

  fz <- firststage_z(x)
  cat("\nFirst stage (Z on D): z = ", format(fz, digits = 4), sep = "")
  if (is.finite(fz) && abs(fz) < 2) {
    cat("  [weak: the LATE ratio may be unstable;",
        "see confint(., method = \"fieller\")]")
  }
  cat("\n")
  invisible(x)
}

#' @noRd
coeftable <- function(x) {
  b <- x$coefficients
  se <- sqrt(diag(x$vcov3))
  zstat <- b / se
  p <- 2 * stats::pnorm(-abs(zstat))
  ci <- cbind(b - stats::qnorm(0.975) * se, b + stats::qnorm(0.975) * se)
  out <- cbind(b, se, zstat, p, ci)
  dimnames(out) <- list(names(b),
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)",
      "[95% conf.", "interval]"))
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

#' @export
confint.drlate <- function(object, parm, level = 0.95, ...) {
  b <- object$coefficients
  se <- sqrt(diag(object$vcov3))
  if (missing(parm)) parm <- names(b)
  a <- (1 - level) / 2
  q <- stats::qnorm(1 - a)
  out <- cbind(b - q * se, b + q * se)[parm, , drop = FALSE]
  colnames(out) <- sprintf("%.1f %%", 100 * c(a, 1 - a))
  out
}
