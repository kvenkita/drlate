# Fieller confidence sets for the ratio late = num / denom.
#
# Inverts the test of H0: num - t * denom = 0 using the JOINT covariance
# of (num, denom) from the stacked sandwich (vcov_full). Unlike the
# delta-method (Wald) interval, the Fieller set has correct coverage even
# when the first stage is weak, at the price of possibly being unbounded.

#' Compute the Fieller confidence set.
#' @param num,denom Point estimates of the numerator and denominator.
#' @param V2 2x2 covariance matrix of (num, denom).
#' @param level Confidence level.
#' @return List with `lower`, `upper`, `type` in {"bounded",
#'   "complement", "whole-line"}, and the inputs.
#' @noRd
fieller_ci <- function(num, denom, V2, level = 0.95) {
  q <- stats::qchisq(level, df = 1)
  vnn <- V2[1, 1]; vnd <- V2[1, 2]; vdd <- V2[2, 2]

  # {t : (num - t*denom)^2 <= q * (vnn - 2 t vnd + t^2 vdd)}
  # <=> a t^2 + b t + c <= 0
  a <- denom^2 - q * vdd
  b <- -2 * (num * denom - q * vnd)
  cc <- num^2 - q * vnn
  disc <- b^2 - 4 * a * cc

  if (disc <= 0) {
    # No real roots: the inequality holds everywhere (a < 0) or nowhere;
    # "nowhere" cannot occur at the point estimate (t = num/denom always
    # satisfies it), so the set is the whole line.
    return(list(lower = -Inf, upper = Inf, type = "whole-line",
                level = level))
  }
  r <- sort((-b + c(-1, 1) * sqrt(disc)) / (2 * a))
  if (a > 0) {
    list(lower = r[1], upper = r[2], type = "bounded", level = level)
  } else {
    # a < 0: the set is the complement (-Inf, r1] U [r2, Inf)
    list(lower = r[1], upper = r[2], type = "complement", level = level)
  }
}

#' Format a Fieller set for printing
#' @noRd
format_fieller <- function(f, digits = 4) {
  fmt <- function(v) format(v, digits = digits)
  switch(f$type,
    bounded    = paste0("[", fmt(f$lower), ", ", fmt(f$upper), "]"),
    complement = paste0("(-Inf, ", fmt(f$lower), "] U [",
                        fmt(f$upper), ", Inf)"),
    `whole-line` = "(-Inf, Inf) - the first stage is uninformative"
  )
}

#' Fieller set from a fitted drlate object
#' @noRd
fieller_from_fit <- function(object, level = 0.95) {
  idx <- c(object$layout$num, object$layout$denom)
  V2 <- object$vcov_full[idx, idx]
  fieller_ci(unname(object$coefficients[2]),
             unname(object$coefficients[3]), V2, level)
}
