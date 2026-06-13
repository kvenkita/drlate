# Imai-Ratkovic (2014) overidentification test for covariate balance after a
# drlate fit (the latebalance overid postestimation feature of lateffects).
#
# The instrument propensity-score coefficients solve the just-identified MLE
# score equations. The covariate-balancing moments (the CBPS conditions) are
# then overidentifying restrictions: if the propensity-score model balances the
# covariates, they hold at the fitted coefficients. The test statistic is the
# conditional-moment (Hansen-J-type) quadratic form in the balancing moments,
# evaluated at the fitted coefficients, with a variance that accounts for the
# first-stage estimation through the influence function.

#' Imai-Ratkovic covariate-balance test
#'
#' Tests whether the estimated instrument propensity score balances the
#' covariates, using the overidentification test of Imai and Ratkovic (2014).
#' The propensity-score MLE score equations identify the coefficients; the
#' covariate-balancing (CBPS) moments are the overidentifying restrictions. A
#' large statistic is evidence that the propensity-score model does not balance
#' the covariates --- a misspecification diagnostic. This is the
#' \proglang{Stata} \code{latebalance overid} postestimation feature.
#'
#' @param object A fitted [drlate()] object (with `keep_data = TRUE`) using a
#'   logistic or probit instrument propensity score.
#' @return An object of class `drlate_balance_test`: a list with `statistic`
#'   (Hansen's J), `df`, `p.value`, `ivmodel`, and `n`, with a `print` method.
#' @references Imai, K. and Ratkovic, M. (2014). Covariate Balancing Propensity
#'   Score. \emph{Journal of the Royal Statistical Society B} 76(1), 243--263.
#' @seealso [balance()] for the standardized-mean-difference diagnostics.
#' @export
balance_test <- function(object) {
  ctx <- need_ctx(object)
  ps <- need_ps(object)
  if (ctx$ivmodel == "cbps") {
    stop("balance_test() is not available after ivmodel = \"cbps\": the ",
         "covariate-balancing propensity score balances the covariates by ",
         "construction, so the overidentifying restrictions hold mechanically ",
         "(cf. Stata's latebalance, unavailable after lateffects balancing).",
         call. = FALSE)
  }
  if (!ctx$ivmodel %in% c("logit", "probit")) {
    stop("balance_test() requires ivmodel = \"logit\" or \"probit\".",
         call. = FALSE)
  }

  X <- ctx$Xz; z <- ctx$z; w <- ctx$w; n <- ctx$n
  beta <- ps$bips
  probit <- ctx$ivmodel == "probit"
  linkinv <- if (probit) stats::pnorm else stats::plogis

  # Propensity-score MLE score (matches make_ps_logit_block / probit block)
  score_mom <- function(b) {
    xb <- drop(X %*% b)
    p <- linkinv(xb)
    if (probit) ((z - p) / (p * (1 - p)) * stats::dnorm(xb)) * X
    else (z - p) * X
  }
  # Covariate-balancing moment (matches make_ps_cbps_block)
  bal_mom <- function(b) {
    p <- linkinv(drop(X %*% b))
    (z / p - (1 - z) / (1 - p)) * X
  }
  sbar <- function(b) colMeans(w * score_mom(b))
  bbar <- function(b) colMeans(w * bal_mom(b))

  A_s <- numDeriv::jacobian(sbar, beta)
  G_b <- numDeriv::jacobian(bbar, beta)
  Ainv <- tryCatch(solve(A_s), error = function(e) {
    stop("singular score Jacobian in the balance test.", call. = FALSE)
  })

  Sw <- w * score_mom(beta)
  Bw <- w * bal_mom(beta)
  # Influence function of the balancing moments, net of first-stage estimation:
  #   psi_i = b_i - G_b A_s^{-1} s_i
  Psi <- Bw - Sw %*% t(Ainv) %*% t(G_b)

  if (is.null(ctx$cluster)) {
    Omega <- crossprod(Psi) / n
  } else {
    Omega <- crossprod(rowsum(Psi, group = ctx$cluster)) / n
  }
  Oinv <- tryCatch(solve(Omega), error = function(e) {
    stop("singular balance-moment variance; the test could not be computed.",
         call. = FALSE)
  })

  bm <- colMeans(Bw)
  stat <- n * drop(t(bm) %*% Oinv %*% bm)
  df <- ncol(X)
  structure(
    list(statistic = stat, df = df,
         p.value = stats::pchisq(stat, df, lower.tail = FALSE),
         ivmodel = ctx$ivmodel, n = n),
    class = "drlate_balance_test"
  )
}

#' @export
print.drlate_balance_test <- function(x, ...) {
  cat("Imai-Ratkovic covariate-balance test (overidentification)\n\n")
  cat(sprintf("  Hansen J = %.4f   df = %d   p-value = %.4f\n",
              x$statistic, x$df, x$p.value))
  cat(sprintf("  Instrument propensity score: %s (n = %d)\n",
              x$ivmodel, x$n))
  cat("\n  H0: the propensity-score model balances the covariates.\n")
  invisible(x)
}
