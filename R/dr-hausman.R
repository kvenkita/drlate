#' Doubly robust Hausman test of unconfoundedness
#'
#' Tests whether the treatment is unconfounded given the covariates, using
#' the comparison proposed by Słoczyński, Uysal, and Wooldridge (2022,
#' Section 5), building on Donald, Hsu, and Lieli (2014). Under
#' **one-sided noncompliance** (nobody takes the treatment without the
#' instrument: \eqn{\Pr(D = 1 \mid Z = 0) = 0}), the LATT identified
#' through the instrument equals the ATT identified through
#' unconfoundedness of the treatment — so a significant difference between
#' the doubly robust LATT estimate (which uses the instrument) and the
#' doubly robust ATT estimate (which does not) is evidence against
#' unconfoundedness. Unlike the textbook OLS-vs-IV Hausman test, this
#' comparison is robust to treatment effect heterogeneity.
#'
#' The DR ATT estimator follows the paper's equation (33): a treatment
#' propensity score \eqn{\Pr(D = 1 \mid X)} is fitted by logit QMLE on the
#' treatment-equation covariates; the outcome model is fitted on the
#' untreated sample weighted by the odds \eqn{\hat p/(1-\hat p)}; and
#' \eqn{\hat\tau_{ATT}} is the treated-sample mean outcome minus the mean
#' imputed counterfactual. The standard error of the difference comes from
#' stacking the moment conditions of *both* estimators (and the difference)
#' into one M-estimation system, so the covariance between them is
#' accounted for analytically — the analytic option suggested in the paper.
#'
#' Note that the two halves adjust on their respective formulas: the LATT
#' half's propensity score uses the *instrument*-equation covariates,
#' while the ATT half's uses the *treatment*-equation covariates (both
#' share the outcome model). Supply the same covariate set to all three
#' formulas unless you intend them to differ.
#'
#' @inheritParams drlate
#' @param ivmodel Instrument propensity score model for the LATT half:
#'   `"logit"` (default) or `"ipt"`.
#'
#' @return An object of class `"htest"` with the z statistic, p-value, and
#'   the DR LATT, DR ATT, and difference estimates.
#'
#' @references
#' Słoczyński, T., S. D. Uysal, and J. M. Wooldridge (2022). "Doubly Robust
#' Estimation of Local Average Treatment Effects Using Inverse Probability
#' Weighted Regression Adjustment." \doi{10.48550/arXiv.2208.01300}
#'
#' Donald, S. G., Y.-C. Hsu, and R. P. Lieli (2014). "Testing the
#' Unconfoundedness Assumption via Inverse Probability Weighted Estimators
#' of (L)ATT." *Journal of Business & Economic Statistics* 32(3), 395-415.
#'
#' @examples
#' d <- drlate_sim
#' d$nvstat[d$rsncode == 0] <- 0L   # impose one-sided noncompliance
#' dr_hausman(lwage ~ age + educ, nvstat ~ age + educ,
#'            rsncode ~ age + educ, data = d)
#'
#' @export
dr_hausman <- function(outcome, treatment, instrument, data,
                       omodel = c("linear", "logit", "poisson"),
                       tmodel = c("logit", "linear", "poisson"),
                       ivmodel = c("logit", "ipt"),
                       weights = NULL, cluster = NULL,
                       pstolerance = 1e-5, subset = NULL) {
  omodel <- match.arg(omodel)
  tmodel <- match.arg(tmodel)
  ivmodel <- match.arg(ivmodel)
  fam <- c(linear = "gaussian", logit = "binomial", poisson = "poisson")
  dname <- paste(deparse(substitute(data)), collapse = "")

  if (!is.null(subset)) data <- data[subset, , drop = FALSE]
  if (is.character(weights) && length(weights) == 1L) {
    weights <- data[[weights]]
  }
  if (is.character(cluster) && length(cluster) == 1L) {
    cluster <- data[[cluster]]
  }

  ctx <- build_ctx(outcome, treatment, instrument, data,
                   omodel = fam[[omodel]], tmodel = fam[[tmodel]],
                   ivmodel = ivmodel, method = "ipwra", estimand = "latt",
                   normalized = TRUE, weights = weights, cluster = cluster,
                   pstolerance = pstolerance, osample = FALSE)

  if (ctx$dmeanz0 != 0) {
    stop("dr_hausman() requires one-sided noncompliance: nobody may take ",
         "the treatment without the instrument (mean of D in the Z = 0 ",
         "arm must be exactly 0; observed ",
         format(ctx$dmeanz0, digits = 4), ").", call. = FALSE)
  }

  # ---- DR LATT half: the existing machinery, blocks reused verbatim ----
  pt <- compute_point(ctx)
  latt_blocks <- pt$est$blocks
  latt_hat <- unname(pt$est$estimates["late"])

  # ---- DR ATT half (paper eq. 33): treatment PS, odds-weighted outcome --
  w <- ctx$w; d <- ctx$d; y <- ctx$y
  fitd <- suppressWarnings(
    stats::glm.fit(ctx$Xt, d, weights = w, family = stats::quasibinomial())
  )
  if (!fitd$converged) {
    stop("convergence not achieved for the treatment propensity score.",
         call. = FALSE)
  }
  gam <- stats::setNames(fitd$coefficients, colnames(ctx$Xt))
  dps <- stats::plogis(drop(ctx$Xt %*% gam))
  # Overlap for ATT (paper eq. 32): P(D=1|X) < 1
  if (any(dps > 1 - pstolerance)) {
    stop(sum(dps > 1 - pstolerance), " observation(s) have an estimated ",
         "treatment propensity score of (almost) 1; the ATT overlap ",
         "assumption fails.", call. = FALSE)
  }

  attw <- w * dps / (1 - dps)
  fya0 <- fit_wglm(y, ctx$Xo, ctx$omodel, attw, d == 0)
  mu0a <- predict_glm(fya0, ctx$Xo, ctx$omodel)
  ya1 <- wmean(y[d == 1], w[d == 1])
  att_hat <- ya1 - wmean(mu0a[d == 1], w[d == 1])
  diff_hat <- latt_hat - att_hat

  ones <- latt_ones(ctx)
  rwd <- rw_odds(ctx, eq = "dhat", X = ctx$Xt)
  att_blocks <- list(
    make_score_logit_block(ctx, "dhat", ctx$Xt, d, gam),
    make_glm_block(ctx, "ya0", ctx$omodel, ctx$Xo, y, 0, rwd, fya0,
                   arm_var = d),
    make_glm_block(ctx, "ya1", "gaussian", ones, y, 1, rw_one(ctx),
                   stats::setNames(ya1, "(Intercept)"), arm_var = d),
    make_contrast_block(ctx, "att",
                        pred_fun(ctx, "ya1", "gaussian", ones),
                        pred_fun(ctx, "ya0", ctx$omodel, ctx$Xo),
                        att_hat, arm = 1, arm_var = d)
  )
  diff_block <- make_contrast_block(ctx, "diff", term_param("late"),
                                    term_param("att"), diff_hat)

  sys <- assemble_system(c(latt_blocks, att_blocks, list(diff_block)))
  V <- drlate_vcov(sys, sys$theta0, ctx$w, ctx$cluster)
  se_diff <- sqrt(V[sys$layout$diff, sys$layout$diff])

  zstat <- diff_hat / se_diff
  pval <- 2 * stats::pnorm(-abs(zstat))

  structure(list(
    statistic = c(z = zstat),
    p.value = pval,
    estimate = c("DR LATT" = latt_hat, "DR ATT" = att_hat,
                 "difference" = diff_hat),
    stderr = se_diff,
    alternative = "two.sided",
    method = paste0("Doubly robust Hausman test of unconfoundedness\n",
                    "(Sloczynski-Uysal-Wooldridge 2022, ",
                    "one-sided noncompliance)"),
    data.name = dname
  ), class = "htest")
}
