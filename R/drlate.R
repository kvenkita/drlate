#' Doubly robust estimation of the LATE and LATT
#'
#' Estimates the local average treatment effect (LATE) or the local average
#' treatment effect on the treated (LATT) with a binary instrument, following
#' Słoczyński, Uysal, and Wooldridge (2022). A faithful R port of the Stata
#' package `drlate` (SSC S459708): point estimates come from sequential
#' weighted regressions, and standard errors are computed jointly for the
#' instrument propensity score, the outcome regression, the treatment
#' regression, and the causal estimand by stacking all moment conditions
#' into a single M-estimation system.
#'
#' @param outcome A formula `y ~ covariates` for the outcome model. Use
#'   `y ~ 1` for no covariates (required when `method = "ipw"`).
#' @param treatment A formula `d ~ covariates` for the treatment model.
#' @param instrument A formula `z ~ covariates` for the instrument
#'   propensity score model; `z` must be binary 0/1. Use `z ~ 1` when
#'   `method = "ra"`.
#' @param data A data frame containing all variables.
#' @param omodel Outcome model family: `"linear"` (default), `"logit"`
#'   (outcome must be 0/1), or `"poisson"` (outcome must be non-negative).
#' @param tmodel Treatment model family: `"logit"` (default; treatment must
#'   be 0/1), `"linear"`, or `"poisson"`.
#' @param ivmodel Instrument propensity score model: `"logit"` (maximum
#'   likelihood; default), `"cbps"` (covariate balancing, Imai and Ratkovic
#'   2014; not available with `estimand = "latt"`), or `"ipt"` (inverse
#'   probability tilting, Graham, Pinto, and Egel 2012).
#' @param method Estimator: `"ipwra"` (inverse-probability-weighted
#'   regression adjustment; default), `"ipw"`, `"aipw"`, or `"ra"`.
#' @param estimand `"late"` (default) or `"latt"`.
#' @param normalized Logical; use normalized moment conditions (default
#'   `TRUE`). Only relevant for `method = "ipw"` and `method = "aipw"`.
#' @param weights Optional sampling weights (a numeric vector, or a column
#'   name in `data` given as a string).
#' @param cluster Optional cluster identifier for clustered standard errors
#'   (a vector, or a column name in `data` given as a string).
#' @param pstolerance Overlap tolerance: estimation stops with an error if
#'   any estimated instrument propensity score is below `pstolerance` or
#'   above `1 - pstolerance`. Default `1e-5`.
#' @param osample Logical; if `TRUE`, overlap violations do not stop
#'   estimation with an error. Instead `drlate()` returns (invisibly) a
#'   logical vector marking the violating observations.
#' @param subset Optional logical or integer vector selecting rows of `data`.
#'
#' @return An object of class `"drlate"`, a list with components including
#'   `coefficients` (the causal estimate, the numerator effect of Z on Y,
#'   and the denominator effect of Z on D), `vcov3` (their variance matrix,
#'   diagonal by construction, as in the Stata package), `vcov_full` (the
#'   joint variance matrix of all stacked parameters), `theta` (all stacked
#'   parameter estimates), `N`, `dmeanz1`, `dmeanz0`, and the call.
#'
#' @references
#' Słoczyński, T., S. D. Uysal, and J. M. Wooldridge (2022). "Doubly Robust
#' Estimation of Local Average Treatment Effects Using Inverse Probability
#' Weighted Regression Adjustment." \doi{10.48550/arXiv.2208.01300}
#'
#' @examples
#' data(drlate_sim)
#' fit <- drlate(lwage ~ age + educ, nvstat ~ age + educ,
#'               rsncode ~ age + educ, data = drlate_sim)
#' summary(fit)
#'
#' @export
drlate <- function(outcome, treatment, instrument, data,
                   omodel = c("linear", "logit", "poisson"),
                   tmodel = c("logit", "linear", "poisson"),
                   ivmodel = c("logit", "cbps", "ipt"),
                   method = c("ipwra", "ipw", "aipw", "ra"),
                   estimand = c("late", "latt"),
                   normalized = TRUE,
                   weights = NULL, cluster = NULL,
                   pstolerance = 1e-5, osample = FALSE,
                   subset = NULL) {
  cl <- match.call()
  omodel <- match.arg(omodel)
  tmodel <- match.arg(tmodel)
  ivmodel <- match.arg(ivmodel)
  method <- match.arg(method)
  estimand <- match.arg(estimand)

  # Map user-facing names to internal family strings
  fam <- c(linear = "gaussian", logit = "binomial", poisson = "poisson")

  if (!is.null(subset)) data <- data[subset, , drop = FALSE]
  if (is.character(weights) && length(weights) == 1L) {
    weights <- data[[weights]]
  }
  if (is.character(cluster) && length(cluster) == 1L) {
    cluster <- data[[cluster]]
  }

  ctx <- build_ctx(outcome, treatment, instrument, data,
                   omodel = fam[[omodel]], tmodel = fam[[tmodel]],
                   ivmodel = ivmodel, method = method, estimand = estimand,
                   normalized = normalized, weights = weights,
                   cluster = cluster, pstolerance = pstolerance,
                   osample = osample)

  # Instrument propensity score (not used by RA) + overlap check
  ps <- NULL
  if (method != "ra") {
    ps <- fit_ps(ctx)
    viol <- check_overlap(ps$ps, pstolerance, osample)
    if (osample && any(viol)) {
      message(sum(viol), " observation(s) violate the overlap assumption; ",
              "returning the violation indicator.")
      return(invisible(viol))
    }
  }

  # Normalize check (drlate_estimate.ado section 7): when the IPW weights
  # already average to one within rounding (IPT weights do by construction),
  # the normalized and unnormalized moments coincide and Stata switches to
  # the unnormalized system.
  if (!is.null(ps) && ctx$statnorm == "nrm" && method != "ra") {
    wt1m <- round(wmean(ps$wt1, ctx$w), 6)
    wt0m <- round(wmean(ps$wt0, ctx$w), 6)
    if (wt1m == 1 && wt0m == 1) ctx$statnorm <- "unnrm"
  }
  # IPT weights are ex-ante normalized (drlate_estimate_late.ado 658-661)
  if (method == "aipw" && ivmodel == "ipt" && ctx$statnorm == "nrm") {
    message("IPT weights are ex-ante normalized; switching to unnormalized ",
            "moments.")
    ctx$statnorm <- "unnrm"
  }

  est <- if (estimand == "late") estimate_late(ctx, ps)
         else estimate_latt(ctx, ps)

  sys <- assemble_system(est$blocks)
  V <- drlate_vcov(sys, sys$theta0, ctx$w, ctx$cluster)

  # Extract the reported triple by name (Stata: late/num/denom _cons)
  labels <- if (estimand == "late") {
    c("LATE: D on Y", "ATE: Z on Y", "ATE: Z on D")
  } else {
    c("LATT: D on Y", "ATT: Z on Y", "ATT: Z on D")
  }
  idx <- c(sys$layout$late, sys$layout$num, sys$layout$denom)
  b3 <- stats::setNames(sys$theta0[idx], labels)
  # Stata reports a diagonal 3x3 V (drlate_estimate.ado section 10)
  V3 <- diag(diag(V[idx, idx, drop = FALSE]), nrow = 3L)
  dimnames(V3) <- list(labels, labels)

  structure(list(
    coefficients = b3,
    vcov3 = V3,
    vcov_full = V,
    theta = sys$theta0,
    N = ctx$n,
    N_clust = if (is.null(ctx$cluster)) NULL
              else length(unique(ctx$cluster)),
    dmeanz1 = ctx$dmeanz1,
    dmeanz0 = ctx$dmeanz0,
    estimand = estimand, method = method,
    omodel = omodel, tmodel = tmodel, ivmodel = ivmodel,
    statnorm = ctx$statnorm,
    case = ctx$case,
    call = cl
  ), class = "drlate")
}
