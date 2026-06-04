# Moment-block architecture.
#
# Every estimator in drlate is a just-identified stacked M-estimation system:
# the per-observation moment conditions of the instrument propensity score,
# the outcome and treatment regressions, and the scalar aggregates
# (num, denom, late, plus normalization constants where applicable).
# Point estimates come from sequential fits (fit-nuisance.R); the stacked
# system exists solely to compute the joint sandwich variance, exactly like
# Stata's `gmm ..., onestep from(...) iterate(0)`.
#
# A *block* owns a slice of the global parameter vector theta and contributes
# a set of per-observation moment columns:
#   list(eq, parnames, k, start, g(theta, layout) -> n x k matrix)
# `layout` maps equation names to index ranges in theta, so blocks can
# reference other blocks' parameters (e.g. the PS linear index `zhat`).

#' @noRd
new_block <- function(eq, parnames, start, g) {
  stopifnot(length(parnames) == length(start))
  list(eq = eq, parnames = parnames, k = length(start),
       start = unname(start), g = g)
}

#' Linear combination of an equation's parameters (Stata gmm `{eq:}`)
#' @noRd
lc <- function(theta, layout, eq, X) drop(X %*% theta[layout[[eq]]])

# ---------------------------------------------------------------------------
# Instrument propensity score blocks
# ---------------------------------------------------------------------------

#' Logit MLE score: (z - p) * x
#' @noRd
make_ps_logit_block <- function(ctx, coefs, eq = "zhat") {
  X <- ctx$Xz; z <- ctx$z
  new_block(eq, paste0(eq, ":", colnames(X)), coefs,
    function(theta, layout) (z - stats::plogis(lc(theta, layout, eq, X))) * X)
}

#' CBPS balancing moment: (z/p - (1-z)/(1-p)) * x
#' @noRd
make_ps_cbps_block <- function(ctx, coefs, eq = "zhat") {
  X <- ctx$Xz; z <- ctx$z
  new_block(eq, paste0(eq, ":", colnames(X)), coefs,
    function(theta, layout) {
      p <- stats::plogis(lc(theta, layout, eq, X))
      (z / p - (1 - z) / (1 - p)) * X
    })
}

#' IPT tilting moment for the Z=1 arm: (z * (1 + e^xb) / e^xb - 1) * x
#' @noRd
make_ps_ipt1_block <- function(ctx, coefs, eq = "zhat1") {
  X <- ctx$Xz; z <- ctx$z
  new_block(eq, paste0(eq, ":", colnames(X)), coefs,
    function(theta, layout) {
      xb <- lc(theta, layout, eq, X)
      (z * (1 + exp(xb)) / exp(xb) - 1) * X
    })
}

#' IPT tilting moment for the Z=0 arm: ((1-z) * (1 + e^xb) - 1) * x
#' @noRd
make_ps_ipt0_block <- function(ctx, coefs, eq = "zhat0") {
  X <- ctx$Xz; z <- ctx$z
  new_block(eq, paste0(eq, ":", colnames(X)), coefs,
    function(theta, layout) {
      xb <- lc(theta, layout, eq, X)
      ((1 - z) * (1 + exp(xb)) - 1) * X
    })
}

# ---------------------------------------------------------------------------
# Reweight factors (the IPW weights inlined into the moment conditions,
# expressed through the PS linear index so first-stage uncertainty
# propagates through the Jacobian)
# ---------------------------------------------------------------------------

#' No reweighting (RA, AIPW regression scores)
#' @noRd
rw_one <- function(ctx) function(theta, layout) 1

#' 1/p = 1 + exp(-zhat)   (Z=1 arm of IPWRA / IPW)
#' @noRd
rw_invp <- function(ctx, eq = "zhat") {
  X <- ctx$Xz
  function(theta, layout) 1 + exp(-lc(theta, layout, eq, X))
}

#' 1/(1-p) = 1 + exp(zhat)   (Z=0 arm of IPWRA / IPW)
#' @noRd
rw_inv1mp <- function(ctx, eq = "zhat") {
  X <- ctx$Xz
  function(theta, layout) 1 + exp(lc(theta, layout, eq, X))
}

#' Odds p/(1-p) = exp(zhat)   (Z=0 arm of LATT estimators)
#' @noRd
rw_odds <- function(ctx, eq = "zhat") {
  X <- ctx$Xz
  function(theta, layout) exp(lc(theta, layout, eq, X))
}

# ---------------------------------------------------------------------------
# Regression blocks (outcome and treatment, all families, all methods)
# ---------------------------------------------------------------------------

#' Weighted GLM score restricted to one instrument arm:
#'   arm_indicator * reweight(theta) * (yvar - mu(X beta)) * X
#' One constructor covers eqy0/eqy1/eqd0/eqd1 across IPWRA, RA, AIPW and the
#' regression pieces of LATT; the (family, arm, reweight) triple is the
#' entire variation between Stata's ~40 hand-written gmm equation strings.
#' @noRd
make_glm_block <- function(ctx, eq, family, X, yvar, arm, reweight, coefs) {
  linkinv <- fam_linkinv(family)
  a <- if (is.null(arm)) rep(1, ctx$n) else if (arm == 1) ctx$z else 1 - ctx$z
  new_block(eq, paste0(eq, ":", colnames(X)), coefs,
    function(theta, layout) {
      mu <- linkinv(lc(theta, layout, eq, X))
      (a * reweight(theta, layout) * (yvar - mu)) * X
    })
}

# ---------------------------------------------------------------------------
# Prediction closures used by the scalar aggregate blocks
# ---------------------------------------------------------------------------

#' Per-observation prediction of a fitted regression equation, or a constant
#' when the corresponding arm is degenerate (one-sided noncompliance).
#' @noRd
pred_fun <- function(ctx, eq, family, X, degenerate_value = NULL) {
  if (!is.null(degenerate_value)) {
    force(degenerate_value)
    return(function(theta, layout) rep(degenerate_value, ctx$n))
  }
  linkinv <- fam_linkinv(family)
  function(theta, layout) linkinv(lc(theta, layout, eq, X))
}

# ---------------------------------------------------------------------------
# Scalar aggregate blocks
# ---------------------------------------------------------------------------

#' Scalar moment: param - (term1(theta) - term0(theta)), optionally restricted
#' to an instrument arm (LATT aggregates are means over the Z=1 subsample).
#' @noRd
make_contrast_block <- function(ctx, eq, term1, term0, start, arm = NULL) {
  a <- if (is.null(arm)) rep(1, ctx$n) else if (arm == 1) ctx$z else 1 - ctx$z
  new_block(eq, eq, start,
    function(theta, layout) {
      cbind(a * (theta[layout[[eq]]] - (term1(theta, layout) -
                                          term0(theta, layout))))
    })
}

#' Scalar moment: param - term(theta)
#' @noRd
make_scalar_block <- function(ctx, eq, term, start, arm = NULL) {
  a <- if (is.null(arm)) rep(1, ctx$n) else if (arm == 1) ctx$z else 1 - ctx$z
  new_block(eq, eq, start,
    function(theta, layout) {
      cbind(a * (theta[layout[[eq]]] - term(theta, layout)))
    })
}

#' The ratio: late - num/denom
#' @noRd
make_late_block <- function(ctx, start) {
  n <- ctx$n
  new_block("late", "late", start,
    function(theta, layout) {
      cbind(rep(theta[layout$late] -
                  theta[layout$num] / theta[layout$denom], n))
    })
}

# ---------------------------------------------------------------------------
# Assembler
# ---------------------------------------------------------------------------

#' Concatenate blocks into a square moment system.
#' Returns theta0 (named starting values = the reported point estimates),
#' g(theta) -> n x p matrix of per-observation moments, and the layout.
#' @noRd
assemble_system <- function(blocks) {
  layout <- list(); pnames <- character(0); theta0 <- numeric(0); pos <- 0L
  for (b in blocks) {
    if (!is.null(layout[[b$eq]])) {
      stop("internal error: duplicate equation name ", b$eq)
    }
    layout[[b$eq]] <- pos + seq_len(b$k)
    pos <- pos + b$k
    pnames <- c(pnames, b$parnames)
    theta0 <- c(theta0, b$start)
  }
  gfun <- function(theta) {
    do.call(cbind, lapply(blocks, function(b) {
      m <- b$g(theta, layout)
      if (is.null(dim(m))) m <- cbind(m)
      m
    }))
  }
  list(theta0 = stats::setNames(theta0, pnames), g = gfun,
       layout = layout, p = pos)
}
