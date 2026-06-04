# LATE estimators: block dispatch per (method, normalization, ivmodel, case).
# Mirrors drlate_estimate_late.ado; block order matches Stata's gmm equation
# order for each case so parameter layouts are directly comparable.

#' @noRd
estimate_late <- function(ctx, ps) {
  switch(ctx$method,
    ipwra = late_ipwra(ctx, ps),
    stop("method = \"", ctx$method, "\" is not implemented yet.",
         call. = FALSE)
  )
}

#' LATE via inverse-probability-weighted regression adjustment.
#' Stata: drlate_estimate_late.ado, METHOD: IPWRA (lines 53-323).
#' @noRd
late_ipwra <- function(ctx, ps) {
  if (ctx$ivmodel == "ipt") {
    stop("ivmodel = \"ipt\" is not implemented yet.", call. = FALSE)
  }
  w <- ctx$w; z <- ctx$z

  # Fitting weights: w/p on the Z=1 arm, w/(1-p) on the Z=0 arm
  fw1 <- w / ps$ps
  fw0 <- w / (1 - ps$ps)

  # Outcome regressions (never degenerate)
  fy1 <- fit_wglm(ctx$y, ctx$Xo, ctx$omodel, fw1, z == 1)
  fy0 <- fit_wglm(ctx$y, ctx$Xo, ctx$omodel, fw0, z == 0)
  mu_y1 <- predict_glm(fy1, ctx$Xo, ctx$omodel)
  mu_y0 <- predict_glm(fy0, ctx$Xo, ctx$omodel)

  # Treatment regressions (skipped for a degenerate arm)
  a1 <- fit_arm(ctx, ctx$d, ctx$Xt, ctx$tmodel, 1, fw1, ctx$dmeanz1)
  a0 <- fit_arm(ctx, ctx$d, ctx$Xt, ctx$tmodel, 0, fw0, ctx$dmeanz0)

  num   <- wmean(mu_y1, w) - wmean(mu_y0, w)
  denom <- wmean(a1$mu, w) - wmean(a0$mu, w)
  late  <- num / denom

  # PS block first (eqips), then eqy0, eqy1, num, eqd0, eqd1, denom, late
  ps_block <- switch(ctx$ivmodel,
    logit = make_ps_logit_block(ctx, ps$bips),
    cbps  = make_ps_cbps_block(ctx, ps$bips)
  )

  blocks <- list(
    ps_block,
    make_glm_block(ctx, "y0", ctx$omodel, ctx$Xo, ctx$y, 0,
                   rw_inv1mp(ctx), fy0),
    make_glm_block(ctx, "y1", ctx$omodel, ctx$Xo, ctx$y, 1,
                   rw_invp(ctx), fy1),
    make_contrast_block(ctx, "num",
                        pred_fun(ctx, "y1", ctx$omodel, ctx$Xo),
                        pred_fun(ctx, "y0", ctx$omodel, ctx$Xo), num)
  )
  if (is.null(a0$degenerate_value)) {
    blocks <- c(blocks, list(
      make_glm_block(ctx, "d0", ctx$tmodel, ctx$Xt, ctx$d, 0,
                     rw_inv1mp(ctx), a0$coefs)))
  }
  if (is.null(a1$degenerate_value)) {
    blocks <- c(blocks, list(
      make_glm_block(ctx, "d1", ctx$tmodel, ctx$Xt, ctx$d, 1,
                     rw_invp(ctx), a1$coefs)))
  }
  blocks <- c(blocks, list(
    make_contrast_block(ctx, "denom",
                        pred_fun(ctx, "d1", ctx$tmodel, ctx$Xt,
                                 a1$degenerate_value),
                        pred_fun(ctx, "d0", ctx$tmodel, ctx$Xt,
                                 a0$degenerate_value), denom),
    make_late_block(ctx, late)
  ))

  list(blocks = blocks,
       estimates = c(late = late, num = num, denom = denom))
}
