# LATE estimators: block dispatch per (method, normalization, ivmodel, case).
# Mirrors drlate_estimate_late.ado; block order matches Stata's gmm equation
# order for each case so parameter layouts are directly comparable.

#' @noRd
estimate_late <- function(ctx, ps) {
  if (ctx$case == "bothdeg" && ctx$method != "ipwra") {
    stop("treatment is degenerate in both instrument arms; the effect of ",
         "Z on D is constant and no estimation is needed.", call. = FALSE)
  }
  switch(ctx$method,
    ipwra  = late_ipwra(ctx, ps),
    ra     = late_ra(ctx),
    ipw    = late_ipw(ctx, ps),
    aipw   = late_aipw(ctx, ps),
    kappa   = late_kappa(ctx, ps),
    kappa0  = late_kappa0(ctx, ps),
    kappa10 = late_kappa10(ctx, ps),
    stop("method = \"", ctx$method, "\" is not implemented yet.",
         call. = FALSE)
  )
}

#' Shared: PS moment block(s) and the arm reweight functions.
#' For logit/cbps a single eqips block and reweights through `zhat`;
#' for ipt two tilt blocks (eqips1, eqips0) and reweights through
#' `zhat1`/`zhat0` (drlate_estimate_late.ado lines 25-38).
#' @noRd
late_ps_setup <- function(ctx, ps) {
  if (ctx$ivmodel == "ipt") {
    list(blocks = list(make_ps_ipt1_block(ctx, ps$bips1),
                       make_ps_ipt0_block(ctx, ps$bips0)),
         rw1 = rw_invp(ctx, "zhat1"),
         rw0 = rw_inv1mp(ctx, "zhat0"))
  } else {
    blk <- switch(ctx$ivmodel,
      logit = make_ps_logit_block(ctx, ps$bips),
      cbps  = make_ps_cbps_block(ctx, ps$bips))
    list(blocks = list(blk), rw1 = rw_invp(ctx), rw0 = rw_inv1mp(ctx))
  }
}

#' LATE via inverse-probability-weighted regression adjustment.
#' Stata: drlate_estimate_late.ado, METHOD: IPWRA (lines 53-323).
#' @noRd
late_ipwra <- function(ctx, ps) {
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

  # PS block(s) first, then eqy0, eqy1, num, eqd0, eqd1, denom, late
  setup <- late_ps_setup(ctx, ps)
  blocks <- c(setup$blocks, list(
    make_glm_block(ctx, "y0", ctx$omodel, ctx$Xo, ctx$y, 0,
                   setup$rw0, fy0),
    make_glm_block(ctx, "y1", ctx$omodel, ctx$Xo, ctx$y, 1,
                   setup$rw1, fy1),
    make_contrast_block(ctx, "num",
                        pred_fun(ctx, "y1", ctx$omodel, ctx$Xo),
                        pred_fun(ctx, "y0", ctx$omodel, ctx$Xo), num)
  ))
  if (is.null(a0$degenerate_value)) {
    blocks <- c(blocks, list(
      make_glm_block(ctx, "d0", ctx$tmodel, ctx$Xt, ctx$d, 0,
                     setup$rw0, a0$coefs)))
  }
  if (is.null(a1$degenerate_value)) {
    blocks <- c(blocks, list(
      make_glm_block(ctx, "d1", ctx$tmodel, ctx$Xt, ctx$d, 1,
                     setup$rw1, a1$coefs)))
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

#' LATE via regression adjustment (no instrument propensity score).
#' Stata: drlate_estimate_late.ado, METHOD: RA (lines 1255-1426).
#' RA introduces separate scalar parameters num0/num1 (and denom0/denom1
#' for non-degenerate arms) for the averaged predictions.
#' @noRd
late_ra <- function(ctx) {
  w <- ctx$w; z <- ctx$z

  fy1 <- fit_wglm(ctx$y, ctx$Xo, ctx$omodel, w, z == 1)
  fy0 <- fit_wglm(ctx$y, ctx$Xo, ctx$omodel, w, z == 0)
  mu_y1 <- predict_glm(fy1, ctx$Xo, ctx$omodel)
  mu_y0 <- predict_glm(fy0, ctx$Xo, ctx$omodel)
  a1 <- fit_arm(ctx, ctx$d, ctx$Xt, ctx$tmodel, 1, w, ctx$dmeanz1)
  a0 <- fit_arm(ctx, ctx$d, ctx$Xt, ctx$tmodel, 0, w, ctx$dmeanz0)

  num1s <- wmean(mu_y1, w); num0s <- wmean(mu_y0, w)
  den1s <- wmean(a1$mu, w); den0s <- wmean(a0$mu, w)
  num   <- num1s - num0s
  denom <- den1s - den0s
  late  <- num / denom

  d0deg <- !is.null(a0$degenerate_value)
  d1deg <- !is.null(a1$degenerate_value)

  # Order (interior): eqy0 eqy1 num0 num1 num eqd0 eqd1 denom0 denom1
  #                   denom late
  blocks <- list(
    make_glm_block(ctx, "y0", ctx$omodel, ctx$Xo, ctx$y, 0, rw_one(ctx), fy0),
    make_glm_block(ctx, "y1", ctx$omodel, ctx$Xo, ctx$y, 1, rw_one(ctx), fy1),
    make_scalar_block(ctx, "num0",
                      pred_fun(ctx, "y0", ctx$omodel, ctx$Xo), num0s),
    make_scalar_block(ctx, "num1",
                      pred_fun(ctx, "y1", ctx$omodel, ctx$Xo), num1s),
    make_contrast_block(ctx, "num", term_param("num1"), term_param("num0"),
                        num)
  )
  if (!d0deg) {
    blocks <- c(blocks, list(
      make_glm_block(ctx, "d0", ctx$tmodel, ctx$Xt, ctx$d, 0, rw_one(ctx),
                     a0$coefs)))
  }
  if (!d1deg) {
    blocks <- c(blocks, list(
      make_glm_block(ctx, "d1", ctx$tmodel, ctx$Xt, ctx$d, 1, rw_one(ctx),
                     a1$coefs)))
  }
  if (!d0deg) {
    blocks <- c(blocks, list(
      make_scalar_block(ctx, "denom0",
                        pred_fun(ctx, "d0", ctx$tmodel, ctx$Xt), den0s)))
  }
  if (!d1deg) {
    blocks <- c(blocks, list(
      make_scalar_block(ctx, "denom1",
                        pred_fun(ctx, "d1", ctx$tmodel, ctx$Xt), den1s)))
  }
  t1 <- if (d1deg) term_const(a1$degenerate_value) else term_param("denom1")
  t0 <- if (d0deg) term_const(a0$degenerate_value) else term_param("denom0")
  blocks <- c(blocks, list(
    make_contrast_block(ctx, "denom", t1, t0, denom),
    make_late_block(ctx, late)
  ))

  list(blocks = blocks,
       estimates = c(late = late, num = num, denom = denom))
}

#' LATE via inverse probability weighting.
#' Stata: drlate_estimate_late.ado, METHOD: IPW (lines 329-651).
#' Normalized: Hajek means via intercept-only weighted regressions (always
#' the linear `regress` form, regardless of tmodel). Unnormalized: raw IPW
#' moments; only y0/d0 carry (decorative) parameters, following Stata.
#' @noRd
late_ipw <- function(ctx, ps) {
  w <- ctx$w; z <- ctx$z
  ones <- matrix(1, ctx$n, 1L, dimnames = list(NULL, "(Intercept)"))
  setup <- late_ps_setup(ctx, ps)
  rw1 <- setup$rw1   # z-arm reweight 1/p, via the PS linear index
  rw0 <- setup$rw0   # 1/(1-p)
  d0deg <- ctx$dmeanz0 %in% c(0, 1)
  d1deg <- ctx$dmeanz1 %in% c(0, 1)

  # With IPT the weights are self-normalized, so Stata's "unnormalized"
  # IPT branch uses the same Hajek/parameter moment structure as the
  # normalized branch (drlate_estimate_late.ado lines 593-650).
  if (ctx$statnorm == "nrm" || ctx$ivmodel == "ipt") {
    fw1 <- w / ps$ps
    fw0 <- w / (1 - ps$ps)
    # Intercept-only `regress` fits = Hajek-weighted means
    fy1 <- fit_wglm(ctx$y, ones, "gaussian", fw1, z == 1)
    fy0 <- fit_wglm(ctx$y, ones, "gaussian", fw0, z == 0)
    fd1 <- if (d1deg) NULL else fit_wglm(ctx$d, ones, "gaussian", fw1, z == 1)
    fd0 <- if (d0deg) NULL else fit_wglm(ctx$d, ones, "gaussian", fw0, z == 0)

    num <- unname(fy1 - fy0)
    den1s <- if (d1deg) ctx$dmeanz1 else unname(fd1)
    den0s <- if (d0deg) ctx$dmeanz0 else unname(fd0)
    denom <- den1s - den0s
    late  <- num / denom

    # Order: eqips* eqy0 eqy1 num [eqd0] [eqd1] denom late
    blocks <- c(setup$blocks, list(
      make_glm_block(ctx, "y0", "gaussian", ones, ctx$y, 0, rw0, fy0),
      make_glm_block(ctx, "y1", "gaussian", ones, ctx$y, 1, rw1, fy1),
      make_contrast_block(ctx, "num",
                          pred_fun(ctx, "y1", "gaussian", ones),
                          pred_fun(ctx, "y0", "gaussian", ones), num)
    ))
    if (!d0deg) {
      blocks <- c(blocks, list(
        make_glm_block(ctx, "d0", "gaussian", ones, ctx$d, 0, rw0, fd0)))
    }
    if (!d1deg) {
      blocks <- c(blocks, list(
        make_glm_block(ctx, "d1", "gaussian", ones, ctx$d, 1, rw1, fd1)))
    }
    blocks <- c(blocks, list(
      make_contrast_block(ctx, "denom",
                          pred_fun(ctx, "d1", "gaussian", ones,
                                   if (d1deg) ctx$dmeanz1 else NULL),
                          pred_fun(ctx, "d0", "gaussian", ones,
                                   if (d0deg) ctx$dmeanz0 else NULL), denom),
      make_late_block(ctx, late)
    ))
  } else {
    # Unnormalized: raw IPW means (drlate_estimate_late.ado lines 495-591)
    w1raw <- ps$wt1
    w0raw <- ps$wt0
    y0s <- wmean(w0raw * ctx$y, w)
    num <- wmean(w1raw * ctx$y, w) - y0s
    den1s <- if (d1deg) ctx$dmeanz1 else wmean(w1raw * ctx$d, w)
    den0s <- if (d0deg) ctx$dmeanz0 else wmean(w0raw * ctx$d, w)
    denom <- den1s - den0s
    late  <- num / denom

    # Per-observation raw weights through the PS linear index
    ew1 <- function(theta, layout) z * rw1(theta, layout)
    ew0 <- function(theta, layout) (1 - z) * rw0(theta, layout)
    y <- ctx$y; d <- ctx$d

    # Order: eqips eqy0ipw num [eqd0ipw] denom late
    blocks <- c(setup$blocks, list(
      make_custom_block(ctx, "y0", y0s, function(theta, layout)
        ew0(theta, layout) * y - theta[layout$y0]),
      make_custom_block(ctx, "num", num, function(theta, layout)
        theta[layout$num] -
          (ew1(theta, layout) * y - ew0(theta, layout) * y))
    ))
    if (!d0deg) {
      blocks <- c(blocks, list(
        make_custom_block(ctx, "d0", den0s, function(theta, layout)
          ew0(theta, layout) * d - theta[layout$d0])))
    }
    denom_mf <- if (d0deg) {
      function(theta, layout)
        theta[layout$denom] - (ew1(theta, layout) * d - ctx$dmeanz0)
    } else if (d1deg) {
      function(theta, layout)
        theta[layout$denom] - (ctx$dmeanz1 - ew0(theta, layout) * d)
    } else {
      function(theta, layout)
        theta[layout$denom] -
          (ew1(theta, layout) * d - ew0(theta, layout) * d)
    }
    blocks <- c(blocks, list(
      make_custom_block(ctx, "denom", denom, denom_mf),
      make_late_block(ctx, late)
    ))
  }

  list(blocks = blocks,
       estimates = c(late = late, num = num, denom = denom))
}

#' LATE via augmented inverse probability weighting.
#' Stata: drlate_estimate_late.ado, METHOD: AIPW (lines 657-1250).
#' Outcome/treatment models are fitted without IPW weights; the augmentation
#' enters through scalar moment blocks. The normalized variant carries
#' explicit w1/w0 normalization parameters; the unnormalized variant has
#' per-case differences in which denom blocks exist (faithful to Stata).
#' @noRd
late_aipw <- function(ctx, ps) {
  w <- ctx$w; z <- ctx$z; y <- ctx$y; d <- ctx$d
  setup <- late_ps_setup(ctx, ps)
  rw1 <- setup$rw1
  rw0 <- setup$rw0
  ew1 <- function(theta, layout) z * rw1(theta, layout)        # z/p
  ew0 <- function(theta, layout) (1 - z) * rw0(theta, layout)  # (1-z)/(1-p)

  # Plain (non-IPW-weighted) regressions
  fy1 <- fit_wglm(y, ctx$Xo, ctx$omodel, w, z == 1)
  fy0 <- fit_wglm(y, ctx$Xo, ctx$omodel, w, z == 0)
  mu_y1 <- predict_glm(fy1, ctx$Xo, ctx$omodel)
  mu_y0 <- predict_glm(fy0, ctx$Xo, ctx$omodel)
  a1 <- fit_arm(ctx, d, ctx$Xt, ctx$tmodel, 1, w, ctx$dmeanz1)
  a0 <- fit_arm(ctx, d, ctx$Xt, ctx$tmodel, 0, w, ctx$dmeanz0)
  d0deg <- !is.null(a0$degenerate_value)
  d1deg <- !is.null(a1$degenerate_value)

  pred_y1 <- pred_fun(ctx, "y1", ctx$omodel, ctx$Xo)
  pred_y0 <- pred_fun(ctx, "y0", ctx$omodel, ctx$Xo)
  pred_d1 <- pred_fun(ctx, "d1", ctx$tmodel, ctx$Xt, a1$degenerate_value)
  pred_d0 <- pred_fun(ctx, "d0", ctx$tmodel, ctx$Xt, a0$degenerate_value)

  if (ctx$statnorm == "unnrm") {
    # AIPW influence terms at the fitted values
    num1s <- wmean((z * y - (z - ps$ps) * mu_y1) / ps$ps, w)
    num0s <- wmean(((1 - z) * y + (z - ps$ps) * mu_y0) / (1 - ps$ps), w)
    den1s <- wmean((z * d - (z - ps$ps) * a1$mu) / ps$ps, w)
    den0s <- wmean(((1 - z) * d + (z - ps$ps) * a0$mu) / (1 - ps$ps), w)
    num   <- num1s - num0s
    denom <- den1s - den0s
    late  <- num / denom

    # eqips* eqy0 eqy1 num0 num1 num then d-blocks per case, denom, late
    blocks <- c(setup$blocks, list(
      make_glm_block(ctx, "y0", ctx$omodel, ctx$Xo, y, 0, rw_one(ctx), fy0),
      make_glm_block(ctx, "y1", ctx$omodel, ctx$Xo, y, 1, rw_one(ctx), fy1),
      make_custom_block(ctx, "num0", num0s, function(theta, layout)
        theta[layout$num0] -
          ew0(theta, layout) * (y - pred_y0(theta, layout)) -
          pred_y0(theta, layout)),
      make_custom_block(ctx, "num1", num1s, function(theta, layout)
        theta[layout$num1] -
          ew1(theta, layout) * (y - pred_y1(theta, layout)) -
          pred_y1(theta, layout)),
      make_contrast_block(ctx, "num", term_param("num1"),
                          term_param("num0"), num)
    ))
    if (d0deg) {
      # eqd1 denom1 denom  (denom = denom1 - dmeanz0)
      blocks <- c(blocks, list(
        make_glm_block(ctx, "d1", ctx$tmodel, ctx$Xt, d, 1, rw_one(ctx),
                       a1$coefs),
        make_custom_block(ctx, "denom1", den1s, function(theta, layout)
          theta[layout$denom1] -
            ew1(theta, layout) * (d - pred_d1(theta, layout)) -
            pred_d1(theta, layout)),
        make_contrast_block(ctx, "denom", term_param("denom1"),
                            term_const(ctx$dmeanz0), denom)))
    } else if (d1deg) {
      blocks <- c(blocks, list(
        make_glm_block(ctx, "d0", ctx$tmodel, ctx$Xt, d, 0, rw_one(ctx),
                       a0$coefs),
        make_custom_block(ctx, "denom0", den0s, function(theta, layout)
          theta[layout$denom0] -
            ew0(theta, layout) * (d - pred_d0(theta, layout)) -
            pred_d0(theta, layout)),
        make_contrast_block(ctx, "denom", term_const(ctx$dmeanz1),
                            term_param("denom0"), denom)))
    } else if (ctx$ivmodel == "ipt") {
      # Interior with IPT: denom0/denom1 ARE separate parameters
      # (drlate_estimate_late.ado lines 954-963), unlike the non-IPT
      # interior below — replicated faithfully.
      blocks <- c(blocks, list(
        make_glm_block(ctx, "d0", ctx$tmodel, ctx$Xt, d, 0, rw_one(ctx),
                       a0$coefs),
        make_glm_block(ctx, "d1", ctx$tmodel, ctx$Xt, d, 1, rw_one(ctx),
                       a1$coefs),
        make_custom_block(ctx, "denom0", den0s, function(theta, layout)
          theta[layout$denom0] -
            ew0(theta, layout) * (d - pred_d0(theta, layout)) -
            pred_d0(theta, layout)),
        make_custom_block(ctx, "denom1", den1s, function(theta, layout)
          theta[layout$denom1] -
            ew1(theta, layout) * (d - pred_d1(theta, layout)) -
            pred_d1(theta, layout)),
        make_contrast_block(ctx, "denom", term_param("denom1"),
                            term_param("denom0"), denom)))
    } else {
      # Interior: no denom0/denom1 parameters; the denom moment is inline
      blocks <- c(blocks, list(
        make_glm_block(ctx, "d0", ctx$tmodel, ctx$Xt, d, 0, rw_one(ctx),
                       a0$coefs),
        make_glm_block(ctx, "d1", ctx$tmodel, ctx$Xt, d, 1, rw_one(ctx),
                       a1$coefs),
        make_custom_block(ctx, "denom", denom, function(theta, layout)
          theta[layout$denom] -
            ((ew1(theta, layout) * (d - pred_d1(theta, layout)) +
                pred_d1(theta, layout)) -
             (ew0(theta, layout) * (d - pred_d0(theta, layout)) +
                pred_d0(theta, layout))))))
    }
    blocks <- c(blocks, list(make_late_block(ctx, late)))
  } else {
    # Normalized: explicit w1/w0 normalization parameters
    w1s <- wmean(z / ps$ps, w)
    w0s <- wmean((1 - z) / (1 - ps$ps), w)
    omega1 <- (z / ps$ps) / w1s
    omega0 <- ((1 - z) / (1 - ps$ps)) / w0s
    num1s <- wmean(omega1 * (y - mu_y1) + mu_y1, w)
    num0s <- wmean(omega0 * (y - mu_y0) + mu_y0, w)
    den1s <- wmean(omega1 * (d - a1$mu) + a1$mu, w)
    den0s <- wmean(omega0 * (d - a0$mu) + a0$mu, w)
    num   <- num1s - num0s
    denom <- den1s - den0s
    late  <- num / denom

    # eqips eqy0 eqy1 eqw1 eqw0 num1 num0 num then d-blocks, denom, late
    blocks <- c(setup$blocks, list(
      make_glm_block(ctx, "y0", ctx$omodel, ctx$Xo, y, 0, rw_one(ctx), fy0),
      make_glm_block(ctx, "y1", ctx$omodel, ctx$Xo, y, 1, rw_one(ctx), fy1),
      make_custom_block(ctx, "w1", w1s, function(theta, layout)
        ew1(theta, layout) - theta[layout$w1]),
      make_custom_block(ctx, "w0", w0s, function(theta, layout)
        ew0(theta, layout) - theta[layout$w0]),
      make_custom_block(ctx, "num1", num1s, function(theta, layout)
        theta[layout$num1] -
          (ew1(theta, layout) / theta[layout$w1]) *
            (y - pred_y1(theta, layout)) -
          pred_y1(theta, layout)),
      make_custom_block(ctx, "num0", num0s, function(theta, layout)
        theta[layout$num0] -
          (ew0(theta, layout) / theta[layout$w0]) *
            (y - pred_y0(theta, layout)) -
          pred_y0(theta, layout)),
      make_contrast_block(ctx, "num", term_param("num1"),
                          term_param("num0"), num)
    ))
    if (d0deg && ctx$dmeanz1 != 1) {
      blocks <- c(blocks, list(
        make_glm_block(ctx, "d1", ctx$tmodel, ctx$Xt, d, 1, rw_one(ctx),
                       a1$coefs),
        make_custom_block(ctx, "denom1", den1s, function(theta, layout)
          theta[layout$denom1] -
            (ew1(theta, layout) / theta[layout$w1]) *
              (d - pred_d1(theta, layout)) -
            pred_d1(theta, layout)),
        make_contrast_block(ctx, "denom", term_param("denom1"),
                            term_const(ctx$dmeanz0), denom)))
    } else if (!d0deg && d1deg) {
      blocks <- c(blocks, list(
        make_glm_block(ctx, "d0", ctx$tmodel, ctx$Xt, d, 0, rw_one(ctx),
                       a0$coefs),
        make_custom_block(ctx, "denom0", den0s, function(theta, layout)
          theta[layout$denom0] -
            (ew0(theta, layout) / theta[layout$w0]) *
              (d - pred_d0(theta, layout)) -
            pred_d0(theta, layout)),
        make_contrast_block(ctx, "denom", term_const(ctx$dmeanz1),
                            term_param("denom0"), denom)))
    } else if (!d0deg && !d1deg) {
      blocks <- c(blocks, list(
        make_glm_block(ctx, "d0", ctx$tmodel, ctx$Xt, d, 0, rw_one(ctx),
                       a0$coefs),
        make_glm_block(ctx, "d1", ctx$tmodel, ctx$Xt, d, 1, rw_one(ctx),
                       a1$coefs),
        make_custom_block(ctx, "denom0", den0s, function(theta, layout)
          theta[layout$denom0] -
            (ew0(theta, layout) / theta[layout$w0]) *
              (d - pred_d0(theta, layout)) -
            pred_d0(theta, layout)),
        make_custom_block(ctx, "denom1", den1s, function(theta, layout)
          theta[layout$denom1] -
            (ew1(theta, layout) / theta[layout$w1]) *
              (d - pred_d1(theta, layout)) -
            pred_d1(theta, layout)),
        make_contrast_block(ctx, "denom", term_param("denom1"),
                            term_param("denom0"), denom)))
    } else {
      stop("unsupported degeneracy pattern for method = \"aipw\" with ",
           "normalized moments.", call. = FALSE)
    }
    blocks <- c(blocks, list(make_late_block(ctx, late)))
  }

  list(blocks = blocks,
       estimates = c(late = late, num = num, denom = denom))
}

#' Shared kappa numerator block: delta = E(ZY/p - (1-Z)Y/(1-p))
#' (kappalate eq_delta), expressed through the PS linear index.
#' @noRd
make_kappa_num_block <- function(ctx, setup, start) {
  z <- ctx$z; y <- ctx$y
  rw1 <- setup$rw1; rw0 <- setup$rw0
  make_custom_block(ctx, "num", start, function(theta, layout)
    z * rw1(theta, layout) * y - (1 - z) * rw0(theta, layout) * y -
      theta[layout$num])
}

#' LATE via unnormalized Abadie kappa weighting (kappalate tau_a).
#' Stata: kappalate.ado eq_delta + eq_gamma + eq_tau_a, gmm onestep
#' iterate(0) — point estimates are closed-form means at the fitted PS,
#' the stack exists only for the joint sandwich.
#' @noRd
late_kappa <- function(ctx, ps) {
  w <- ctx$w; z <- ctx$z; y <- ctx$y; d <- ctx$d
  setup <- late_ps_setup(ctx, ps)
  rw1 <- setup$rw1; rw0 <- setup$rw0

  num   <- wmean(ps$wt1 * y - ps$wt0 * y, w)
  denom <- wmean(1 - d * (1 - z) / (1 - ps$ps) - (1 - d) * z / ps$ps, w)
  late  <- num / denom

  # Order: eqips, num (delta), denom (gamma), late — as in kappalate.ado
  blocks <- c(setup$blocks, list(
    make_kappa_num_block(ctx, setup, num),
    make_custom_block(ctx, "denom", denom, function(theta, layout)
      1 - d * (1 - z) * rw0(theta, layout) -
        (1 - d) * z * rw1(theta, layout) - theta[layout$denom]),
    make_late_block(ctx, late)
  ))
  list(blocks = blocks,
       estimates = c(late = late, num = num, denom = denom))
}

#' LATE via the (1-D)-arm unnormalized kappa weighting (kappalate tau_a,0).
#' Stata: kappalate.ado eq_delta + eq_gamma0 + eq_tau_a0; gamma0 uses the
#' (D-1) contrast form, identical pointwise to (1-D)((1-Z)-(1-p))/(p(1-p)).
#' @noRd
late_kappa0 <- function(ctx, ps) {
  w <- ctx$w; z <- ctx$z; y <- ctx$y; d <- ctx$d
  setup <- late_ps_setup(ctx, ps)
  rw1 <- setup$rw1; rw0 <- setup$rw0

  num   <- wmean(ps$wt1 * y - ps$wt0 * y, w)
  denom <- wmean((d - 1) * (ps$wt1 - ps$wt0), w)
  late  <- num / denom

  blocks <- c(setup$blocks, list(
    make_kappa_num_block(ctx, setup, num),
    make_custom_block(ctx, "denom", denom, function(theta, layout)
      (d - 1) * (z * rw1(theta, layout) - (1 - z) * rw0(theta, layout)) -
        theta[layout$denom]),
    make_late_block(ctx, late)
  ))
  list(blocks = blocks,
       estimates = c(late = late, num = num, denom = denom))
}

#' LATE via normalized Abadie kappa weighting (kappalate tau_a,10):
#' delta1/gamma1 - delta0/gamma0, the contrast of kappa-weighted complier
#' potential-outcome means. Stata: kappalate.ado eq_delta1, eq_gamma1,
#' eq_delta0, eq_gamma0, eq_tau_a10 (same block order).
#' @noRd
late_kappa10 <- function(ctx, ps) {
  w <- ctx$w; z <- ctx$z; y <- ctx$y; d <- ctx$d
  setup <- late_ps_setup(ctx, ps)
  rw1 <- setup$rw1; rw0 <- setup$rw0

  kap1 <- ps$wt1 - ps$wt0       # z/p - (1-z)/(1-p) at the fitted PS
  num1s   <- wmean(d * kap1 * y, w)
  denom1s <- wmean(d * kap1, w)
  num0s   <- wmean((d - 1) * kap1 * y, w)
  denom0s <- wmean((d - 1) * kap1, w)
  late <- num1s / denom1s - num0s / denom0s

  contrast <- function(theta, layout)
    z * rw1(theta, layout) - (1 - z) * rw0(theta, layout)

  blocks <- c(setup$blocks, list(
    make_custom_block(ctx, "num1", num1s, function(theta, layout)
      d * contrast(theta, layout) * y - theta[layout$num1]),
    make_custom_block(ctx, "denom1", denom1s, function(theta, layout)
      d * contrast(theta, layout) - theta[layout$denom1]),
    make_custom_block(ctx, "num0", num0s, function(theta, layout)
      (d - 1) * contrast(theta, layout) * y - theta[layout$num0]),
    make_custom_block(ctx, "denom0", denom0s, function(theta, layout)
      (d - 1) * contrast(theta, layout) - theta[layout$denom0]),
    make_late_diff_block(ctx, late)
  ))
  list(blocks = blocks, estimates = c(late = late))
}
