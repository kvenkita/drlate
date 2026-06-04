# LATT estimators: block dispatch per (method, normalization, ivmodel, case).
# Mirrors drlate_estimate_latt.ado. Distinguishing features relative to LATE:
# * Z=1 outcome/treatment "regressions" are intercept-only weighted means
#   (plain `regress`), regardless of omodel/tmodel.
# * Z=0 regressions are weighted by the ATT odds weight w * p/(1-p) (IPWRA,
#   IPW-normalized) or fitted plain (AIPW), and their moment conditions are
#   reweighted by exp(zhat) = p/(1-p).
# * Aggregates are means over the Z=1 subsample.
# Several Stata quirks are replicated deliberately and flagged inline.

#' @noRd
estimate_latt <- function(ctx, ps) {
  if (ctx$case == "bothdeg" && ctx$method != "ipwra") {
    stop("treatment is degenerate in both instrument arms; the effect of ",
         "Z on D is constant and no estimation is needed.", call. = FALSE)
  }
  switch(ctx$method,
    ipwra = latt_ipwra(ctx, ps),
    ra    = latt_ra(ctx),
    ipw   = latt_ipw(ctx, ps),
    aipw  = latt_aipw(ctx, ps),
    stop("method = \"", ctx$method, "\" is not implemented yet.",
         call. = FALSE)
  )
}

#' Instrument-PS block for LATT: logit score, or the single Z=0 IPT tilt
#' with its parameter equation named `zhat` so that the odds reweight
#' exp(zhat) refers to the tilted coefficients
#' (drlate_estimate_latt.ado lines 25-33: bips = bips0).
#' @noRd
latt_ps_block <- function(ctx, ps) {
  switch(ctx$ivmodel,
    logit = make_ps_logit_block(ctx, ps$bips),
    ipt   = make_ps_ipt0_block(ctx, ps$bips0, eq = "zhat"),
    stop("ivmodel = \"", ctx$ivmodel, "\" is not available with ",
         "estimand = \"latt\".", call. = FALSE)
  )
}

#' @noRd
latt_ones <- function(ctx) {
  matrix(1, ctx$n, 1L, dimnames = list(NULL, "(Intercept)"))
}

#' LATT via IPWRA. Stata: drlate_estimate_latt.ado lines 47-197.
#' @noRd
latt_ipwra <- function(ctx, ps) {
  w <- ctx$w; z <- ctx$z; y <- ctx$y; d <- ctx$d
  ones <- latt_ones(ctx)
  attw <- w * ps$ps / (1 - ps$ps)   # Z=0 fitting weight (ATT odds)
  s1 <- z == 1
  rwodds <- rw_odds(ctx)

  # Z=1: intercept-only weighted means; Z=0: odds-weighted models
  fy1 <- fit_wglm(y, ones, "gaussian", w, s1)
  fy0 <- fit_wglm(y, ctx$Xo, ctx$omodel, attw, !s1)
  mu_y0 <- predict_glm(fy0, ctx$Xo, ctx$omodel)
  a1 <- fit_arm_latt1(ctx, d, w, ctx$dmeanz1)
  a0 <- if (ctx$dmeanz0 %in% c(0, 1)) {
    list(coefs = NULL, mu = rep(ctx$dmeanz0, ctx$n),
         degenerate_value = ctx$dmeanz0)
  } else {
    b <- fit_wglm(d, ctx$Xt, ctx$tmodel, attw, !s1)
    list(coefs = b, mu = predict_glm(b, ctx$Xt, ctx$tmodel),
         degenerate_value = NULL)
  }
  d0deg <- !is.null(a0$degenerate_value)
  d1deg <- !is.null(a1$degenerate_value)

  num   <- unname(fy1) - wmean(mu_y0[s1], w[s1])
  denom <- wmean(a1$mu[s1], w[s1]) - wmean(a0$mu[s1], w[s1])
  late  <- num / denom

  pred_y1 <- pred_fun(ctx, "y1", "gaussian", ones)
  pred_y0 <- pred_fun(ctx, "y0", ctx$omodel, ctx$Xo)
  pred_d1 <- pred_fun(ctx, "d1", "gaussian", ones, a1$degenerate_value)
  pred_d0 <- pred_fun(ctx, "d0", ctx$tmodel, ctx$Xt, a0$degenerate_value)

  # eqips eqy0 eqy1 num [eqd0] [eqd1] denom late
  blocks <- list(
    latt_ps_block(ctx, ps),
    make_glm_block(ctx, "y0", ctx$omodel, ctx$Xo, y, 0, rwodds, fy0),
    make_glm_block(ctx, "y1", "gaussian", ones, y, 1, rw_one(ctx), fy1),
    make_contrast_block(ctx, "num", pred_y1, pred_y0, num, arm = 1)
  )
  if (!d0deg) {
    blocks <- c(blocks, list(
      make_glm_block(ctx, "d0", ctx$tmodel, ctx$Xt, d, 0, rwodds, a0$coefs)))
  }
  if (!d1deg) {
    blocks <- c(blocks, list(
      make_glm_block(ctx, "d1", "gaussian", ones, d, 1, rw_one(ctx),
                     a1$coefs)))
  }
  t1 <- if (d1deg) term_const(ctx$dmeanz1) else pred_d1
  t0 <- if (d0deg) term_const(ctx$dmeanz0) else pred_d0
  blocks <- c(blocks, list(
    make_contrast_block(ctx, "denom", t1, t0, denom, arm = 1),
    make_late_block(ctx, late)
  ))

  list(blocks = blocks,
       estimates = c(late = late, num = num, denom = denom))
}

#' Z=1 arm intercept-only mean fit (plain `regress`), or a constant
#' @noRd
fit_arm_latt1 <- function(ctx, v, w, dmean) {
  if (dmean %in% c(0, 1)) {
    return(list(coefs = NULL, mu = rep(dmean, ctx$n),
                degenerate_value = dmean))
  }
  ones <- latt_ones(ctx)
  b <- fit_wglm(v, ones, "gaussian", w, ctx$z == 1)
  list(coefs = b, mu = rep(unname(b), ctx$n), degenerate_value = NULL)
}

#' LATT via RA. Stata: drlate_estimate_latt.ado lines 871-1016.
#' Quirk replicated: the denominator moment is NOT restricted to the Z=1
#' arm even though the estimate is a Z=1 mean (the num moment IS restricted);
#' with iterate(0) this leaves the point estimate untouched and only the
#' variance reflects the full-sample moment, exactly as in Stata.
#' @noRd
latt_ra <- function(ctx) {
  w <- ctx$w; z <- ctx$z; y <- ctx$y; d <- ctx$d
  ones <- latt_ones(ctx)
  s1 <- z == 1

  fy1 <- fit_wglm(y, ones, "gaussian", w, s1)
  fy0 <- fit_wglm(y, ctx$Xo, ctx$omodel, w, !s1)
  mu_y0 <- predict_glm(fy0, ctx$Xo, ctx$omodel)
  a1 <- fit_arm_latt1(ctx, d, w, ctx$dmeanz1)
  a0 <- fit_arm(ctx, d, ctx$Xt, ctx$tmodel, 0, w, ctx$dmeanz0)
  d0deg <- !is.null(a0$degenerate_value)
  d1deg <- !is.null(a1$degenerate_value)

  num   <- unname(fy1) - wmean(mu_y0[s1], w[s1])
  denom <- wmean(a1$mu[s1], w[s1]) - wmean(a0$mu[s1], w[s1])
  late  <- num / denom

  pred_y1 <- pred_fun(ctx, "y1", "gaussian", ones)
  pred_y0 <- pred_fun(ctx, "y0", ctx$omodel, ctx$Xo)
  pred_d1 <- pred_fun(ctx, "d1", "gaussian", ones, a1$degenerate_value)
  pred_d0 <- pred_fun(ctx, "d0", ctx$tmodel, ctx$Xt, a0$degenerate_value)

  # eqy0 eqy1 num [eqd0] [eqd1] denom late
  blocks <- list(
    make_glm_block(ctx, "y0", ctx$omodel, ctx$Xo, y, 0, rw_one(ctx), fy0),
    make_glm_block(ctx, "y1", "gaussian", ones, y, 1, rw_one(ctx), fy1),
    make_contrast_block(ctx, "num", pred_y1, pred_y0, num, arm = 1)
  )
  if (!d0deg) {
    blocks <- c(blocks, list(
      make_glm_block(ctx, "d0", ctx$tmodel, ctx$Xt, d, 0, rw_one(ctx),
                     a0$coefs)))
  }
  if (!d1deg) {
    blocks <- c(blocks, list(
      make_glm_block(ctx, "d1", "gaussian", ones, d, 1, rw_one(ctx),
                     a1$coefs)))
  }
  t1 <- if (d1deg) term_const(ctx$dmeanz1) else pred_d1
  t0 <- if (d0deg) term_const(ctx$dmeanz0) else pred_d0
  blocks <- c(blocks, list(
    make_contrast_block(ctx, "denom", t1, t0, denom),   # full-sample (quirk)
    make_late_block(ctx, late)
  ))

  list(blocks = blocks,
       estimates = c(late = late, num = num, denom = denom))
}

#' LATT via IPW. Stata: drlate_estimate_latt.ado lines 202-419.
#' @noRd
latt_ipw <- function(ctx, ps) {
  w <- ctx$w; z <- ctx$z; y <- ctx$y; d <- ctx$d
  ones <- latt_ones(ctx)
  s1 <- z == 1
  rwodds <- rw_odds(ctx)
  d0deg <- ctx$dmeanz0 %in% c(0, 1)
  d1deg <- ctx$dmeanz1 %in% c(0, 1)

  if (ctx$statnorm == "nrm") {
    attw <- w * ps$ps / (1 - ps$ps)
    fy1 <- fit_wglm(y, ones, "gaussian", w, s1)
    fy0 <- fit_wglm(y, ones, "gaussian", attw, !s1)
    fd1 <- if (d1deg) NULL else fit_wglm(d, ones, "gaussian", w, s1)
    fd0 <- if (d0deg) NULL else fit_wglm(d, ones, "gaussian", attw, !s1)

    num <- unname(fy1 - fy0)
    den1s <- if (d1deg) ctx$dmeanz1 else unname(fd1)
    den0s <- if (d0deg) ctx$dmeanz0 else unname(fd0)
    denom <- den1s - den0s
    late  <- num / denom

    # eqips eqy0 eqy1 num [eqd0] [eqd1] denom late
    blocks <- list(
      latt_ps_block(ctx, ps),
      make_glm_block(ctx, "y0", "gaussian", ones, y, 0, rwodds, fy0),
      make_glm_block(ctx, "y1", "gaussian", ones, y, 1, rw_one(ctx), fy1),
      make_contrast_block(ctx, "num",
                          pred_fun(ctx, "y1", "gaussian", ones),
                          pred_fun(ctx, "y0", "gaussian", ones), num,
                          arm = 1)
    )
    if (!d0deg) {
      blocks <- c(blocks, list(
        make_glm_block(ctx, "d0", "gaussian", ones, d, 0, rwodds, fd0)))
    }
    if (!d1deg) {
      blocks <- c(blocks, list(
        make_glm_block(ctx, "d1", "gaussian", ones, d, 1, rw_one(ctx), fd1)))
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
    # Unnormalized (lines 308-418): scaled by the treated share w1
    w1s <- wmean(z, w)
    omega1 <- z / w1s
    omega0 <- ((1 - z) * ps$ps / (1 - ps$ps)) / w1s
    num1s <- wmean(omega1 * y, w)
    num0s <- wmean(omega0 * y, w)
    den1s <- if (d1deg) ctx$dmeanz1 else wmean(omega1 * d, w)
    den0s <- if (d0deg) ctx$dmeanz0 else wmean(omega0 * d, w)
    num   <- num1s - num0s
    denom <- den1s - den0s
    late  <- num / denom

    # eqips eqw1 eqy0 eqy1 num [eqd0] [eqd1] denom late
    blocks <- list(
      latt_ps_block(ctx, ps),
      make_custom_block(ctx, "w1", w1s, function(theta, layout)
        z - theta[layout$w1]),
      make_custom_block(ctx, "y0", num0s, function(theta, layout)
        (1 - z) * rwodds(theta, layout) * y / theta[layout$w1] -
          theta[layout$y0]),
      make_custom_block(ctx, "y1", num1s, function(theta, layout)
        z * (y / theta[layout$w1]) - theta[layout$y1]),
      make_contrast_block(ctx, "num", term_param("y1"), term_param("y0"),
                          num)
    )
    if (!d0deg) {
      blocks <- c(blocks, list(
        make_custom_block(ctx, "d0", den0s, function(theta, layout)
          (1 - z) * rwodds(theta, layout) * d / theta[layout$w1] -
            theta[layout$d0])))
    }
    if (!d1deg) {
      blocks <- c(blocks, list(
        make_custom_block(ctx, "d1", den1s, function(theta, layout)
          z * (d / theta[layout$w1]) - theta[layout$d1])))
    }
    t1 <- if (d1deg) term_const(ctx$dmeanz1) else term_param("d1")
    t0 <- if (d0deg) term_const(ctx$dmeanz0) else term_param("d0")
    blocks <- c(blocks, list(
      make_contrast_block(ctx, "denom", t1, t0, denom),
      make_late_block(ctx, late)
    ))
  }

  list(blocks = blocks,
       estimates = c(late = late, num = num, denom = denom))
}

#' LATT via AIPW. Stata: drlate_estimate_latt.ado lines 424-865.
#' Quirk replicated: in the unnormalized variant Stata computes the treated
#' share w1 WITHOUT weights (line 431), unlike the normalized variant.
#' @noRd
latt_aipw <- function(ctx, ps) {
  w <- ctx$w; z <- ctx$z; y <- ctx$y; d <- ctx$d
  ones <- latt_ones(ctx)
  s1 <- z == 1
  rwodds <- rw_odds(ctx)
  attw_raw <- (1 - z) * ps$ps / (1 - ps$ps)

  # Plain-weight regressions on the Z=0 arm
  fy0 <- fit_wglm(y, ctx$Xo, ctx$omodel, w, !s1)
  mu_y0 <- predict_glm(fy0, ctx$Xo, ctx$omodel)
  a0 <- fit_arm(ctx, d, ctx$Xt, ctx$tmodel, 0, w, ctx$dmeanz0)
  d0deg <- !is.null(a0$degenerate_value)
  d1deg <- ctx$dmeanz1 %in% c(0, 1)

  pred_y0 <- pred_fun(ctx, "y0", ctx$omodel, ctx$Xo)
  pred_d0 <- pred_fun(ctx, "d0", ctx$tmodel, ctx$Xt, a0$degenerate_value)

  if (ctx$statnorm == "unnrm") {
    w1s <- mean(z)                      # unweighted (Stata line 431)
    num1s <- wmean((z / w1s) * y, w)
    num0s <- wmean(attw_raw * (y - mu_y0) / w1s, w) +
             wmean(z * mu_y0 / w1s, w)
    den1s <- if (d1deg) ctx$dmeanz1 else wmean((z / w1s) * d, w)
    den0s <- if (d0deg) ctx$dmeanz0 else
             wmean(attw_raw * (d - a0$mu) / w1s, w) +
             wmean(z * a0$mu / w1s, w)
    num   <- num1s - num0s
    denom <- den1s - den0s
    late  <- num / denom

    # eqips eqw1 eqy0 eqy1 eqnum0 eqnum [d-blocks per case] eqdenom late
    blocks <- list(
      latt_ps_block(ctx, ps),
      make_custom_block(ctx, "w1", w1s, function(theta, layout)
        z - theta[layout$w1]),
      make_glm_block(ctx, "y0", ctx$omodel, ctx$Xo, y, 0, rw_one(ctx), fy0),
      make_glm_block(ctx, "y1", "gaussian", ones, y, 1, rw_one(ctx),
                     stats::setNames(num1s, "(Intercept)")),
      make_custom_block(ctx, "num0", num0s, function(theta, layout)
        theta[layout$num0] -
          (rwodds(theta, layout) * (1 - z) * (y - pred_y0(theta, layout)) /
             theta[layout$w1] +
           z * pred_y0(theta, layout) / theta[layout$w1])),
      make_custom_block(ctx, "num", num, function(theta, layout)
        theta[layout$num] - (theta[layout$y1] - theta[layout$num0]))
    )
    if (d0deg && ctx$dmeanz1 != 1) {
      blocks <- c(blocks, list(
        make_glm_block(ctx, "d1", "gaussian", ones, d, 1, rw_one(ctx),
                       stats::setNames(den1s, "(Intercept)")),
        make_contrast_block(ctx, "denom", term_param("d1"),
                            term_const(ctx$dmeanz0), denom)))
    } else if (!d0deg && d1deg) {
      blocks <- c(blocks, list(
        make_glm_block(ctx, "d0", ctx$tmodel, ctx$Xt, d, 0, rw_one(ctx),
                       a0$coefs),
        make_custom_block(ctx, "denom0", den0s, function(theta, layout)
          theta[layout$denom0] -
            (rwodds(theta, layout) * (1 - z) *
               (d - pred_d0(theta, layout)) / theta[layout$w1] +
             z * pred_d0(theta, layout) / theta[layout$w1])),
        make_contrast_block(ctx, "denom", term_const(ctx$dmeanz1),
                            term_param("denom0"), denom)))
    } else if (!d0deg && !d1deg) {
      blocks <- c(blocks, list(
        make_glm_block(ctx, "d0", ctx$tmodel, ctx$Xt, d, 0, rw_one(ctx),
                       a0$coefs),
        make_glm_block(ctx, "d1", "gaussian", ones, d, 1, rw_one(ctx),
                       stats::setNames(den1s, "(Intercept)")),
        make_custom_block(ctx, "denom0", den0s, function(theta, layout)
          theta[layout$denom0] -
            (rwodds(theta, layout) * (1 - z) *
               (d - pred_d0(theta, layout)) / theta[layout$w1] +
             z * pred_d0(theta, layout) / theta[layout$w1])),
        make_contrast_block(ctx, "denom", term_param("d1"),
                            term_param("denom0"), denom)))
    } else {
      stop("unsupported degeneracy pattern for method = \"aipw\" with ",
           "estimand = \"latt\".", call. = FALSE)
    }
    blocks <- c(blocks, list(make_late_block(ctx, late)))
  } else {
    # Normalized (lines 662-865)
    w1s <- wmean(z, w)
    wnorms <- wmean(attw_raw, w)
    num1s <- wmean((z / w1s) * y, w)
    num0s <- wmean(attw_raw * (y - mu_y0) / wnorms, w) +
             wmean(z * mu_y0 / w1s, w)
    den1s <- if (d1deg) ctx$dmeanz1 else wmean((z / w1s) * d, w)
    den0s <- if (d0deg) ctx$dmeanz0 else
             wmean(attw_raw * (d - a0$mu) / wnorms, w) +
             wmean(z * a0$mu / w1s, w)
    num   <- num1s - num0s
    denom <- den1s - den0s
    late  <- num / denom

    # eqips eqw1 eqwnorm eqy0 eqy1 eqnum0 eqnum [d-blocks] eqdenom late
    blocks <- list(
      latt_ps_block(ctx, ps),
      make_custom_block(ctx, "w1", w1s, function(theta, layout)
        z - theta[layout$w1]),
      make_custom_block(ctx, "wnorm", wnorms, function(theta, layout)
        (1 - z) * rwodds(theta, layout) - theta[layout$wnorm]),
      make_glm_block(ctx, "y0", ctx$omodel, ctx$Xo, y, 0, rw_one(ctx), fy0),
      make_glm_block(ctx, "y1", "gaussian", ones, y, 1, rw_one(ctx),
                     stats::setNames(num1s, "(Intercept)")),
      make_custom_block(ctx, "num0", num0s, function(theta, layout)
        theta[layout$num0] -
          (rwodds(theta, layout) * (1 - z) * (y - pred_y0(theta, layout)) /
             theta[layout$wnorm] +
           z * pred_y0(theta, layout) / theta[layout$w1])),
      make_custom_block(ctx, "num", num, function(theta, layout)
        theta[layout$num] - (theta[layout$y1] - theta[layout$num0]))
    )
    if (d0deg && ctx$dmeanz1 != 1) {
      blocks <- c(blocks, list(
        make_glm_block(ctx, "d1", "gaussian", ones, d, 1, rw_one(ctx),
                       stats::setNames(den1s, "(Intercept)")),
        make_contrast_block(ctx, "denom", term_param("d1"),
                            term_const(ctx$dmeanz0), denom)))
    } else if (!d0deg && d1deg) {
      blocks <- c(blocks, list(
        make_glm_block(ctx, "d0", ctx$tmodel, ctx$Xt, d, 0, rw_one(ctx),
                       a0$coefs),
        make_custom_block(ctx, "denom0", den0s, function(theta, layout)
          theta[layout$denom0] -
            (rwodds(theta, layout) * (1 - z) *
               (d - pred_d0(theta, layout)) / theta[layout$wnorm] +
             z * pred_d0(theta, layout) / theta[layout$w1])),
        make_contrast_block(ctx, "denom", term_const(ctx$dmeanz1),
                            term_param("denom0"), denom)))
    } else if (!d0deg && !d1deg) {
      blocks <- c(blocks, list(
        make_glm_block(ctx, "d0", ctx$tmodel, ctx$Xt, d, 0, rw_one(ctx),
                       a0$coefs),
        make_glm_block(ctx, "d1", "gaussian", ones, d, 1, rw_one(ctx),
                       stats::setNames(den1s, "(Intercept)")),
        make_custom_block(ctx, "denom0", den0s, function(theta, layout)
          theta[layout$denom0] -
            (rwodds(theta, layout) * (1 - z) *
               (d - pred_d0(theta, layout)) / theta[layout$wnorm] +
             z * pred_d0(theta, layout) / theta[layout$w1])),
        make_contrast_block(ctx, "denom", term_param("d1"),
                            term_param("denom0"), denom)))
    } else {
      stop("unsupported degeneracy pattern for method = \"aipw\" with ",
           "estimand = \"latt\".", call. = FALSE)
    }
    blocks <- c(blocks, list(make_late_block(ctx, late)))
  }

  list(blocks = blocks,
       estimates = c(late = late, num = num, denom = denom))
}
