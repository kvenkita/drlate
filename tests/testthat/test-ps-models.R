# CBPS and IPT instrument propensity score models.

make_ctx <- function(data, ivmodel, estimand = "late", method = "ipwra") {
  build_ctx(lwage ~ age + educ, nvstat ~ age + educ, rsncode ~ age + educ,
            data, omodel = "linear", tmodel = "logit",
            ivmodel = ivmodel, method = method, estimand = estimand,
            normalized = TRUE)
}

test_that("CBPS solves the balancing moment exactly", {
  ctx <- make_ctx(drlate_sim, "cbps")
  ps <- fit_ps(ctx)
  bal <- colSums(ctx$w * (ctx$z / ps$ps - (1 - ctx$z) / (1 - ps$ps)) *
                   ctx$Xz) / sum(ctx$w)
  expect_lt(max(abs(bal)), 1e-9)
})

test_that("IPT solves both tilting moments exactly", {
  ctx <- make_ctx(drlate_sim, "ipt")
  ps <- fit_ps(ctx)
  t1 <- colSums(ctx$w * (ctx$z / ps$ps1 - 1) * ctx$Xz) / sum(ctx$w)
  t0 <- colSums(ctx$w * ((1 - ctx$z) / (1 - ps$ps0) - 1) * ctx$Xz) /
        sum(ctx$w)
  expect_lt(max(abs(t1)), 1e-9)
  expect_lt(max(abs(t0)), 1e-9)
  # IPT weights are exactly normalized
  expect_equal(wmean(ps$wt1, ctx$w), 1, tolerance = 1e-9)
  expect_equal(wmean(ps$wt0, ctx$w), 1, tolerance = 1e-9)
})

test_that("CBPS works with LATE ipwra/ipw/aipw and is rejected for LATT", {
  d <- drlate_sim
  for (me in c("ipwra", "aipw")) {
    fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                  ivmodel = "cbps", method = me)
    expect_true(all(is.finite(coef(fit))), label = paste("cbps", me))
    expect_true(all(is.finite(sqrt(diag(fit$vcov3)))))
  }
  fit <- drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age, data = d,
                ivmodel = "cbps", method = "ipw")
  expect_true(all(is.finite(coef(fit))))
  expect_error(
    drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
           ivmodel = "cbps", estimand = "latt"),
    "not available"
  )
})

test_that("IPT works for LATE and LATT and switches to unnormalized", {
  d <- drlate_sim
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                ivmodel = "ipt")
  expect_true(all(is.finite(coef(fit))))
  # ipw + ipt: normalize-check switches to unnrm automatically
  fit_ipw <- drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age, data = d,
                    ivmodel = "ipt", method = "ipw")
  expect_equal(fit_ipw$statnorm, "unnrm")
  expect_true(all(is.finite(coef(fit_ipw))))
  fit_aipw <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                     ivmodel = "ipt", method = "aipw")
  expect_equal(fit_aipw$statnorm, "unnrm")
  expect_true(all(is.finite(coef(fit_aipw))))
  fit_latt <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                     ivmodel = "ipt", estimand = "latt")
  expect_true(all(is.finite(coef(fit_latt))))
  expect_true(all(is.finite(sqrt(diag(fit_latt$vcov3)))))
})

test_that("IPT LATE-IPWRA system is square with zero averaged moments", {
  ctx <- make_ctx(drlate_sim, "ipt")
  ps <- fit_ps(ctx)
  est <- estimate_late(ctx, ps)
  sys <- assemble_system(est$blocks)
  gbar <- colMeans(ctx$w * sys$g(sys$theta0))
  expect_lt(max(abs(gbar)), 1e-7)
  expect_equal(length(sys$theta0), ncol(sys$g(sys$theta0)))
  # two PS equations in the stack
  expect_true(all(c("zhat1", "zhat0") %in% names(sys$layout)))
})

# ---------------------------------------------------------------------------
# Golden fixtures
# ---------------------------------------------------------------------------

golden_ps <- function(id, ..., estimand = "late") {
  test_that(paste0("matches Stata: ", id), {
    skip_if_no_fixture(id)
    fit <- drlate(..., data = sipp_data(), estimand = estimand)
    expect_matches_fixture(fit, id)
  })
}

golden_ps("late_ipwra_cbps", lwage ~ age_5, nvstat ~ age_5,
          rsncode ~ age_5, ivmodel = "cbps")
golden_ps("late_ipw_cbps", lwage ~ 1, nvstat ~ 1, rsncode ~ age_5,
          ivmodel = "cbps", method = "ipw")
golden_ps("late_aipw_cbps", lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5,
          ivmodel = "cbps", method = "aipw")
golden_ps("late_ipwra_ipt", lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5,
          ivmodel = "ipt")
golden_ps("late_ipw_ipt", lwage ~ 1, nvstat ~ 1, rsncode ~ age_5,
          ivmodel = "ipt", method = "ipw")
golden_ps("late_aipw_ipt", lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5,
          ivmodel = "ipt", method = "aipw")
golden_ps("latt_ipwra_ipt", lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5,
          ivmodel = "ipt", estimand = "latt")
