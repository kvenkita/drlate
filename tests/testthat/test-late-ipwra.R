# LATE x IPWRA: internal consistency + golden fixtures.

test_that("IPWRA point estimate equals the independent sequential computation", {
  d <- drlate_sim
  fit <- drlate(lwage ~ age + educ, nvstat ~ age + educ, rsncode ~ age + educ,
                data = d)

  ps <- fitted(glm(rsncode ~ age + educ, binomial, data = d))
  m <- model.matrix(~ age + educ, d)
  s1 <- d$rsncode == 1
  da <- transform(d, ww1 = 1 / ps, ww0 = 1 / (1 - ps))
  fy1 <- lm(lwage ~ age + educ, da[s1, ], weights = ww1)
  fy0 <- lm(lwage ~ age + educ, da[!s1, ], weights = ww0)
  fd1 <- suppressWarnings(glm(nvstat ~ age + educ, quasibinomial, da[s1, ],
                              weights = ww1))
  fd0 <- suppressWarnings(glm(nvstat ~ age + educ, quasibinomial, da[!s1, ],
                              weights = ww0))
  num <- mean(m %*% coef(fy1)) - mean(m %*% coef(fy0))
  den <- mean(plogis(m %*% coef(fd1))) - mean(plogis(m %*% coef(fd0)))

  expect_equal(unname(coef(fit)), c(num / den, num, den), tolerance = 1e-10)
})

test_that("IPWRA recovers the simulated LATE within sampling error", {
  fit <- drlate(lwage ~ age + educ, nvstat ~ age + educ, rsncode ~ age + educ,
                data = drlate_sim)
  ci <- confint(fit)["LATE: D on Y", ]
  expect_gt(0.5, ci[1])   # true complier effect in the DGP is 0.5
  expect_lt(0.5, ci[2])
})

test_that("moment system is square and zero at the estimates", {
  d <- drlate_sim
  fit <- drlate(lwage ~ age + educ, nvstat ~ age + educ, rsncode ~ age + educ,
                data = d)
  # All stacked moment conditions must average to (numerically) zero at
  # theta0; this is what makes iterate(0) exact in the Stata original.
  ctx <- build_ctx(lwage ~ age + educ, nvstat ~ age + educ,
                   rsncode ~ age + educ, d,
                   omodel = "gaussian", tmodel = "binomial",
                   ivmodel = "logit", method = "ipwra", estimand = "late",
                   normalized = TRUE)
  ps <- fit_ps(ctx)
  est <- estimate_late(ctx, ps)
  sys <- assemble_system(est$blocks)
  gbar <- colMeans(ctx$w * sys$g(sys$theta0))
  expect_lt(max(abs(gbar)), 1e-8)
  expect_equal(length(sys$theta0), ncol(sys$g(sys$theta0)))
})

test_that("validation errors mirror the Stata package", {
  d <- drlate_sim
  d$badz <- d$rsncode + 1
  expect_error(
    drlate(lwage ~ age, nvstat ~ age, badz ~ age, data = d),
    "must be binary"
  )
  expect_error(
    drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
           omodel = "logit"),
    "must be binary"
  )
  d$negy <- d$lwage - 10
  expect_error(
    drlate(negy ~ age, nvstat ~ age, rsncode ~ age, data = d,
           omodel = "poisson"),
    "non-negative"
  )
  expect_error(
    drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
           method = "ipw"),
    "not allowed"
  )
  expect_error(
    drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
           method = "ra"),
    "not allowed"
  )
  expect_error(
    drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
           ivmodel = "cbps", estimand = "latt"),
    "not available"
  )
})

# ---------------------------------------------------------------------------
# Golden fixtures (require running inst/stata/make-fixtures.do in Stata)
# ---------------------------------------------------------------------------

test_that("matches Stata: late_ipwra_lin_logit_logit", {
  skip_if_no_fixture("late_ipwra_lin_logit_logit")
  d <- sipp_data()
  fit <- drlate(lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5, data = d)
  expect_matches_fixture(fit, "late_ipwra_lin_logit_logit")
})

test_that("matches Stata: late_ipwra_lin_lin_logit", {
  skip_if_no_fixture("late_ipwra_lin_lin_logit")
  d <- sipp_data()
  fit <- drlate(lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5, data = d,
                tmodel = "linear")
  expect_matches_fixture(fit, "late_ipwra_lin_lin_logit")
})

test_that("matches Stata: late_ipwra_lin_pois_logit", {
  skip_if_no_fixture("late_ipwra_lin_pois_logit")
  d <- sipp_data()
  fit <- drlate(lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5, data = d,
                tmodel = "poisson")
  expect_matches_fixture(fit, "late_ipwra_lin_pois_logit")
})

test_that("matches Stata: late_ipwra_logit_logit_logit", {
  skip_if_no_fixture("late_ipwra_logit_logit_logit")
  d <- sipp_data()
  fit <- drlate(hiwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5, data = d,
                omodel = "logit")
  expect_matches_fixture(fit, "late_ipwra_logit_logit_logit")
})

test_that("matches Stata: late_ipwra_pois_logit_logit", {
  skip_if_no_fixture("late_ipwra_pois_logit_logit")
  d <- sipp_data()
  fit <- drlate(kwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5, data = d,
                omodel = "poisson")
  expect_matches_fixture(fit, "late_ipwra_pois_logit_logit")
})

test_that("matches Stata: late_ipwra_multix", {
  skip_if_no_fixture("late_ipwra_multix")
  d <- sipp_data()
  fit <- drlate(lwage ~ age_5 + educ, nvstat ~ age_5 + educ,
                rsncode ~ age_5 + educ, data = d)
  expect_matches_fixture(fit, "late_ipwra_multix")
})
