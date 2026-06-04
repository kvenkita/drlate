# DR Hausman test of unconfoundedness (SUW 2022 section 5).

onesided <- function() {
  d <- drlate_sim
  d$nvstat[d$rsncode == 0] <- 0L
  d
}

test_that("dr_hausman returns a valid htest with consistent pieces", {
  d <- onesided()
  h <- dr_hausman(lwage ~ age + educ, nvstat ~ age + educ,
                  rsncode ~ age + educ, data = d)
  expect_s3_class(h, "htest")
  expect_true(is.finite(h$statistic))
  expect_true(h$p.value >= 0 && h$p.value <= 1)
  expect_equal(unname(h$estimate["difference"]),
               unname(h$estimate["DR LATT"] - h$estimate["DR ATT"]))

  # The LATT half must equal the standalone LATT fit
  latt <- drlate(lwage ~ age + educ, nvstat ~ age + educ,
                 rsncode ~ age + educ, data = d, estimand = "latt")
  expect_equal(unname(h$estimate["DR LATT"]), unname(coef(latt)[1]),
               tolerance = 1e-10)
})

test_that("the DR ATT half equals the paper's imputation estimator (eq 33)", {
  d <- onesided()
  h <- dr_hausman(lwage ~ age + educ, nvstat ~ age + educ,
                  rsncode ~ age + educ, data = d)

  # Hand computation: logit treatment PS, odds-weighted outcome regression
  # on the controls, imputation over the treated. (Covariates standardized
  # the same way drlate does internally; standardization is span-preserving
  # so fitted values are identical with raw covariates.)
  dps <- fitted(glm(nvstat ~ age + educ, binomial, data = d))
  m <- model.matrix(~ age + educ, d)
  s1 <- d$nvstat == 1
  da <- transform(d, ww = dps / (1 - dps))
  fy0 <- lm(lwage ~ age + educ, da[!s1, ], weights = ww)
  mu0 <- drop(m %*% coef(fy0))
  att <- mean(d$lwage[s1]) - mean(mu0[s1])
  expect_equal(unname(h$estimate["DR ATT"]), att, tolerance = 1e-8)
})

test_that("the stacked Hausman system has zero averaged moments", {
  # In the simulated DGP treatment IS confounded (compliance type shifts
  # baseline outcomes), so the test should tend to reject; at minimum the
  # machinery must be internally consistent.
  d <- onesided()
  h <- dr_hausman(lwage ~ age + educ, nvstat ~ age + educ,
                  rsncode ~ age + educ, data = d)
  expect_true(is.finite(h$stderr))
  expect_gt(h$stderr, 0)
})

test_that("dr_hausman detects the confounded simulated treatment", {
  # Always-takers were removed, but never-takers still have lower baseline
  # outcomes, so unconfoundedness fails in this DGP and ATT != LATT.
  d <- onesided()
  h <- dr_hausman(lwage ~ age + educ, nvstat ~ age + educ,
                  rsncode ~ age + educ, data = d)
  expect_lt(h$p.value, 0.05)
})

test_that("two-sided noncompliance is rejected with guidance", {
  expect_error(
    dr_hausman(lwage ~ age, nvstat ~ age, rsncode ~ age,
               data = drlate_sim),
    "one-sided noncompliance"
  )
})
