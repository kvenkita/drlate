# Imai-Ratkovic overidentification balance test (latebalance overid).

test_that("balance_test() returns an overid test object", {
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age + educ,
                data = drlate_sim)
  bt <- balance_test(fit)
  expect_s3_class(bt, "drlate_balance_test")
  expect_true(all(c("statistic", "df", "p.value") %in% names(bt)))
  expect_gte(bt$statistic, 0)
  expect_equal(bt$df, ncol(fit$ctx$Xz))            # one condition per PS column
  expect_equal(bt$p.value,
               pchisq(bt$statistic, bt$df, lower.tail = FALSE))
  out <- paste(capture.output(print(bt)), collapse = "\n")
  expect_match(out, "balance")
})

test_that("balance_test() errors for cbps, ra, and lean fits", {
  d <- drlate_sim
  fit_cbps <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                     ivmodel = "cbps")
  expect_error(balance_test(fit_cbps), "cbps")
  fit_ra <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ 1, data = d,
                   method = "ra")
  expect_error(balance_test(fit_ra), "method = \"ra\"")
  fit_lean <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                     keep_data = FALSE)
  expect_error(balance_test(fit_lean), "keep_data")
})

test_that("balance_test() rejects under a misspecified propensity score", {
  set.seed(123)
  n <- 4000
  x <- runif(n, -2, 2)
  mk <- function(pz) {
    z <- rbinom(n, 1, pz)
    d <- rbinom(n, 1, plogis(-0.2 + 1.0 * z + 0.3 * x))
    y <- 0.5 * d + 0.4 * x + rnorm(n)
    data.frame(y, d, z, x)
  }
  dg <- mk(plogis(1.0 * x))        # well specified for a linear logit
  db <- mk(plogis(0.5 * x^3))      # monotone but nonlinear -> linear logit wrong
  fg <- drlate(y ~ x, d ~ x, z ~ x, data = dg)
  fb <- drlate(y ~ x, d ~ x, z ~ x, data = db)
  bg <- balance_test(fg)
  bb <- balance_test(fb)
  expect_gt(bb$statistic, bg$statistic)
  expect_lt(bb$p.value, 0.05)
})
