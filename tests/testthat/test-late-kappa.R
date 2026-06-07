# Kappa-weighting estimators (Sloczynski, Uysal & Wooldridge 2025, JBES;
# Stata kappalate): tau_a, tau_a,0, tau_a,10.

test_that("kappa methods validate inputs", {
  d <- drlate_sim
  expect_error(
    drlate(lwage ~ age, nvstat ~ 1, rsncode ~ age + educ, data = d,
           method = "kappa"),
    "covariates are not allowed")
  expect_error(
    drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d,
           method = "kappa", estimand = "latt"),
    "\"late\" only")
  expect_error(
    drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d,
           method = "kappa", ivmodel = "ipt"),
    "ipt.*not available")
  expect_error(
    drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d,
           method = "kappa10", ivmodel = "cbps"),
    "cbps.*only with method = \"kappa\"")
  d2 <- d
  d2$nvstat <- d2$nvstat + 1   # treatment in {1, 2}: not binary 0/1
  expect_error(
    drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d2,
           method = "kappa", tmodel = "linear"),
    "binary")
})
