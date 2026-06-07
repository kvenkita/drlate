# Kappa-weighting estimators (Sloczynski, Uysal & Wooldridge 2025, JBES;
# Stata kappalate): tau_a, tau_a,0, tau_a,10.

test_that("kappa (tau_a) equals its closed form; system square, moments zero", {
  d <- drlate_sim
  fit <- expect_valid_system(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, d,
                             method = "kappa")
  ps <- fitted(glm(rsncode ~ age + educ, binomial, data = d))
  z <- d$rsncode; dd <- d$nvstat; y <- d$lwage
  delta <- mean(z * y / ps - (1 - z) * y / (1 - ps))
  gam <- mean(1 - dd * (1 - z) / (1 - ps) - (1 - dd) * z / ps)
  expect_equal(unname(coef(fit)), c(delta / gam, delta, gam),
               tolerance = 1e-8)
})

test_that("kappa works with ivmodel = cbps", {
  fit <- drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ,
                data = drlate_sim, method = "kappa", ivmodel = "cbps")
  expect_true(all(is.finite(coef(fit))))
  expect_true(all(is.finite(sqrt(diag(fit$vcov3)))))
})

test_that("kappa0 (tau_a,0) equals its closed form; system square, moments zero", {
  d <- drlate_sim
  fit <- expect_valid_system(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, d,
                             method = "kappa0")
  ps <- fitted(glm(rsncode ~ age + educ, binomial, data = d))
  z <- d$rsncode; dd <- d$nvstat; y <- d$lwage
  delta <- mean(z * y / ps - (1 - z) * y / (1 - ps))
  gam0 <- mean((dd - 1) * (z / ps - (1 - z) / (1 - ps)))
  expect_equal(unname(coef(fit)), c(delta / gam0, delta, gam0),
               tolerance = 1e-8)
})

test_that("kappa methods validate inputs", {
  d <- drlate_sim
  expect_error(
    drlate(lwage ~ age, nvstat ~ 1, rsncode ~ age + educ, data = d,
           method = "kappa"),
    "covariates are not allowed")
  expect_error(
    drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d,
           method = "kappa", estimand = "latt"),
    "estimand = \"late\" only")
  expect_error(
    drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d,
           method = "kappa", ivmodel = "ipt"),
    "ipt.*not available")
  expect_error(
    drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d,
           method = "kappa0", estimand = "latt"),
    "estimand = \"late\" only")
  expect_error(
    drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d,
           method = "kappa10", ivmodel = "cbps"),
    "cbps.*only with method = \"kappa\"")
  expect_error(
    drlate(lwage ~ age, nvstat ~ 1, rsncode ~ age + educ, data = d,
           method = "kappa0"),
    "covariates are not allowed")
  d2 <- d
  d2$nvstat <- d2$nvstat + 1   # treatment in {1, 2}: not binary 0/1
  expect_error(
    drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d2,
           method = "kappa", tmodel = "linear"),
    "binary")
})

test_that("kappa10 (tau_a,10) equals its closed form; one reported coef", {
  d <- drlate_sim
  fit <- expect_valid_system(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, d,
                             method = "kappa10")
  ps <- fitted(glm(rsncode ~ age + educ, binomial, data = d))
  z <- d$rsncode; dd <- d$nvstat; y <- d$lwage
  kap1 <- z / ps - (1 - z) / (1 - ps)
  tau <- mean(dd * kap1 * y) / mean(dd * kap1) -
         mean((dd - 1) * kap1 * y) / mean((dd - 1) * kap1)
  expect_length(coef(fit), 1L)
  expect_named(coef(fit), "LATE: D on Y")
  expect_equal(unname(coef(fit)), tau, tolerance = 1e-8)
  expect_equal(dim(fit$vcov3), c(1L, 1L))
  expect_true(is.finite(sqrt(fit$vcov3[1, 1])))
})
