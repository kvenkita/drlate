# Probit instrument propensity score (kappalate zmodel(probit)) and the
# drlate_compare normalization metadata.

test_that("probit IPS is restricted to the kappalate-covered methods", {
  d <- drlate_sim
  expect_error(
    drlate(lwage ~ age, nvstat ~ age, rsncode ~ age + educ, data = d,
           method = "ipwra", ivmodel = "probit"),
    "probit")
  expect_error(
    drlate(lwage ~ age + educ, nvstat ~ age + educ, rsncode ~ 1,
           data = d, method = "ra", ivmodel = "probit"),
    "probit")
  expect_error(
    drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d,
           method = "ipw", ivmodel = "probit", estimand = "latt"),
    "probit")
})

test_that("kappa + probit equals its closed form; system square, moments zero", {
  d <- drlate_sim
  fit <- expect_valid_system(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, d,
                             method = "kappa", ivmodel = "probit")
  ps <- fitted(glm(rsncode ~ age + educ, binomial("probit"), data = d,
                   control = glm.control(epsilon = 1e-12, maxit = 100)))
  z <- d$rsncode; dd <- d$nvstat; y <- d$lwage
  delta <- mean(z * y / ps - (1 - z) * y / (1 - ps))
  gam <- mean(1 - dd * (1 - z) / (1 - ps) - (1 - dd) * z / ps)
  expect_equal(unname(coef(fit)), c(delta / gam, delta, gam),
               tolerance = 1e-8)
})

test_that("ipw + probit equals the Hajek and raw closed forms", {
  d <- drlate_sim
  ps <- fitted(glm(rsncode ~ age + educ, binomial("probit"), data = d,
                   control = glm.control(epsilon = 1e-12, maxit = 100)))
  z <- d$rsncode; dd <- d$nvstat; y <- d$lwage

  # normalized (tau_u): Hajek means per instrument arm
  fit_n <- expect_valid_system(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ,
                               d, method = "ipw", ivmodel = "probit")
  hj <- function(v, arm) {
    wts <- if (arm == 1) z / ps else (1 - z) / (1 - ps)
    weighted.mean(v[z == arm], wts[z == arm])
  }
  num <- hj(y, 1) - hj(y, 0)
  den <- hj(dd, 1) - hj(dd, 0)
  expect_equal(unname(coef(fit_n)), c(num / den, num, den),
               tolerance = 1e-8)

  # unnormalized (tau_a,1): raw IPW means
  fit_u <- expect_valid_system(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ,
                               d, method = "ipw", ivmodel = "probit",
                               normalized = FALSE)
  num_u <- mean(z * y / ps - (1 - z) * y / (1 - ps))
  den_u <- mean(z * dd / ps - (1 - z) * dd / (1 - ps))
  expect_equal(unname(coef(fit_u)), c(num_u / den_u, num_u, den_u),
               tolerance = 1e-8)
})

test_that("kappa0 and kappa10 work with probit and print its name", {
  d <- drlate_sim
  f0 <- expect_valid_system(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, d,
                            method = "kappa0", ivmodel = "probit")
  expect_true(all(is.finite(sqrt(diag(f0$vcov3)))))
  f10 <- expect_valid_system(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, d,
                             method = "kappa10", ivmodel = "probit")
  expect_length(coef(f10), 1L)
  expect_true(is.finite(sqrt(f10$vcov3[1, 1])))
  expect_output(print(f10), "probit (MLE)", fixed = TRUE)
})

test_that("probit estimators match Stata kappalate (probit, SIPP)", {
  skip_if_no_fixture("kappalate_probit_all")
  d <- sipp_data()
  fx <- read_fixture("kappalate_probit_all")
  spec <- function(...) drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age_5,
                               data = d, ivmodel = "probit", ...)
  # tol_b = 1e-5 (vs 1e-6 for logit): the probit ML optimum is found by
  # different algorithms (IRLS at epsilon 1e-12 here, Newton-Raphson with
  # Stata's nrtolerance there), and Stata's stopping rule leaves ~3e-7
  # absolute slack in these weight-sensitive unnormalized estimators
  # (~5e-6 relative at SIPP's tau_a ~ 0.058; ~1.5e-6 of a standard error).
  # SEs still agree at 1e-4.
  expect_kappa_fixture(spec(method = "kappa"), fx, 1, tol_b = 1e-5)
  expect_kappa_fixture(spec(method = "ipw", normalized = FALSE), fx, 2,
                       tol_b = 1e-5)
  expect_kappa_fixture(spec(method = "kappa0"), fx, 3, tol_b = 1e-5)
  expect_kappa_fixture(spec(method = "kappa10"), fx, 4, tol_b = 1e-5)
  expect_kappa_fixture(spec(method = "ipw"), fx, 5, tol_b = 1e-5)
})

test_that("drlate_compare reports the estimator's own normalization", {
  cmp <- drlate_compare(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ,
                        data = drlate_sim,
                        methods = c("ipw", "kappa", "kappa0", "kappa10"))
  norm <- setNames(cmp$normalized, cmp$method)
  expect_true(norm[["ipw"]])        # tau_u, normalized
  expect_false(norm[["kappa"]])     # tau_a, unnormalized
  expect_false(norm[["kappa0"]])    # tau_a,0, unnormalized
  expect_true(norm[["kappa10"]])    # tau_a,10, normalized
})
