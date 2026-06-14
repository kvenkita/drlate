# Probit and fractional (flogit/fprobit) outcome and treatment families
# (parity with Stata lateffects omodel/tmodel options).

# Binary-outcome DGP (for probit and logit/probit-equivalence tests)
bin_data <- function(n = 1500, seed = 3) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.4 * x))
  d <- rbinom(n, 1, plogis(-0.2 + 1.2 * z + 0.3 * x))
  y <- rbinom(n, 1, plogis(0.1 + 0.7 * d + 0.4 * x))
  data.frame(y, d, z, x)
}

# Fractional-outcome DGP (outcome strictly in (0, 1))
frac_data <- function(n = 1500, seed = 5) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.4 * x))
  d <- rbinom(n, 1, plogis(-0.2 + 1.2 * z + 0.3 * x))
  y <- plogis(0.2 + 0.8 * d + 0.5 * x + rnorm(n))
  data.frame(y, d, z, x)
}

test_that("flogit equals logit and fprobit equals probit on a binary outcome", {
  dat <- bin_data()
  f_logit  <- drlate(y ~ x, d ~ x, z ~ x, data = dat, omodel = "logit")
  f_flogit <- drlate(y ~ x, d ~ x, z ~ x, data = dat, omodel = "flogit")
  expect_equal(coef(f_flogit), coef(f_logit), tolerance = 1e-12)
  expect_equal(vcov(f_flogit), vcov(f_logit), tolerance = 1e-12)

  f_probit  <- drlate(y ~ x, d ~ x, z ~ x, data = dat, omodel = "probit")
  f_fprobit <- drlate(y ~ x, d ~ x, z ~ x, data = dat, omodel = "fprobit")
  expect_equal(coef(f_fprobit), coef(f_probit), tolerance = 1e-12)
  expect_equal(vcov(f_fprobit), vcov(f_probit), tolerance = 1e-12)
})

test_that("the stacked system is valid (score root) for the new families", {
  bd <- bin_data(); fd <- frac_data()
  expect_valid_system(y ~ x, d ~ x, z ~ x, data = bd, omodel = "probit")
  expect_valid_system(y ~ x, d ~ x, z ~ x, data = bd, tmodel = "probit")
  expect_valid_system(y ~ x, d ~ x, z ~ x, data = fd, omodel = "flogit")
  expect_valid_system(y ~ x, d ~ x, z ~ x, data = fd, omodel = "fprobit")
})

test_that("probit outcome arm coefficients match a direct weighted probit glm", {
  dat <- bin_data()
  fit <- drlate(y ~ x, d ~ x, z ~ x, data = dat, omodel = "probit")
  Xo <- fit$ctx$Xo
  z <- fit$ctx$z
  w1 <- 1 / fit$ps$ps
  g1 <- suppressWarnings(stats::glm.fit(
    Xo[z == 1, ], dat$y[z == 1], weights = w1[z == 1],
    family = stats::quasibinomial("probit"),
    control = stats::glm.control(epsilon = 1e-12, maxit = 100L)))
  expect_equal(unname(fit$theta[fit$layout$y1]), unname(stats::coef(g1)),
               tolerance = 1e-7)
})

test_that("fractional outcomes are estimated and reported with their label", {
  fd <- frac_data()
  fit <- drlate(y ~ x, d ~ x, z ~ x, data = fd, omodel = "flogit")
  expect_true(is.finite(coef(fit)[[1]]))   # the LATE estimate
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "fractional logit")
})

test_that("outcome-domain validation fires for binary and fractional families", {
  d <- drlate_sim                         # lwage is continuous
  expect_error(drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                      omodel = "logit"), "binary")
  expect_error(drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                      omodel = "probit"), "binary")
  expect_error(drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                      omodel = "flogit"), "0, 1")
  # A fractional outcome out of range is rejected
  fd <- frac_data(); fd$y[1] <- 1.5
  expect_error(drlate(y ~ x, d ~ x, z ~ x, data = fd, omodel = "fprobit"),
               "0, 1")
})

test_that("probit IPWRA analytic SE agrees with the bootstrap", {
  dat <- bin_data(n = 1200)
  fit  <- drlate(y ~ x, d ~ x, z ~ x, data = dat, omodel = "probit")
  fitb <- drlate(y ~ x, d ~ x, z ~ x, data = dat, omodel = "probit",
                 vcov = "bootstrap", boot_reps = 299L, boot_seed = 1L)
  se_a <- sqrt(diag(vcov(fit)))[[1]]    # SE of the LATE
  se_b <- sqrt(diag(vcov(fitb)))[[1]]
  expect_equal(se_b, se_a, tolerance = 0.2)   # within 20%
})
