# Complier profiling: kappa_weights() and complier_means().

test_that("kappa_weights() returns Abadie kappa matching its formula", {
  d <- drlate_sim
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age + educ, data = d)
  k <- kappa_weights(fit, normalize = FALSE)
  z <- fit$ctx$z; dd <- fit$ctx$d; ps <- fit$ps$ps
  kap <- 1 - dd * (1 - z) / (1 - ps) - (1 - dd) * z / ps
  expect_equal(unname(k), unname(kap), tolerance = 1e-12)
  expect_length(k, fit$ctx$n)
})

test_that("kappa_weights(normalize = TRUE) sums to one", {
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age + educ,
                data = drlate_sim)
  k <- kappa_weights(fit)
  expect_equal(sum(k), 1, tolerance = 1e-12)
})

test_that("kappa_weights() errors for RA and lean fits", {
  d <- drlate_sim
  fit_ra <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ 1, data = d,
                   method = "ra")
  expect_error(kappa_weights(fit_ra), "method = \"ra\"")
  fit_lean <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                     keep_data = FALSE)
  expect_error(kappa_weights(fit_lean), "keep_data")
})

test_that("complier_means() reports population vs complier means (original scale)", {
  d <- drlate_sim
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age + educ, data = d)
  cm <- complier_means(fit)
  expect_s3_class(cm, "data.frame")
  expect_named(cm, c("variable", "population_mean", "complier_mean", "difference"))
  # educ is a factor, so it enters as indicator columns (educcollege, ...)
  expect_true(all(c("age", "educcollege") %in% cm$variable))
  expect_false("(Intercept)" %in% cm$variable)

  # Hand computation on the ORIGINAL covariate scale (drlate_sim has no
  # missing data, so the estimation sample is the full data in order).
  kap <- kappa_weights(fit, normalize = FALSE)
  w <- fit$ctx$w
  pop  <- sum(w * d$age) / sum(w)
  comp <- sum(w * kap * d$age) / sum(w * kap)
  row <- cm[cm$variable == "age", ]
  expect_equal(row$population_mean, pop,  tolerance = 1e-8)
  expect_equal(row$complier_mean,   comp, tolerance = 1e-8)
  expect_equal(row$difference,      comp - pop, tolerance = 1e-8)
})

test_that("complier_means() accepts a covariate subset and errors for RA", {
  d <- drlate_sim
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age + educ, data = d)
  cm <- complier_means(fit, vars = "educcollege")
  expect_equal(cm$variable, "educcollege")
  fit_ra <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ 1, data = d,
                   method = "ra")
  expect_error(complier_means(fit_ra), "method = \"ra\"")
})
