# Fieller confidence sets for the LATE ratio.

test_that("fieller quadratic solution matches polyroot", {
  num <- 0.3; denom <- 0.6
  V2 <- matrix(c(0.004, 0.001, 0.001, 0.002), 2)
  f <- fieller_ci(num, denom, V2, level = 0.95)
  q <- qchisq(0.95, 1)
  a <- denom^2 - q * V2[2, 2]
  b <- -2 * (num * denom - q * V2[1, 2])
  cc <- num^2 - q * V2[1, 1]
  roots <- sort(Re(polyroot(c(cc, b, a))))
  expect_equal(c(f$lower, f$upper), roots, tolerance = 1e-10)
  expect_equal(f$type, "bounded")
})

test_that("strong first stage: Fieller is close to (wider than) Wald", {
  fit <- drlate(lwage ~ age + educ, nvstat ~ age + educ,
                rsncode ~ age + educ, data = drlate_sim)
  f <- confint(fit, method = "fieller")
  expect_s3_class(f, "drlate_fieller")
  expect_equal(f$type, "bounded")
  wald <- confint(fit)["LATE: D on Y", ]
  # Same neighborhood (note: Fieller uses the joint (num, denom) covariance
  # while the Wald row uses the delta-method variance of the ratio)
  expect_lt(abs(f$lower - wald[1]), 0.05)
  expect_lt(abs(f$upper - wald[2]), 0.05)
  # The point estimate is inside the set
  expect_gt(coef(fit)[1], f$lower)
  expect_lt(coef(fit)[1], f$upper)
})

test_that("weak first stage yields complement or whole-line sets", {
  num <- 0.3; denom <- 0.05
  V2 <- matrix(c(0.01, 0, 0, 0.01), 2)   # denominator z = 0.5
  f <- fieller_ci(num, denom, V2)
  expect_true(f$type %in% c("complement", "whole-line"))
  out <- paste(capture.output(
    print(structure(c(f, estimand = "LATE"), class = "drlate_fieller"))),
    collapse = "\n")
  expect_match(out, "Inf")
})

test_that("weak-instrument print path shows the Fieller set", {
  set.seed(42)
  d <- drlate_sim[1:300, ]
  d$zweak <- sample(d$rsncode)
  fit <- drlate(lwage ~ age, nvstat ~ age, zweak ~ age, data = d)
  expect_lt(abs(firststage_z(fit)), 2)
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "Fieller")
})
