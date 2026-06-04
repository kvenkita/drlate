# Balance diagnostics and plots.

test_that("balance() returns SMDs matching a hand computation", {
  d <- drlate_sim
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age + educ, data = d)
  b <- balance(fit)
  expect_s3_class(b, "data.frame")
  expect_named(b, c("variable", "smd_unweighted", "smd_weighted"))
  expect_true("age" %in% b$variable)
  expect_false("(Intercept)" %in% b$variable)

  # Hand computation for age (standardized by the fit, so recompute on
  # the standardized column retained in ctx)
  X <- fit$ctx$Xz
  z <- fit$ctx$z
  ps <- fit$ps$ps
  agec <- X[, "age"]
  pooled <- sqrt((var(agec[z == 1]) + var(agec[z == 0])) / 2)
  smd_u <- (mean(agec[z == 1]) - mean(agec[z == 0])) / pooled
  smd_w <- (weighted.mean(agec[z == 1], 1 / ps[z == 1]) -
            weighted.mean(agec[z == 0], 1 / (1 - ps[z == 0]))) / pooled
  expect_equal(b$smd_unweighted[b$variable == "age"], smd_u,
               tolerance = 1e-10)
  expect_equal(b$smd_weighted[b$variable == "age"], smd_w,
               tolerance = 1e-10)

  # Weighting should improve balance markedly in this DGP
  expect_lt(max(abs(b$smd_weighted)), max(abs(b$smd_unweighted)))
})

test_that("balance() errors for RA and keep_data = FALSE fits", {
  d <- drlate_sim
  fit_ra <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ 1, data = d,
                   method = "ra")
  expect_error(balance(fit_ra), "method = \"ra\"")
  fit_lean <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                     keep_data = FALSE)
  expect_error(balance(fit_lean), "keep_data")
})

test_that("plot.drlate returns ggplot objects for each type", {
  skip_if_not_installed("ggplot2")
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age + educ,
                data = drlate_sim)
  for (tp in c("overlap", "balance", "weights")) {
    p <- plot(fit, type = tp)
    expect_s3_class(p, "ggplot")
  }
  # LATT branch of the weights plot
  fit_latt <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age,
                     data = drlate_sim, estimand = "latt")
  expect_s3_class(plot(fit_latt, type = "weights"), "ggplot")
})

test_that("plot.drlate errors informatively for RA fits", {
  skip_if_not_installed("ggplot2")
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ 1, data = drlate_sim,
                method = "ra")
  expect_error(plot(fit, type = "overlap"), "method = \"ra\"")
  expect_error(plot(fit, type = "weights"), "method = \"ra\"")
})

test_that("print shows the first-stage z and flags weakness", {
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age,
                data = drlate_sim)
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "First stage \\(Z on D\\): z =")
  expect_no_match(out, "weak")  # strong instrument in the sim

  # Construct a weak-instrument fit: shuffle the instrument so the first
  # stage is ~0 on a small subsample
  set.seed(42)
  d <- drlate_sim[1:300, ]
  d$zweak <- sample(d$rsncode)
  fit_w <- drlate(lwage ~ age, nvstat ~ age, zweak ~ age, data = d)
  out_w <- paste(capture.output(print(fit_w)), collapse = "\n")
  expect_match(out_w, "weak")
})
