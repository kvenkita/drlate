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

test_that("balance(detail = TRUE) adds variance ratios and weighted arm means", {
  d <- drlate_sim
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age + educ, data = d)

  # Default call is unchanged (backward compatible)
  b0 <- balance(fit)
  expect_named(b0, c("variable", "smd_unweighted", "smd_weighted"))

  b <- balance(fit, detail = TRUE)
  expect_true(all(c("smd_unweighted", "smd_weighted",
                    "mean_weighted_z1", "mean_weighted_z0",
                    "vratio_unweighted", "vratio_weighted") %in% names(b)))

  # Hand checks for age on the original scale (complete data -> full sample)
  z <- fit$ctx$z; w <- fit$ctx$w; ps <- fit$ps$ps
  age <- d$age
  vr_u <- var(age[z == 1]) / var(age[z == 0])
  expect_equal(b$vratio_unweighted[b$variable == "age"], vr_u,
               tolerance = 1e-8)
  mw1 <- weighted.mean(age[z == 1], (w / ps)[z == 1])
  expect_equal(b$mean_weighted_z1[b$variable == "age"], mw1,
               tolerance = 1e-8)
})

test_that("plot.drlate supports density overlap and balance-density types", {
  skip_if_not_installed("ggplot2")
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age + educ,
                data = drlate_sim)
  has_density <- function(p) {
    any(vapply(p$layers, function(l) inherits(l$geom, "GeomDensity"),
               logical(1)))
  }
  expect_false(has_density(plot(fit, type = "overlap")))
  p_dens <- plot(fit, type = "overlap", geom = "density")
  expect_s3_class(p_dens, "ggplot")
  expect_true(has_density(p_dens))

  expect_s3_class(plot(fit, type = "balance_density"), "ggplot")
  expect_s3_class(plot(fit, type = "balance_density", var = "age"), "ggplot")
})
