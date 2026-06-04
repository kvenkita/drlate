# Multi-estimator comparison.

test_that("drlate_compare rows equal standalone fits", {
  d <- drlate_sim
  suppressMessages(
    cmp <- drlate_compare(lwage ~ age, nvstat ~ age, rsncode ~ age,
                          data = d)
  )
  expect_s3_class(cmp, "drlate_compare")
  expect_equal(cmp$method, c("ipwra", "ipw", "aipw", "ra"))

  fit_ipwra <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d)
  expect_equal(cmp$estimate[1], unname(coef(fit_ipwra)[1]))

  fit_ipw <- drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age, data = d,
                    method = "ipw")
  expect_equal(cmp$estimate[2], unname(coef(fit_ipw)[1]))

  fit_ra <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ 1, data = d,
                   method = "ra")
  expect_equal(cmp$estimate[4], unname(coef(fit_ra)[1]))
})

test_that("formula auto-adjustment messages are emitted", {
  d <- drlate_sim
  expect_message(
    drlate_compare(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                   methods = "ipw"),
    "dropping outcome/treatment covariates"
  )
  expect_message(
    drlate_compare(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                   methods = "ra"),
    "dropping instrument covariates"
  )
})

test_that("both_norms adds unnormalized ipw/aipw rows", {
  d <- drlate_sim
  suppressMessages(
    cmp <- drlate_compare(lwage ~ age, nvstat ~ age, rsncode ~ age,
                          data = d, methods = c("ipw", "aipw"),
                          both_norms = TRUE)
  )
  expect_equal(nrow(cmp), 4L)
  expect_equal(sum(!cmp$normalized), 2L)
})

test_that("print and plot work", {
  d <- drlate_sim
  suppressMessages(
    cmp <- drlate_compare(lwage ~ age, nvstat ~ age, rsncode ~ age,
                          data = d, methods = c("ipwra", "ra"))
  )
  out <- paste(capture.output(print(cmp)), collapse = "\n")
  expect_match(out, "Estimator comparison")
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(cmp), "ggplot")
})

test_that("estimand passes through (latt comparison)", {
  d <- drlate_sim
  suppressMessages(
    cmp <- drlate_compare(lwage ~ age, nvstat ~ age, rsncode ~ age,
                          data = d, methods = c("ipwra", "aipw"),
                          estimand = "latt")
  )
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                estimand = "latt")
  expect_equal(cmp$estimate[1], unname(coef(fit)[1]))
})
