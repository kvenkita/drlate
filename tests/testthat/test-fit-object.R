# Fit-object internals retained for diagnostics/bootstrap/Fieller.

test_that("keep_data retains ctx, ps, and layout; FALSE drops them", {
  d <- drlate_sim
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d)
  expect_false(is.null(fit$ctx))
  expect_false(is.null(fit$ps))
  expect_false(is.null(fit$layout))
  expect_true(all(c("late", "num", "denom") %in% names(fit$layout)))
  expect_equal(length(fit$ps$ps), fit$N)

  lean <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                 keep_data = FALSE)
  expect_null(lean$ctx)
  expect_null(lean$ps)
  expect_false(is.null(lean$layout))  # layout is tiny; always kept
})

test_that("vcov_full has nonzero num-denom covariance (vcov3 is diagonal)", {
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = drlate_sim)
  idx <- c(fit$layout$num, fit$layout$denom)
  V2 <- fit$vcov_full[idx, idx]
  expect_true(abs(V2[1, 2]) > 0)
  expect_equal(fit$vcov3[1, 2], 0)
  # diagonal of vcov_full at the reported indices equals vcov3 diagonal
  expect_equal(unname(V2[1, 1]), unname(fit$vcov3[2, 2]))
  expect_equal(unname(V2[2, 2]), unname(fit$vcov3[3, 3]))
})

test_that("ra fits carry no ps", {
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ 1, data = drlate_sim,
                method = "ra")
  expect_null(fit$ps)
  expect_false(is.null(fit$ctx))
})
