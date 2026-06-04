# Nonparametric bootstrap inference.

test_that("bootstrap SEs are finite, reproducible, and near analytic", {
  d <- drlate_sim
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                vcov = "bootstrap", boot_reps = 199, boot_seed = 7)
  expect_equal(fit$vcov_method, "bootstrap")
  expect_true(all(is.finite(fit$boot$se)))
  expect_lte(fit$boot$reps_ok, 199)

  # Same seed reproduces; analytic SE is in the same ballpark
  fit2 <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                 vcov = "bootstrap", boot_reps = 199, boot_seed = 7)
  expect_equal(fit$boot$se, fit2$boot$se)
  ratio <- fit$boot$se[1] / sqrt(fit$vcov3[1, 1])
  expect_gt(ratio, 0.7)
  expect_lt(ratio, 1.4)
})

test_that("bootstrap point estimates are unchanged from analytic fit", {
  d <- drlate_sim
  fa <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d)
  fb <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
               vcov = "bootstrap", boot_reps = 50, boot_seed = 1)
  expect_equal(coef(fa), coef(fb))
  expect_equal(fa$vcov3, fb$vcov3)  # analytic V still stored
})

test_that("cluster bootstrap resamples whole clusters", {
  d <- drlate_sim
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                cluster = d$educ, vcov = "bootstrap", boot_reps = 25,
                boot_seed = 3)
  # With only 3 clusters, draws are coarse but must be finite often enough
  expect_gt(fit$boot$reps_ok, 0)
  expect_true(all(is.finite(fit$boot$se)))
})

test_that("failed draws are skipped and counted", {
  # Nearly-empty Z=0 arm: a noticeable share of resamples lose the arm
  # entirely; those draws must be skipped (NA), counted, and warned about,
  # not crash the bootstrap.
  d1 <- drlate_sim[drlate_sim$rsncode == 1, ][1:38, ]
  d0 <- drlate_sim[drlate_sim$rsncode == 0, ][1:2, ]
  d <- rbind(d1, d0)

  warnings_seen <- character()
  fit <- withCallingHandlers(
    drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
           vcov = "bootstrap", boot_reps = 99, boot_seed = 11),
    warning = function(cnd) {
      warnings_seen <<- c(warnings_seen, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("draws failed", warnings_seen)))
  expect_lt(fit$boot$reps_ok, 99)
})

test_that("print/confint reflect bootstrap inference", {
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age,
                data = drlate_sim, vcov = "bootstrap", boot_reps = 99,
                boot_seed = 5)
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "nonparametric bootstrap")
  ci <- confint(fit)
  expect_equal(dim(ci), c(3L, 2L))
  # Percentile CI must come from draw quantiles
  expect_equal(unname(ci[1, 1]),
               unname(quantile(fit$boot$draws[, "late"], 0.025)))
})
