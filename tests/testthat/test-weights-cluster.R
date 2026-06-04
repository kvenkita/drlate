# Sampling weights and clustered standard errors.

test_that("weighted IPWRA matches the closed-form weighted computation", {
  d <- drlate_sim
  d$wpw <- 1 + (d$kwage - floor(d$kwage))
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                weights = "wpw")
  wv <- d$wpw
  ps <- fitted(suppressWarnings(
    glm(rsncode ~ age, quasibinomial, data = d, weights = wv)))
  m <- model.matrix(~ age, d)
  s1 <- d$rsncode == 1
  da <- transform(d, w1 = wv / ps, w0 = wv / (1 - ps))
  fy1 <- lm(lwage ~ age, da[s1, ], weights = w1)
  fy0 <- lm(lwage ~ age, da[!s1, ], weights = w0)
  fd1 <- suppressWarnings(glm(nvstat ~ age, quasibinomial, da[s1, ],
                              weights = w1))
  fd0 <- suppressWarnings(glm(nvstat ~ age, quasibinomial, da[!s1, ],
                              weights = w0))
  num <- weighted.mean(m %*% coef(fy1), wv) -
         weighted.mean(m %*% coef(fy0), wv)
  den <- weighted.mean(plogis(m %*% coef(fd1)), wv) -
         weighted.mean(plogis(m %*% coef(fd0)), wv)
  expect_equal(unname(coef(fit)), c(num / den, num, den), tolerance = 1e-8)
})

test_that("weights as a vector and as a column name agree", {
  d <- drlate_sim
  d$wpw <- 1 + (d$kwage - floor(d$kwage))
  f1 <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
               weights = "wpw")
  f2 <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
               weights = d$wpw)
  expect_equal(coef(f1), coef(f2))
  expect_equal(f1$vcov3, f2$vcov3)
})

test_that("singleton clusters reproduce robust SEs up to M/(M-1)", {
  d <- drlate_sim
  fit_r <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d)
  fit_c <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                  cluster = seq_len(nrow(d)))
  n <- nrow(d)
  expect_equal(fit_c$N_clust, n)
  expect_equal(diag(fit_c$vcov3), diag(fit_r$vcov3) * n / (n - 1),
               tolerance = 1e-10)
})

test_that("clustered SEs run for every estimand/method combination", {
  d <- drlate_sim
  cl <- d$educ
  for (es in c("late", "latt")) {
    for (me in c("ipwra", "aipw", "ra", "ipw")) {
      yf <- if (me == "ipw") lwage ~ 1 else lwage ~ age
      df <- if (me == "ipw") nvstat ~ 1 else nvstat ~ age
      zf <- if (me == "ra") rsncode ~ 1 else rsncode ~ age
      fit <- drlate(yf, df, zf, data = d, estimand = es, method = me,
                    cluster = cl)
      expect_equal(fit$N_clust, 3L, label = paste(es, me, "N_clust"))
      expect_true(all(is.finite(sqrt(diag(fit$vcov3)))),
                  label = paste(es, me, "clustered SEs"))
    }
  }
})

test_that("weights and cluster combine", {
  d <- drlate_sim
  d$wpw <- 1 + (d$kwage - floor(d$kwage))
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
                weights = "wpw", cluster = d$educ)
  expect_true(all(is.finite(coef(fit))))
  expect_true(all(is.finite(sqrt(diag(fit$vcov3)))))
})

# ---------------------------------------------------------------------------
# Golden fixtures
# ---------------------------------------------------------------------------

test_that("matches Stata: late_ipwra_pw", {
  skip_if_no_fixture("late_ipwra_pw")
  d <- sipp_data()
  fit <- drlate(lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5, data = d,
                weights = "wpw")
  expect_matches_fixture(fit, "late_ipwra_pw")
})

test_that("matches Stata: late_ipwra_cluster", {
  skip_if_no_fixture("late_ipwra_cluster")
  d <- sipp_data()
  fit <- drlate(lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5, data = d,
                cluster = "cluvar")
  expect_matches_fixture(fit, "late_ipwra_cluster")
})

test_that("matches Stata: late_ipwra_pw_cluster", {
  skip_if_no_fixture("late_ipwra_pw_cluster")
  d <- sipp_data()
  fit <- drlate(lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5, data = d,
                weights = "wpw", cluster = "cluvar")
  expect_matches_fixture(fit, "late_ipwra_pw_cluster")
})

test_that("matches Stata: late_aipw_pw", {
  skip_if_no_fixture("late_aipw_pw")
  d <- sipp_data()
  fit <- drlate(lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5, data = d,
                method = "aipw", weights = "wpw")
  expect_matches_fixture(fit, "late_aipw_pw")
})

test_that("matches Stata: latt_ipwra_pw", {
  skip_if_no_fixture("latt_ipwra_pw")
  d <- sipp_data()
  fit <- drlate(lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5, data = d,
                estimand = "latt", weights = "wpw")
  expect_matches_fixture(fit, "latt_ipwra_pw")
})

test_that("matches Stata: latt_ipwra_cluster", {
  skip_if_no_fixture("latt_ipwra_cluster")
  d <- sipp_data()
  fit <- drlate(lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5, data = d,
                estimand = "latt", cluster = "cluvar")
  expect_matches_fixture(fit, "latt_ipwra_cluster")
})
