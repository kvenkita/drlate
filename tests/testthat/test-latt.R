# LATT estimators: internal consistency + golden fixtures.

test_that("LATT IPWRA matches the closed-form sequential computation", {
  d <- drlate_sim
  fit <- drlate(lwage ~ age + educ, nvstat ~ age + educ, rsncode ~ age + educ,
                data = d, estimand = "latt")
  ps <- fitted(glm(rsncode ~ age + educ, binomial, data = d))
  z <- d$rsncode
  s1 <- z == 1
  m <- model.matrix(~ age + educ, d)
  da <- transform(d, attw = ps / (1 - ps))
  y1c <- mean(d$lwage[s1])
  fy0 <- lm(lwage ~ age + educ, da[!s1, ], weights = attw)
  fd1 <- mean(d$nvstat[s1])
  fd0 <- suppressWarnings(glm(nvstat ~ age + educ, quasibinomial, da[!s1, ],
                              weights = attw))
  mu_y0 <- drop(m %*% coef(fy0))
  mu_d0 <- plogis(drop(m %*% coef(fd0)))
  num <- y1c - mean(mu_y0[s1])
  den <- fd1 - mean(mu_d0[s1])
  expect_equal(unname(coef(fit)), c(num / den, num, den), tolerance = 1e-8)
})

test_that("LATT IPW normalized equals odds-weighted Hajek means", {
  d <- drlate_sim
  fit <- drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d,
                estimand = "latt", method = "ipw")
  ps <- fitted(glm(rsncode ~ age + educ, binomial, data = d))
  z <- d$rsncode
  s1 <- z == 1
  odds <- ps / (1 - ps)
  num <- mean(d$lwage[s1]) - weighted.mean(d$lwage[!s1], odds[!s1])
  den <- mean(d$nvstat[s1]) - weighted.mean(d$nvstat[!s1], odds[!s1])
  expect_equal(unname(coef(fit)), c(num / den, num, den), tolerance = 1e-8)
})

test_that("LATT IPW unnormalized equals raw treated-share-scaled means", {
  d <- drlate_sim
  fit <- drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d,
                estimand = "latt", method = "ipw", normalized = FALSE)
  ps <- fitted(glm(rsncode ~ age + educ, binomial, data = d))
  z <- d$rsncode
  w1 <- mean(z)
  om1 <- z / w1
  om0 <- (1 - z) * ps / (1 - ps) / w1
  num <- mean(om1 * d$lwage) - mean(om0 * d$lwage)
  den <- mean(om1 * d$nvstat) - mean(om0 * d$nvstat)
  expect_equal(unname(coef(fit)), c(num / den, num, den), tolerance = 1e-8)
})

test_that("LATT AIPW (both normalizations) matches closed form", {
  d <- drlate_sim
  ps <- fitted(glm(rsncode ~ age + educ, binomial, data = d))
  z <- d$rsncode
  s1 <- z == 1
  m <- model.matrix(~ age + educ, d)
  mu_y0 <- drop(m %*% coef(lm(lwage ~ age + educ, d[!s1, ])))
  mu_d0 <- plogis(drop(m %*% coef(glm(nvstat ~ age + educ, binomial,
                                      d[!s1, ]))))
  attw <- (1 - z) * ps / (1 - ps)

  # Unnormalized
  fit_u <- drlate(lwage ~ age + educ, nvstat ~ age + educ,
                  rsncode ~ age + educ, data = d,
                  estimand = "latt", method = "aipw", normalized = FALSE)
  w1 <- mean(z)
  num_u <- mean(z / w1 * d$lwage) -
           (mean(attw * (d$lwage - mu_y0) / w1) + mean(z * mu_y0 / w1))
  den_u <- mean(z / w1 * d$nvstat) -
           (mean(attw * (d$nvstat - mu_d0) / w1) + mean(z * mu_d0 / w1))
  expect_equal(unname(coef(fit_u)), c(num_u / den_u, num_u, den_u),
               tolerance = 1e-8)

  # Normalized
  fit_n <- drlate(lwage ~ age + educ, nvstat ~ age + educ,
                  rsncode ~ age + educ, data = d,
                  estimand = "latt", method = "aipw")
  wn <- mean(attw)
  num_n <- mean(z / w1 * d$lwage) -
           (mean(attw * (d$lwage - mu_y0) / wn) + mean(z * mu_y0 / w1))
  den_n <- mean(z / w1 * d$nvstat) -
           (mean(attw * (d$nvstat - mu_d0) / wn) + mean(z * mu_d0 / w1))
  expect_equal(unname(coef(fit_n)), c(num_n / den_n, num_n, den_n),
               tolerance = 1e-8)
})

test_that("LATT RA matches closed form", {
  d <- drlate_sim
  fit <- drlate(lwage ~ age + educ, nvstat ~ age + educ, rsncode ~ 1,
                data = d, estimand = "latt", method = "ra")
  z <- d$rsncode
  s1 <- z == 1
  m <- model.matrix(~ age + educ, d)
  mu_y0 <- drop(m %*% coef(lm(lwage ~ age + educ, d[!s1, ])))
  mu_d0 <- plogis(drop(m %*% coef(glm(nvstat ~ age + educ, binomial,
                                      d[!s1, ]))))
  num <- mean(d$lwage[s1]) - mean(mu_y0[s1])
  den <- mean(d$nvstat[s1]) - mean(mu_d0[s1])
  expect_equal(unname(coef(fit)), c(num / den, num, den), tolerance = 1e-8)
})

test_that("all LATT methods run across families and produce finite SEs", {
  d <- drlate_sim
  grid <- list(
    list(method = "ipwra", zf = rsncode ~ age),
    list(method = "ra",    zf = rsncode ~ 1),
    list(method = "aipw",  zf = rsncode ~ age)
  )
  for (om in c("linear", "logit", "poisson")) {
    yf <- switch(om, linear = lwage ~ age, logit = hijob ~ age,
                 poisson = kwage ~ age)
    for (spec in grid) {
      fit <- drlate(yf, nvstat ~ age, spec$zf, data = d,
                    estimand = "latt", method = spec$method, omodel = om)
      expect_true(all(is.finite(coef(fit))),
                  label = paste("latt", om, spec$method, "estimates"))
      expect_true(all(is.finite(sqrt(diag(fit$vcov3)))),
                  label = paste("latt", om, spec$method, "SEs"))
    }
  }
})

test_that("LATT degenerate compliance cases run", {
  d0 <- drlate_sim; d0$nvstat[d0$rsncode == 0] <- 0L
  d1 <- drlate_sim; d1$nvstat[d1$rsncode == 1] <- 1L
  for (me in c("ipwra", "aipw", "ra", "ipw")) {
    yf <- if (me == "ipw") lwage ~ 1 else lwage ~ age
    df <- if (me == "ipw") nvstat ~ 1 else nvstat ~ age
    zf <- if (me == "ra") rsncode ~ 1 else rsncode ~ age
    f0 <- drlate(yf, df, zf, data = d0, estimand = "latt", method = me)
    expect_true(all(is.finite(coef(f0))), label = paste("latt z0deg", me))
    f1 <- drlate(yf, df, zf, data = d1, estimand = "latt", method = me)
    expect_true(all(is.finite(coef(f1))), label = paste("latt z1deg", me))
    expect_equal(f1$dmeanz1, 1)
  }
})

# ---------------------------------------------------------------------------
# Golden fixtures
# ---------------------------------------------------------------------------

golden_latt <- function(id, ...) {
  test_that(paste0("matches Stata: ", id), {
    skip_if_no_fixture(id)
    fit <- drlate(..., data = sipp_data(), estimand = "latt")
    expect_matches_fixture(fit, id)
  })
}

golden_latt("latt_ipwra", lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5)
golden_latt("latt_ipw_nrm", lwage ~ 1, nvstat ~ 1, rsncode ~ age_5,
            method = "ipw")
golden_latt("latt_ipw_unnrm", lwage ~ 1, nvstat ~ 1, rsncode ~ age_5,
            method = "ipw", normalized = FALSE)
golden_latt("latt_aipw_nrm", lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5,
            method = "aipw")
golden_latt("latt_aipw_unnrm", lwage ~ age_5, nvstat ~ age_5,
            rsncode ~ age_5, method = "aipw", normalized = FALSE)
golden_latt("latt_ra", lwage ~ age_5, nvstat ~ age_5, rsncode ~ 1,
            method = "ra")
