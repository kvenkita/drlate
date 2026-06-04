# LATE x {RA, IPW, AIPW}: internal consistency + golden fixtures.

# Helper: assert the stacked system is square with zero averaged moments
expect_valid_system <- function(formY, formD, formZ, data, ...) {
  args <- list(...)
  fit <- drlate(formY, formD, formZ, data = data, ...)
  fam <- c(linear = "gaussian", logit = "binomial", poisson = "poisson")
  ctx <- build_ctx(formY, formD, formZ, data,
                   omodel = fam[[args$omodel %||% "linear"]],
                   tmodel = fam[[args$tmodel %||% "logit"]],
                   ivmodel = args$ivmodel %||% "logit",
                   method = args$method %||% "ipwra",
                   estimand = "late",
                   normalized = args$normalized %||% TRUE)
  ps <- if (ctx$method != "ra") fit_ps(ctx) else NULL
  est <- estimate_late(ctx, ps)
  sys <- assemble_system(est$blocks)
  gbar <- colMeans(ctx$w * sys$g(sys$theta0))
  expect_lt(max(abs(gbar)), 1e-7)
  expect_equal(length(sys$theta0), ncol(sys$g(sys$theta0)))
  fit
}

`%||%` <- function(a, b) if (is.null(a)) b else a

test_that("RA system is square with zero moments; estimate matches closed form", {
  d <- drlate_sim
  fit <- expect_valid_system(lwage ~ age + educ, nvstat ~ age + educ,
                             rsncode ~ 1, d, method = "ra")
  # Closed form: plain regressions per arm, average predictions
  m <- model.matrix(~ age + educ, d)
  s1 <- d$rsncode == 1
  fy1 <- lm(lwage ~ age + educ, d[s1, ])
  fy0 <- lm(lwage ~ age + educ, d[!s1, ])
  fd1 <- glm(nvstat ~ age + educ, binomial, d[s1, ])
  fd0 <- glm(nvstat ~ age + educ, binomial, d[!s1, ])
  num <- mean(m %*% coef(fy1)) - mean(m %*% coef(fy0))
  den <- mean(plogis(m %*% coef(fd1))) - mean(plogis(m %*% coef(fd0)))
  expect_equal(unname(coef(fit)), c(num / den, num, den), tolerance = 1e-8)
})

test_that("IPW normalized equals Hajek-weighted means", {
  d <- drlate_sim
  fit <- expect_valid_system(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, d,
                             method = "ipw")
  ps <- fitted(glm(rsncode ~ age + educ, binomial, data = d))
  z <- d$rsncode
  hj <- function(v, arm) {
    wts <- if (arm == 1) (z / ps) else ((1 - z) / (1 - ps))
    sel <- z == arm
    weighted.mean(v[sel], wts[sel])
  }
  num <- hj(d$lwage, 1) - hj(d$lwage, 0)
  den <- hj(d$nvstat, 1) - hj(d$nvstat, 0)
  expect_equal(unname(coef(fit)), c(num / den, num, den), tolerance = 1e-8)
})

test_that("IPW unnormalized equals raw IPW means", {
  d <- drlate_sim
  fit <- expect_valid_system(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, d,
                             method = "ipw", normalized = FALSE)
  ps <- fitted(glm(rsncode ~ age + educ, binomial, data = d))
  z <- d$rsncode
  num <- mean(z / ps * d$lwage) - mean((1 - z) / (1 - ps) * d$lwage)
  den <- mean(z / ps * d$nvstat) - mean((1 - z) / (1 - ps) * d$nvstat)
  expect_equal(unname(coef(fit)), c(num / den, num, den), tolerance = 1e-8)
})

test_that("AIPW unnormalized equals the influence-function formula", {
  d <- drlate_sim
  fit <- expect_valid_system(lwage ~ age + educ, nvstat ~ age + educ,
                             rsncode ~ age + educ, d,
                             method = "aipw", normalized = FALSE)
  ps <- fitted(glm(rsncode ~ age + educ, binomial, data = d))
  z <- d$rsncode
  m <- model.matrix(~ age + educ, d)
  s1 <- z == 1
  mu_y1 <- drop(m %*% coef(lm(lwage ~ age + educ, d[s1, ])))
  mu_y0 <- drop(m %*% coef(lm(lwage ~ age + educ, d[!s1, ])))
  mu_d1 <- plogis(drop(m %*% coef(glm(nvstat ~ age + educ, binomial, d[s1, ]))))
  mu_d0 <- plogis(drop(m %*% coef(glm(nvstat ~ age + educ, binomial, d[!s1, ]))))
  num <- mean((z * d$lwage - (z - ps) * mu_y1) / ps) -
         mean(((1 - z) * d$lwage + (z - ps) * mu_y0) / (1 - ps))
  den <- mean((z * d$nvstat - (z - ps) * mu_d1) / ps) -
         mean(((1 - z) * d$nvstat + (z - ps) * mu_d0) / (1 - ps))
  expect_equal(unname(coef(fit)), c(num / den, num, den), tolerance = 1e-8)
})

test_that("AIPW normalized system is square with zero moments", {
  d <- drlate_sim
  fit <- expect_valid_system(lwage ~ age + educ, nvstat ~ age + educ,
                             rsncode ~ age + educ, d, method = "aipw")
  expect_s3_class(fit, "drlate")
  # Hajek-normalized AIPW closed form
  ps <- fitted(glm(rsncode ~ age + educ, binomial, data = d))
  z <- d$rsncode
  m <- model.matrix(~ age + educ, d)
  s1 <- z == 1
  mu_y1 <- drop(m %*% coef(lm(lwage ~ age + educ, d[s1, ])))
  mu_y0 <- drop(m %*% coef(lm(lwage ~ age + educ, d[!s1, ])))
  om1 <- (z / ps) / mean(z / ps)
  om0 <- ((1 - z) / (1 - ps)) / mean((1 - z) / (1 - ps))
  num <- mean(om1 * (d$lwage - mu_y1) + mu_y1) -
         mean(om0 * (d$lwage - mu_y0) + mu_y0)
  expect_equal(unname(coef(fit)[2]), num, tolerance = 1e-8)
})

test_that("all omodel/tmodel families run for each method", {
  d <- drlate_sim
  for (om in c("linear", "logit", "poisson")) {
    yf <- switch(om, linear = lwage ~ age, logit = hijob ~ age,
                 poisson = kwage ~ age)
    for (tm in c("logit", "linear", "poisson")) {
      for (me in c("ipwra", "aipw", "ra")) {
        zf <- if (me == "ra") rsncode ~ 1 else rsncode ~ age
        fit <- drlate(yf, nvstat ~ age, zf, data = d,
                      omodel = om, tmodel = tm, method = me)
        expect_true(is.finite(coef(fit)[1]),
                    label = paste(om, tm, me, "finite estimate"))
        expect_true(all(is.finite(sqrt(diag(fit$vcov3)))),
                    label = paste(om, tm, me, "finite SEs"))
      }
    }
  }
})

# ---------------------------------------------------------------------------
# Golden fixtures
# ---------------------------------------------------------------------------

golden <- function(id, ...) {
  test_that(paste0("matches Stata: ", id), {
    skip_if_no_fixture(id)
    fit <- drlate(..., data = sipp_data())
    expect_matches_fixture(fit, id)
  })
}

golden("late_ipw_nrm", lwage ~ 1, nvstat ~ 1, rsncode ~ age_5,
       method = "ipw")
golden("late_ipw_unnrm", lwage ~ 1, nvstat ~ 1, rsncode ~ age_5,
       method = "ipw", normalized = FALSE)
golden("late_aipw_nrm", lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5,
       method = "aipw")
golden("late_aipw_unnrm", lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5,
       method = "aipw", normalized = FALSE)
golden("late_aipw_nrm_logit_y", hiwage ~ age_5, nvstat ~ age_5,
       rsncode ~ age_5, method = "aipw", omodel = "logit")
golden("late_aipw_nrm_pois_y", kwage ~ age_5, nvstat ~ age_5,
       rsncode ~ age_5, method = "aipw", omodel = "poisson")
golden("late_ra", lwage ~ age_5, nvstat ~ age_5, rsncode ~ 1,
       method = "ra")
golden("late_ra_pois_y", kwage ~ age_5, nvstat ~ age_5, rsncode ~ 1,
       method = "ra", omodel = "poisson")
