# Kappa-weighting estimators (Sloczynski, Uysal & Wooldridge 2025, JBES;
# Stata kappalate): tau_a, tau_a,0, tau_a,10.

test_that("kappa (tau_a) equals its closed form; system square, moments zero", {
  d <- drlate_sim
  fit <- expect_valid_system(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, d,
                             method = "kappa")
  ps <- fitted(glm(rsncode ~ age + educ, binomial, data = d))
  z <- d$rsncode; dd <- d$nvstat; y <- d$lwage
  delta <- mean(z * y / ps - (1 - z) * y / (1 - ps))
  gam <- mean(1 - dd * (1 - z) / (1 - ps) - (1 - dd) * z / ps)
  expect_equal(unname(coef(fit)), c(delta / gam, delta, gam),
               tolerance = 1e-8)
})

test_that("kappa works with ivmodel = cbps", {
  fit <- drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ,
                data = drlate_sim, method = "kappa", ivmodel = "cbps")
  expect_true(all(is.finite(coef(fit))))
  expect_true(all(is.finite(sqrt(diag(fit$vcov3)))))
})

test_that("kappa0 (tau_a,0) equals its closed form; system square, moments zero", {
  d <- drlate_sim
  fit <- expect_valid_system(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, d,
                             method = "kappa0")
  ps <- fitted(glm(rsncode ~ age + educ, binomial, data = d))
  z <- d$rsncode; dd <- d$nvstat; y <- d$lwage
  delta <- mean(z * y / ps - (1 - z) * y / (1 - ps))
  gam0 <- mean((dd - 1) * (z / ps - (1 - z) / (1 - ps)))
  expect_equal(unname(coef(fit)), c(delta / gam0, delta, gam0),
               tolerance = 1e-8)
})

test_that("kappa methods validate inputs", {
  d <- drlate_sim
  expect_error(
    drlate(lwage ~ age, nvstat ~ 1, rsncode ~ age + educ, data = d,
           method = "kappa"),
    "covariates are not allowed")
  expect_error(
    drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d,
           method = "kappa", estimand = "latt"),
    "estimand = \"late\" only")
  expect_error(
    drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d,
           method = "kappa", ivmodel = "ipt"),
    "ipt.*not available")
  expect_error(
    drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d,
           method = "kappa0", estimand = "latt"),
    "estimand = \"late\" only")
  expect_error(
    drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d,
           method = "kappa10", ivmodel = "cbps"),
    "cbps.*only with method = \"kappa\"")
  expect_error(
    drlate(lwage ~ age, nvstat ~ 1, rsncode ~ age + educ, data = d,
           method = "kappa0"),
    "covariates are not allowed")
  d2 <- d
  d2$nvstat <- d2$nvstat + 1   # treatment in {1, 2}: not binary 0/1
  expect_error(
    drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d2,
           method = "kappa", tmodel = "linear"),
    "binary")
})

test_that("kappa10 (tau_a,10) equals its closed form; one reported coef", {
  d <- drlate_sim
  fit <- expect_valid_system(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, d,
                             method = "kappa10")
  ps <- fitted(glm(rsncode ~ age + educ, binomial, data = d))
  z <- d$rsncode; dd <- d$nvstat; y <- d$lwage
  kap1 <- z / ps - (1 - z) / (1 - ps)
  tau <- mean(dd * kap1 * y) / mean(dd * kap1) -
         mean((dd - 1) * kap1 * y) / mean((dd - 1) * kap1)
  expect_length(coef(fit), 1L)
  expect_named(coef(fit), "LATE: D on Y")
  expect_equal(unname(coef(fit)), tau, tolerance = 1e-8)
  expect_equal(dim(fit$vcov3), c(1L, 1L))
  expect_true(is.finite(sqrt(fit$vcov3[1, 1])))
})

test_that("Fieller works for kappa/kappa0 and errors for kappa10", {
  d <- drlate_sim
  for (m in c("kappa", "kappa0")) {
    f <- drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d,
                method = m)
    ci <- confint(f, method = "fieller")
    expect_s3_class(ci, "drlate_fieller")
    expect_identical(ci$type, "bounded")
    expect_true(ci$lower < coef(f)[[1]] && ci$upper > coef(f)[[1]])
  }
  f10 <- drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d,
                method = "kappa10")
  expect_error(confint(f10, method = "fieller"), "kappa10")
})

test_that("kappa10 prints (first-stage z from gamma1) without error", {
  f10 <- drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ,
                data = drlate_sim, method = "kappa10")
  expect_output(print(f10), "First stage")
  expect_true(is.finite(firststage_z(f10)))
})

test_that("bootstrap works for kappa10 (single reported coefficient)", {
  fit <- drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ,
                data = drlate_sim, method = "kappa10",
                vcov = "bootstrap", boot_reps = 25L, boot_seed = 42)
  expect_identical(colnames(fit$boot$draws), "late")
  expect_true(is.finite(fit$boot$se[["late"]]))
  ci <- confint(fit)
  expect_identical(nrow(ci), 1L)
})

test_that("print shows the kappalate estimator aliases", {
  d <- drlate_sim
  spec <- function(...) drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ,
                               data = d, ...)
  expect_output(print(spec(method = "kappa")), "tau_a;", fixed = TRUE)
  expect_output(print(spec(method = "kappa0")), "tau_a,0", fixed = TRUE)
  expect_output(print(spec(method = "kappa10")), "tau_a,10", fixed = TRUE)
  expect_output(print(spec(method = "kappa")), "none (kappa weighting)",
                fixed = TRUE)
  expect_output(print(spec(method = "ipw")), "tau_u", fixed = TRUE)
  expect_output(print(spec(method = "ipw", normalized = FALSE)),
                "tau_a,1", fixed = TRUE)
})

test_that("drlate_compare accepts the kappa methods", {
  cmp <- drlate_compare(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ,
                        data = drlate_sim,
                        methods = c("ipw", "kappa", "kappa0", "kappa10"))
  expect_identical(nrow(cmp), 4L)
  expect_true(all(is.finite(cmp$estimate)))
  expect_true(all(is.finite(cmp$se)))
})

test_that("kappa methods accept sampling weights and cluster", {
  d <- drlate_sim
  d$wt <- 1 + (seq_len(nrow(d)) %% 3) / 10
  d$cl <- rep_len(1:50, nrow(d))
  for (m in c("kappa", "kappa0", "kappa10")) {
    f <- drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d,
                method = m, weights = "wt", cluster = "cl")
    expect_true(all(is.finite(coef(f))))
    expect_true(all(is.finite(sqrt(diag(f$vcov3)))))
    expect_identical(f$N_clust, 50L)
    # weights actually move the estimate
    f0 <- drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ, data = d,
                 method = m)
    expect_false(isTRUE(all.equal(coef(f)[[1]], coef(f0)[[1]])))
  }
})
