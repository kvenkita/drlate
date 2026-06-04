# One-sided noncompliance: degenerate compliance means drop the
# corresponding treatment regression from the stacked system
# (drlate_estimate_late.ado per-case gmm variants).

# Construct one-sided noncompliance variants of the simulated data
onesided_z0 <- function() {
  d <- drlate_sim
  d$nvstat[d$rsncode == 0] <- 0L   # no always-takers: D = 0 when Z = 0
  d
}
onesided_z1 <- function() {
  d <- drlate_sim
  d$nvstat[d$rsncode == 1] <- 1L   # no never-takers: D = 1 when Z = 1
  d
}

methods_grid <- list(
  list(method = "ipwra", yf = lwage ~ age, df = nvstat ~ age,
       zf = rsncode ~ age),
  list(method = "ra",    yf = lwage ~ age, df = nvstat ~ age,
       zf = rsncode ~ 1),
  list(method = "ipw",   yf = lwage ~ 1,   df = nvstat ~ 1,
       zf = rsncode ~ age),
  list(method = "aipw",  yf = lwage ~ age, df = nvstat ~ age,
       zf = rsncode ~ age)
)

for (deg in c("z0deg", "z1deg")) {
  dat <- if (deg == "z0deg") onesided_z0() else onesided_z1()
  for (spec in methods_grid) {
    for (nrm in if (spec$method %in% c("ipw", "aipw")) c(TRUE, FALSE)
                else TRUE) {
      label <- paste(deg, spec$method, if (nrm) "nrm" else "unnrm")
      test_that(paste("degenerate case runs:", label), {
        fit <- drlate(spec$yf, spec$df, spec$zf, data = dat,
                      method = spec$method, normalized = nrm)
        expect_equal(fit$case, deg)
        expect_true(all(is.finite(coef(fit))))
        expect_true(all(is.finite(sqrt(diag(fit$vcov3)))))
        # the dropped arm's regression parameters must be absent
        dropped <- if (deg == "z0deg") "^d0" else "^d1"
        expect_false(any(grepl(dropped, names(fit$theta))))
        if (deg == "z0deg") expect_equal(fit$dmeanz0, 0)
        if (deg == "z1deg") expect_equal(fit$dmeanz1, 1)
      })
    }
  }
}

test_that("degenerate IPWRA estimate matches closed form (z0deg)", {
  d <- onesided_z0()
  fit <- drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d)
  ps <- fitted(glm(rsncode ~ age, binomial, data = d))
  m <- model.matrix(~ age, d)
  s1 <- d$rsncode == 1
  da <- transform(d, ww1 = 1 / ps, ww0 = 1 / (1 - ps))
  fy1 <- lm(lwage ~ age, da[s1, ], weights = ww1)
  fy0 <- lm(lwage ~ age, da[!s1, ], weights = ww0)
  fd1 <- suppressWarnings(glm(nvstat ~ age, quasibinomial, da[s1, ],
                              weights = ww1))
  num <- mean(m %*% coef(fy1)) - mean(m %*% coef(fy0))
  den <- mean(plogis(m %*% coef(fd1))) - 0
  expect_equal(unname(coef(fit)), c(num / den, num, den), tolerance = 1e-8)
})

test_that("both-degenerate (perfect compliance) errors informatively", {
  d <- drlate_sim
  d$nvstat <- d$rsncode
  expect_error(
    drlate(lwage ~ age, nvstat ~ age, rsncode ~ age, data = d,
           method = "aipw"),
    "degenerate in both"
  )
})
