# Dev smoke test (not part of the package build)
setwd("C:/Users/kyle/Documents/Projects/Personal/drlate")
devtools::load_all(quiet = TRUE)

fit <- drlate(lwage ~ age + educ, nvstat ~ age + educ, rsncode ~ age + educ,
              data = drlate_sim)
print(fit)
cat("\ncase:", fit$case, " dmeanz1:", round(fit$dmeanz1, 4),
    " dmeanz0:", round(fit$dmeanz0, 4), "\n")

# Internal consistency check 1: with z ~ 1 the IPWRA collapses toward the
# simple Wald/2SLS estimand; compare against the raw Wald ratio.
wald <- with(drlate_sim,
  (mean(lwage[rsncode == 1]) - mean(lwage[rsncode == 0])) /
  (mean(nvstat[rsncode == 1]) - mean(nvstat[rsncode == 0])))
cat("raw Wald ratio:", round(wald, 6), "\n")

# Internal consistency check 2: point estimate must equal the closed-form
# sequential computation done independently here.
ps <- fitted(glm(rsncode ~ age + educ, binomial, data = drlate_sim))
w1 <- 1 / ps; w0 <- 1 / (1 - ps)
m  <- model.matrix(~ age + educ, drlate_sim)
s1 <- drlate_sim$rsncode == 1
da <- transform(drlate_sim, ww1 = w1, ww0 = w0)
fy1 <- lm(lwage ~ age + educ, da[s1, ], weights = ww1)
fy0 <- lm(lwage ~ age + educ, da[!s1, ], weights = ww0)
fd1 <- suppressWarnings(glm(nvstat ~ age + educ, quasibinomial, da[s1, ],
                            weights = ww1))
fd0 <- suppressWarnings(glm(nvstat ~ age + educ, quasibinomial, da[!s1, ],
                            weights = ww0))
num <- mean(m %*% coef(fy1)) - mean(m %*% coef(fy0))
den <- mean(plogis(m %*% coef(fd1))) - mean(plogis(m %*% coef(fd0)))
cat("independent IPWRA late:", num / den, "  drlate:",
    unname(coef(fit)[1]), "\n")
stopifnot(abs(num / den - coef(fit)[1]) < 1e-10)
cat("OK: point estimate matches independent computation\n")
