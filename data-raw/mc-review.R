# Monte Carlo evidence for the peer review of the added features.
# (Dev script, not part of the package build.)
setwd("C:/Users/kyle/Documents/Projects/Personal/drlate")
devtools::load_all(quiet = TRUE)
set.seed(20260605)

## ---------------------------------------------------------------------
## 1. Size of the DR Hausman test under a TRUE null
##    (one-sided noncompliance, treatment unconfounded given X:
##    compliance depends only on X, baseline outcome independent of type)
## ---------------------------------------------------------------------
gen_h0 <- function(n = 1500) {
  age <- rnorm(n)
  educ <- rbinom(n, 1, 0.4)
  z <- rbinom(n, 1, plogis(-0.3 + 0.6 * age + 0.5 * educ))
  ctype <- rbinom(n, 1, plogis(0.2 + 0.5 * age))     # complier given X
  d <- z * ctype                                      # no always-takers
  y <- 1 + 0.8 * age + 0.5 * educ + 0.5 * d + rnorm(n)
  data.frame(y, d, z, age, educ)
}

R1 <- 300
rej <- rep(NA, R1)
for (r in seq_len(R1)) {
  dat <- gen_h0()
  h <- tryCatch(
    dr_hausman(y ~ age + educ, d ~ age + educ, z ~ age + educ, data = dat),
    error = function(e) NULL)
  if (!is.null(h)) rej[r] <- h$p.value < 0.05
}
cat(sprintf("HAUSMAN size at nominal 5%%: %.3f (ok reps: %d/%d)\n",
            mean(rej, na.rm = TRUE), sum(!is.na(rej)), R1))

## ---------------------------------------------------------------------
## 2. Coverage of Fieller vs Wald under a WEAK instrument (true LATE 0.5)
## ---------------------------------------------------------------------
gen_weak <- function(n = 800, pcomp = 0.12) {
  age <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.3 * age))
  comp <- rbinom(n, 1, pcomp)                # few compliers -> weak stage
  at <- rbinom(n, 1, 0.10)                   # some always-takers
  d <- ifelse(at == 1, 1L, z * comp)
  y <- 1 + 0.5 * age + 0.5 * d + rnorm(n)
  data.frame(y, d, z, age)
}

in_fieller <- function(f, t) {
  switch(f$type,
    bounded = t >= f$lower && t <= f$upper,
    complement = t <= f$lower || t >= f$upper,
    `whole-line` = TRUE)
}

R2 <- 300
cov_f <- cov_w <- fz <- rep(NA, R2)
types <- character(R2)
for (r in seq_len(R2)) {
  dat <- gen_weak()
  fit <- tryCatch(
    drlate(y ~ age, d ~ age, z ~ age, data = dat),
    error = function(e) NULL)
  if (is.null(fit)) next
  f <- confint(fit, method = "fieller")
  wald <- confint(fit)[1, ]
  cov_f[r] <- in_fieller(f, 0.5)
  cov_w[r] <- wald[1] <= 0.5 && 0.5 <= wald[2]
  fz[r] <- drlate:::firststage_z(fit)
  types[r] <- f$type
}
cat(sprintf("WEAK-IV (mean |first-stage z| = %.2f):\n",
            mean(abs(fz), na.rm = TRUE)))
cat(sprintf("  Fieller 95%% coverage: %.3f   Wald 95%% coverage: %.3f\n",
            mean(cov_f, na.rm = TRUE), mean(cov_w, na.rm = TRUE)))
cat("  Fieller set types: ")
print(table(types[types != ""]))

## ---------------------------------------------------------------------
## 3. Analytic SE accuracy and CI coverage at moderate strength
##    (sanity check of the joint sandwich on a fresh DGP, plus bootstrap)
## ---------------------------------------------------------------------
gen_mid <- function(n = 1000) {
  age <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.4 * age))
  type <- sample(c("c", "n", "a"), n, TRUE, c(0.5, 0.35, 0.15))
  d <- ifelse(type == "a", 1L, ifelse(type == "n", 0L, z))
  y <- 1 + 0.6 * age + 0.3 * (type == "a") + 0.5 * d * (type == "c") +
       0.2 * d * (type != "c") + rnorm(n)
  data.frame(y, d, z, age)
}

R3 <- 300
est <- se <- covr <- rep(NA, R3)
for (r in seq_len(R3)) {
  dat <- gen_mid()
  fit <- tryCatch(drlate(y ~ age, d ~ age, z ~ age, data = dat),
                  error = function(e) NULL)
  if (is.null(fit)) next
  est[r] <- coef(fit)[1]
  se[r] <- sqrt(fit$vcov3[1, 1])
  ci <- confint(fit)[1, ]
  covr[r] <- ci[1] <= 0.5 && 0.5 <= ci[2]
}
cat(sprintf(paste0("SANDWICH check: MC sd of estimate = %.4f, ",
                   "mean analytic SE = %.4f, 95%% CI coverage = %.3f\n"),
            sd(est, na.rm = TRUE), mean(se, na.rm = TRUE),
            mean(covr, na.rm = TRUE)))

# Bootstrap SE vs analytic on one dataset
dat <- gen_mid(2000)
fb <- drlate(y ~ age, d ~ age, z ~ age, data = dat,
             vcov = "bootstrap", boot_reps = 499, boot_seed = 1)
cat(sprintf("BOOTSTRAP: boot SE = %.4f vs analytic SE = %.4f (n=2000)\n",
            fb$boot$se[1], sqrt(fb$vcov3[1, 1])))
cat("MC DONE\n")
