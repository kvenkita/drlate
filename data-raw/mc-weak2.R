# Genuinely weak instrument: Fieller vs Wald coverage (true LATE = 0.5)
setwd("C:/Users/kyle/Documents/Projects/Personal/drlate")
devtools::load_all(quiet = TRUE)
set.seed(20260606)

gen_weak <- function(n = 500, pcomp = 0.05) {
  age <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.3 * age))
  comp <- rbinom(n, 1, pcomp)
  at <- rbinom(n, 1, 0.10)
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

R <- 400
cov_f <- cov_w <- fz <- rep(NA, R)
types <- character(R)
for (r in seq_len(R)) {
  dat <- gen_weak()
  fit <- tryCatch(drlate(y ~ age, d ~ age, z ~ age, data = dat),
                  error = function(e) NULL)
  if (is.null(fit)) next
  f <- confint(fit, method = "fieller")
  wald <- confint(fit)[1, ]
  cov_f[r] <- in_fieller(f, 0.5)
  cov_w[r] <- wald[1] <= 0.5 && 0.5 <= wald[2]
  fz[r] <- drlate:::firststage_z(fit)
  types[r] <- f$type
}
cat(sprintf("WEAK-IV-2 (mean |z| = %.2f, ok %d/%d):\n",
            mean(abs(fz), na.rm = TRUE), sum(!is.na(cov_f)), R))
cat(sprintf("  Fieller 95%% coverage: %.3f   Wald 95%% coverage: %.3f\n",
            mean(cov_f, na.rm = TRUE), mean(cov_w, na.rm = TRUE)))
print(table(types[types != ""]))
cat("DONE\n")
