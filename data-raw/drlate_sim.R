# Generate the simulated example dataset bundled with the package.
# Mimics the structure of the SIPP extract used in the Stata drlate help file:
# a binary instrument, a binary treatment with two-sided noncompliance, and
# continuous / binary / count outcomes so every model family can be exercised.
set.seed(20260604)

n <- 2000
age   <- round(rnorm(n, 35, 8))
educ  <- sample(c("hs", "college", "graduate"), n, replace = TRUE,
                prob = c(0.5, 0.35, 0.15))
x     <- (age - 35) / 8

# Binary instrument with covariate-dependent propensity (bounded well inside (0,1))
pz    <- plogis(-0.2 + 0.5 * x + 0.4 * (educ == "college") + 0.6 * (educ == "graduate"))
z     <- rbinom(n, 1, pz)

# Compliance types: 60% compliers, 25% never-takers, 15% always-takers
type  <- sample(c("c", "n", "a"), n, replace = TRUE, prob = c(0.6, 0.25, 0.15))
d     <- ifelse(type == "a", 1L, ifelse(type == "n", 0L, z))

# Potential outcomes: LATE (complier effect) = 0.5
y0    <- 1 + 0.8 * x + 0.3 * (educ != "hs") + rnorm(n)
tau   <- ifelse(type == "c", 0.5, 0.2)
y     <- y0 + tau * d

drlate_sim <- data.frame(
  lwage  = y,
  kwage  = exp(y / 2),                          # positive outcome for poisson
  hijob  = as.integer(y > median(y)),           # binary outcome for logit
  nvstat = d,
  rsncode = z,
  age    = age,
  educ   = factor(educ, levels = c("hs", "college", "graduate"))
)

usethis::use_data(drlate_sim, overwrite = TRUE)
