# Helpers shared by the estimator test files.

#' Assert the stacked system is square with zero averaged moments,
#' then return the fitted drlate object.
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
