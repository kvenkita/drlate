# Complier profiling: Abadie's kappa weights and complier covariate means
# (the estat compliers postestimation feature of Stata's lateffects).

#' Abadie's kappa weights
#'
#' Returns the per-observation Abadie kappa weight implied by a fitted
#' [drlate()] object,
#' \deqn{\kappa = 1 - \frac{D(1 - Z)}{1 - p(X)} - \frac{(1 - D) Z}{p(X)},}
#' where \eqn{p(X)} is the estimated instrument propensity score. The kappa
#' weights identify the complier subpopulation: for any function \eqn{g} of the
#' data, \eqn{E[g \mid \mathrm{complier}] = E[\kappa g] / E[\kappa]}
#' (Abadie 2003). They are the weights used by [complier_means()] and are the
#' \proglang{Stata} \code{estat compliers, genkappa()} object.
#'
#' @param object A fitted [drlate()] object (with `keep_data = TRUE`) using an
#'   instrument propensity score (any `method` except `"ra"`).
#' @param normalize Logical. If `TRUE` (default), the returned weights are the
#'   sampling-weighted, normalized weights \eqn{w\kappa / \sum w\kappa} that sum
#'   to one (the form used to compute complier averages). If `FALSE`, the raw
#'   kappa values are returned.
#' @return A numeric vector with one entry per estimation-sample observation.
#' @seealso [complier_means()]
#' @examples
#' fit <- drlate(lwage ~ age + educ, nvstat ~ age + educ,
#'               rsncode ~ age + educ, data = drlate_sim)
#' head(kappa_weights(fit))
#' @export
kappa_weights <- function(object, normalize = TRUE) {
  ctx <- need_ctx(object)
  ps <- need_ps(object)
  z <- ctx$z; d <- ctx$d
  p <- ps$ps
  kap <- 1 - d * (1 - z) / (1 - p) - (1 - d) * z / p
  if (normalize) {
    wk <- ctx$w * kap
    return(wk / sum(wk))
  }
  kap
}

#' Complier covariate means
#'
#' Compares the average of each covariate in the full estimation sample with
#' its average in the complier subpopulation, the latter computed with the
#' normalized Abadie kappa weights of [kappa_weights()]. Because the local
#' average treatment effect is a causal effect for compliers, knowing how
#' compliers differ from the population aids interpretation. This is the
#' \proglang{Stata} \code{estat compliers} postestimation feature.
#'
#' Covariate values are reported on their original scale.
#'
#' @param object A fitted [drlate()] object (with `keep_data = TRUE`) using an
#'   instrument propensity score (any `method` except `"ra"`).
#' @param vars Optional character vector selecting a subset of the model
#'   covariates. Defaults to all covariates across the three model formulas.
#' @return A data frame with one row per covariate and columns `variable`,
#'   `population_mean`, `complier_mean`, and `difference`
#'   (`complier_mean - population_mean`).
#' @seealso [kappa_weights()]
#' @examples
#' fit <- drlate(lwage ~ age + educ, nvstat ~ age + educ,
#'               rsncode ~ age + educ, data = drlate_sim)
#' complier_means(fit)
#' @export
complier_means <- function(object, vars = NULL) {
  ctx <- need_ctx(object)
  need_ps(object)
  X <- diag_covariates(ctx)
  if (!is.null(vars)) {
    missing <- setdiff(vars, colnames(X))
    if (length(missing)) {
      stop("unknown covariate(s): ", paste(missing, collapse = ", "),
           ". Available: ", paste(colnames(X), collapse = ", "), ".",
           call. = FALSE)
    }
    X <- X[, vars, drop = FALSE]
  }
  w <- ctx$w
  wk <- w * kappa_weights(object, normalize = FALSE)
  pop  <- colSums(w * X)  / sum(w)
  comp <- colSums(wk * X) / sum(wk)
  data.frame(
    variable = colnames(X),
    population_mean = unname(pop),
    complier_mean = unname(comp),
    difference = unname(comp - pop),
    row.names = NULL
  )
}
