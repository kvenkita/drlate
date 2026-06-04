#' Simulated example data for drlate
#'
#' A simulated dataset with a binary instrument, a binary treatment with
#' two-sided noncompliance, and continuous, positive, and binary outcome
#' variables, designed to exercise every model family supported by
#' [drlate()]. The complier average treatment effect (LATE) used in the
#' data-generating process is 0.5.
#'
#' @format A data frame with 2,000 rows and 7 variables:
#' \describe{
#'   \item{lwage}{continuous outcome}
#'   \item{kwage}{positive outcome (for Poisson models), `exp(lwage / 2)`}
#'   \item{hijob}{binary outcome (for logit models)}
#'   \item{nvstat}{binary treatment}
#'   \item{rsncode}{binary instrument}
#'   \item{age}{continuous covariate}
#'   \item{educ}{factor covariate with levels `hs`, `college`, `graduate`}
#' }
#'
#' @source Simulated; see `data-raw/drlate_sim.R` in the package sources.
"drlate_sim"
