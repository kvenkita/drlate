#' @keywords internal
"_PACKAGE"

#' @importFrom stats glm lm.wfit binomial quasibinomial quasipoisson gaussian
#'   plogis qnorm pnorm model.frame model.matrix model.response model.weights
#'   na.omit terms coef vcov confint nobs printCoefmat setNames weighted.mean
NULL

# The .data pronoun used in ggplot2 aes() mappings (ggplot2 is in Suggests)
utils::globalVariables(".data")
