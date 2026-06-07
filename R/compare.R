#' Compare drlate estimators in one call
#'
#' Runs several estimators on the same specification and collects the
#' causal estimates with their confidence intervals — the sensitivity
#' comparison applied papers routinely report. Formula restrictions are
#' handled automatically: `method = "ipw"` drops the outcome/treatment
#' covariates and `method = "ra"` drops the instrument covariates (each
#' with a message), matching the requirements of those estimators.
#'
#' @details
#' Because IPW carries no outcome/treatment regressions and RA carries no
#' instrument propensity score, the automatic formula adjustment means
#' the rows do not share a single adjustment specification: differences
#' between the IPW or RA row and the doubly robust rows reflect both the
#' estimator *and* the reduced specification. Read the comparison as a
#' robustness display, not as a test that isolates estimator choice; the
#' doubly robust rows (IPWRA, AIPW) are the like-for-like pair.
#'
#' @inheritParams drlate
#' @param methods Estimators to run (any of the `method` values accepted by [drlate()]).
#' @param both_norms Logical; also run the unnormalized variants of
#'   `"ipw"` and `"aipw"` (default `FALSE`).
#' @param ... Passed on to [drlate()] (e.g. `omodel`, `tmodel`, `ivmodel`,
#'   `estimand`, `weights`, `cluster`).
#'
#' @return An object of class `"drlate_compare"`: a data frame with columns
#'   `method`, `normalized`, `estimate`, `se`, `ci_lo`, `ci_hi`, with a
#'   `print` method and a dot-whisker `plot` method.
#'
#' @examples
#' cmp <- drlate_compare(lwage ~ age + educ, nvstat ~ age + educ,
#'                       rsncode ~ age + educ, data = drlate_sim)
#' cmp
#'
#' @export
drlate_compare <- function(outcome, treatment, instrument, data,
                           methods = c("ipwra", "ipw", "aipw", "ra"),
                           both_norms = FALSE, ...) {
  methods <- match.arg(methods, c("ipwra", "ipw", "aipw", "ra",
                                  "kappa", "kappa0", "kappa10"),
                       several.ok = TRUE)
  lhs <- function(f) f[[2]]

  specs <- list()
  for (me in methods) {
    norms <- if (both_norms && me %in% c("ipw", "aipw")) c(TRUE, FALSE)
             else TRUE
    for (nr in norms) specs[[length(specs) + 1L]] <- list(me = me, nr = nr)
  }

  rows <- lapply(specs, function(sp) {
    fo <- outcome; ft <- treatment; fz <- instrument
    if (sp$me %in% c("ipw", "kappa", "kappa0", "kappa10")) {
      if (length(all.vars(fo)) > 1L || length(all.vars(ft)) > 1L) {
        message("method = \"", sp$me, "\": dropping outcome/treatment ",
                "covariates (weighted means only).")
      }
      fo <- stats::as.formula(call("~", lhs(outcome), 1))
      ft <- stats::as.formula(call("~", lhs(treatment), 1))
    }
    if (sp$me == "ra" && length(all.vars(fz)) > 1L) {
      message("method = \"ra\": dropping instrument covariates ",
              "(no propensity score).")
      fz <- stats::as.formula(call("~", lhs(instrument), 1))
    }
    fit <- tryCatch(
      drlate(fo, ft, fz, data = data, method = sp$me,
             normalized = sp$nr, keep_data = FALSE, ...),
      error = function(e) {
        warning("method \"", sp$me, "\" failed: ", conditionMessage(e),
                call. = FALSE)
        NULL
      }
    )
    if (is.null(fit)) {
      return(data.frame(method = sp$me, normalized = sp$nr,
                        estimate = NA_real_, se = NA_real_,
                        ci_lo = NA_real_, ci_hi = NA_real_))
    }
    ci <- confint(fit)[1, ]
    data.frame(method = sp$me, normalized = fit$statnorm == "nrm",
               estimate = unname(coef(fit)[1]),
               se = sqrt(fit$vcov3[1, 1]),
               ci_lo = unname(ci[1]), ci_hi = unname(ci[2]))
  })

  out <- do.call(rbind, rows)
  # The normalize-check can switch a requested normalized fit to
  # unnormalized (e.g. with IPT weights); drop the resulting duplicates.
  out <- out[!duplicated(out[, c("method", "normalized")]), , drop = FALSE]
  rownames(out) <- NULL
  class(out) <- c("drlate_compare", "data.frame")
  attr(out, "estimand") <- list(...)$estimand %||% "late"
  out
}

#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' @export
print.drlate_compare <- function(x, digits = 4, ...) {
  est <- toupper(attr(x, "estimand") %||% "late")
  cat("Estimator comparison (", est, ")\n\n", sep = "")
  lab <- ifelse(x$method %in% c("ipw", "aipw"),
                paste0(x$method, ifelse(x$normalized, " (nrm)", " (unnrm)")),
                x$method)
  tab <- data.frame(estimator = lab,
                    estimate = round(x$estimate, digits),
                    se = round(x$se, digits),
                    `95% CI` = sprintf("[%s, %s]",
                                       format(x$ci_lo, digits = digits),
                                       format(x$ci_hi, digits = digits)),
                    check.names = FALSE)
  print(tab, row.names = FALSE)
  invisible(x)
}

#' @export
plot.drlate_compare <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("package 'ggplot2' is required for this plot.", call. = FALSE)
  }
  df <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  df$label <- ifelse(df$method %in% c("ipw", "aipw"),
                     paste0(toupper(df$method),
                            ifelse(df$normalized, " (nrm)", " (unnrm)")),
                     toupper(df$method))
  df$label <- factor(df$label, levels = rev(unique(df$label)))
  est <- toupper(attr(x, "estimand") %||% "late")
  ggplot2::ggplot(df, ggplot2::aes(x = .data$estimate, y = .data$label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        linewidth = 0.3) +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = .data$ci_lo,
                                          xmax = .data$ci_hi),
                             color = "#2c6e8f") +
    ggplot2::labs(x = paste(est, "estimate (95% CI)"), y = NULL,
                  title = "Sensitivity to estimator choice") +
    ggplot2::theme_minimal()
}
