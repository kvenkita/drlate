# Build the estimation context (`ctx`) from the user's three formulas.
# The ctx is an immutable list passed to every downstream component.

#' @noRd
build_ctx <- function(outcome, treatment, instrument, data,
                      omodel, tmodel, ivmodel, method, estimand, normalized,
                      weights = NULL, cluster = NULL,
                      pstolerance = 1e-5, osample = FALSE) {
  stopifnot(inherits(outcome, "formula"), inherits(treatment, "formula"),
            inherits(instrument, "formula"))

  # --- Evaluate the three model frames on the full data (na.pass), then
  # --- intersect complete cases so all equations share one sample
  # --- (mirrors marksample touse + markout in drlate_estimate.ado).
  mf_y <- stats::model.frame(outcome,    data, na.action = stats::na.pass)
  mf_d <- stats::model.frame(treatment,  data, na.action = stats::na.pass)
  mf_z <- stats::model.frame(instrument, data, na.action = stats::na.pass)
  n0 <- nrow(data)
  if (nrow(mf_y) != n0 || nrow(mf_d) != n0 || nrow(mf_z) != n0) {
    stop("internal error: model frames have differing lengths", call. = FALSE)
  }

  ok <- stats::complete.cases(mf_y) & stats::complete.cases(mf_d) &
        stats::complete.cases(mf_z)
  if (!is.null(weights)) ok <- ok & !is.na(weights)
  if (!is.null(cluster)) ok <- ok & !is.na(cluster)

  keep <- which(ok)
  n <- length(keep)
  if (n == 0L) stop("no complete observations.", call. = FALSE)

  y <- as.numeric(stats::model.response(mf_y))[keep]
  d <- as.numeric(stats::model.response(mf_d))[keep]
  z <- as.numeric(stats::model.response(mf_z))[keep]

  w <- if (is.null(weights)) rep(1, n) else as.numeric(weights)[keep]
  if (any(w < 0)) stop("negative weights are not allowed.", call. = FALSE)
  cl <- if (is.null(cluster)) NULL else cluster[keep]

  # --- Model matrices (with intercept; factors expand via model.matrix) ---
  Xo <- stats::model.matrix(stats::terms(outcome),    mf_y[keep, , drop = FALSE])
  Xt <- stats::model.matrix(stats::terms(treatment),  mf_d[keep, , drop = FALSE])
  Xz <- stats::model.matrix(stats::terms(instrument), mf_z[keep, , drop = FALSE])

  # --- Method-specific covariate restrictions (drlate_estimate.ado) ---
  if (method == "ipw" && (ncol(Xo) > 1L || ncol(Xt) > 1L)) {
    stop("covariates are not allowed in the outcome or treatment equations ",
         "with method = \"ipw\"; use `y ~ 1` and `d ~ 1`.", call. = FALSE)
  }
  if (method == "ra" && ncol(Xz) > 1L) {
    stop("covariates are not allowed in the instrument equation with ",
         "method = \"ra\"; use `z ~ 1`.", call. = FALSE)
  }
  if (method %in% c("kappa", "kappa0", "kappa10")) {
    if (estimand == "latt") {
      stop("the kappa-weighting estimators are available for ",
           "estimand = \"late\" only (Sloczynski, Uysal, and Wooldridge ",
           "2025).", call. = FALSE)
    }
    if (ivmodel == "ipt") {
      stop("ivmodel = \"ipt\" is not available with the kappa-weighting ",
           "estimators.", call. = FALSE)
    }
    if (method %in% c("kappa0", "kappa10") && ivmodel == "cbps") {
      stop("ivmodel = \"cbps\" is available only with method = \"kappa\" ",
           "among the kappa-weighting estimators, following the Stata ",
           "kappalate command.", call. = FALSE)
    }
    if (ncol(Xo) > 1L || ncol(Xt) > 1L) {
      stop("covariates are not allowed in the outcome or treatment ",
           "equations with the kappa-weighting estimators; use `y ~ 1` ",
           "and `d ~ 1`.", call. = FALSE)
    }
    check_binary(d, all.vars(treatment)[1L], "treatment")
  }
  if (ivmodel == "cbps" && estimand == "latt") {
    stop("ivmodel = \"cbps\" is not available with estimand = \"latt\".",
         call. = FALSE)
  }
  if (ivmodel == "probit" &&
      !(method %in% c("ipw", "kappa", "kappa0", "kappa10") &&
        estimand == "late")) {
    stop("ivmodel = \"probit\" is available only for the weighting ",
         "estimators covered by the Stata kappalate command (method ",
         "\"ipw\", \"kappa\", \"kappa0\", or \"kappa10\") with ",
         "estimand = \"late\".", call. = FALSE)
  }

  # --- Input validation (drlate_estimate.ado section 4) ---
  zname <- all.vars(instrument)[1L]
  check_binary(z, zname, "instrument")
  check_family(y, all.vars(outcome)[1L],   omodel, "outcome")
  check_family(d, all.vars(treatment)[1L], tmodel, "treatment")

  # --- Standardize continuous covariate columns (span-preserving) ---
  Xo <- standardize_mm(Xo, w)
  Xt <- standardize_mm(Xt, w)
  Xz <- standardize_mm(Xz, w)

  # --- Compliance means (drlate_estimate.ado section 5) ---
  # Stata uses a plain `summarize` here, i.e. UNWEIGHTED means even under
  # pweights. Harmless for estimation (the values only enter the moment
  # conditions when exactly 0 or 1, where weighting is irrelevant), but
  # replicated so the reported dmeanz1/dmeanz0 match e() exactly.
  dmeanz1 <- mean(d[z == 1])
  dmeanz0 <- mean(d[z == 0])
  case <- if (dmeanz0 %in% c(0, 1) && dmeanz1 %in% c(0, 1)) "bothdeg"
          else if (dmeanz0 %in% c(0, 1)) "z0deg"
          else if (dmeanz1 %in% c(0, 1)) "z1deg"
          else "interior"

  # IPWRA and RA always use normalized moments; IPT weights are ex-ante
  # normalized (drlate_estimate_late.ado lines 54-56, 658-661).
  statnorm <- if (normalized) "nrm" else "unnrm"
  if (method == "ipwra") statnorm <- "nrm"

  list(
    y = y, d = d, z = z, w = w, cluster = cl, n = n,
    Xo = Xo, Xt = Xt, Xz = Xz,
    omodel = omodel, tmodel = tmodel, ivmodel = ivmodel,
    method = method, estimand = estimand, statnorm = statnorm,
    dmeanz1 = dmeanz1, dmeanz0 = dmeanz0, case = case,
    pstolerance = pstolerance, osample = osample
  )
}
