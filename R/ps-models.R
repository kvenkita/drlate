# Instrument propensity score estimation: logit MLE, CBPS, IPT.
# Mirrors _drlate_ps.ado. CBPS/IPT solve just-identified moment systems by
# damped Newton, warm-started from the logit MLE (Stata starts its one-step
# GMM from the logit fit).

#' Fit the instrument propensity score; returns coefficients and fitted ps.
#' @noRd
fit_ps <- function(ctx) {
  z <- ctx$z; Xz <- ctx$Xz; w <- ctx$w

  if (ctx$ivmodel == "probit") {
    # Probit MLE (kappalate zmodel(probit)); no logit warm start needed.
    # Tight IRLS tolerance: the unnormalized weighting estimators amplify
    # propensity-score rounding, and Stata's Newton iterates further than
    # glm.fit's default epsilon.
    fitp <- suppressWarnings(
      stats::glm.fit(Xz, z, weights = w,
                     family = stats::quasibinomial(link = "probit"),
                     control = stats::glm.control(epsilon = 1e-12,
                                                  maxit = 100L))
    )
    if (!fitp$converged) {
      stop("convergence not achieved for probit instrument propensity ",
           "score estimation.", call. = FALSE)
    }
    b <- stats::coef(fitp)
    ps <- stats::pnorm(drop(Xz %*% b))
    return(list(bips = b, ps = ps,
                wt1 = z / ps, wt0 = (1 - z) / (1 - ps)))
  }

  # Logit MLE start (and the full answer when ivmodel == "logit")
  fit <- suppressWarnings(
    stats::glm.fit(Xz, z, weights = w, family = stats::quasibinomial())
  )
  if (!fit$converged) {
    stop("convergence not achieved for logit instrument propensity score ",
         "estimation.", call. = FALSE)
  }
  b_logit <- stats::coef(fit)

  if (ctx$ivmodel == "logit") {
    ps <- stats::plogis(drop(Xz %*% b_logit))
    return(list(bips = b_logit, ps = ps,
                wt1 = z / ps, wt0 = (1 - z) / (1 - ps)))
  }

  if (ctx$ivmodel == "cbps") {
    # Balancing moment: (z/p - (1-z)/(1-p)) * x = 0
    gmom <- function(b) {
      p <- stats::plogis(drop(Xz %*% b))
      colSums(w * (z / p - (1 - z) / (1 - p)) * Xz) / sum(w)
    }
    b <- newton_solve(gmom, b_logit)
    ps <- stats::plogis(drop(Xz %*% b))
    return(list(bips = b, ps = ps,
                wt1 = z / ps, wt0 = (1 - z) / (1 - ps)))
  }

  if (ctx$ivmodel == "ipt") {
    # Tilted PS for the Z=1 arm: (z * (1 + exp(xb)) / exp(xb) - 1) * x = 0
    g1 <- function(b) {
      xb <- drop(Xz %*% b)
      colSums(w * (z * (1 + exp(xb)) / exp(xb) - 1) * Xz) / sum(w)
    }
    # Tilted PS for the Z=0 arm: ((1-z) * (1 + exp(xb)) - 1) * x = 0
    g0 <- function(b) {
      xb <- drop(Xz %*% b)
      colSums(w * ((1 - z) * (1 + exp(xb)) - 1) * Xz) / sum(w)
    }
    b1 <- newton_solve(g1, b_logit)
    b0 <- newton_solve(g0, b_logit)
    ps1 <- stats::plogis(drop(Xz %*% b1))
    ps0 <- stats::plogis(drop(Xz %*% b0))
    # Combined ps for the overlap check (_drlate_ps.ado)
    ps <- ifelse(z == 1, ps1, ps0)
    return(list(bips1 = b1, bips0 = b0, ps1 = ps1, ps0 = ps0, ps = ps,
                wt1 = z / ps1, wt0 = (1 - z) / (1 - ps0)))
  }

  stop("unknown ivmodel: ", ctx$ivmodel)
}

#' Damped Newton solver for small just-identified moment systems.
#' Jacobian by numDeriv; step-halving on non-decrease of the moment norm.
#' @noRd
newton_solve <- function(gmom, start, tol = 1e-10, maxit = 100L) {
  b <- start
  gnorm <- function(b) sqrt(sum(gmom(b)^2))
  f <- gnorm(b)
  for (it in seq_len(maxit)) {
    if (f < tol) return(b)
    J <- numDeriv::jacobian(gmom, b)
    step <- tryCatch(solve(J, gmom(b)), error = function(e) NULL)
    if (is.null(step)) {
      stop("singular Jacobian in moment solver; convergence not achieved.",
           call. = FALSE)
    }
    lambda <- 1
    repeat {
      b_new <- b - lambda * step
      f_new <- gnorm(b_new)
      if (f_new < f || lambda < 1e-8) break
      lambda <- lambda / 2
    }
    if (lambda < 1e-8 && f_new >= f) {
      stop("convergence not achieved in moment solver.", call. = FALSE)
    }
    b <- b_new; f <- f_new
  }
  if (f >= tol) {
    stop("convergence not achieved in moment solver after ", maxit,
         " iterations.", call. = FALSE)
  }
  b
}
