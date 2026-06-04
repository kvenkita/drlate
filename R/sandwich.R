# Joint M-estimation sandwich variance over the stacked moment system.
# Replicates Stata `gmm ..., onestep winitial(I) from(theta0) quickderivatives
# vce(robust|cluster) iterate(0)` on a just-identified system:
#   V = A^{-1} B A^{-T} / n
# with A the (numeric) Jacobian of the weighted averaged moments and B the
# (cluster-)robust outer product of the weighted per-observation moments.

#' @noRd
drlate_vcov <- function(sys, theta, w, cluster = NULL) {
  n <- length(w)
  gbar <- function(t) colMeans(w * sys$g(t))
  A <- numDeriv::jacobian(gbar, theta)
  Ainv <- tryCatch(solve(A), error = function(e) {
    stop("singular Jacobian in the stacked moment system; the variance ",
         "matrix could not be computed.", call. = FALSE)
  })
  G <- w * sys$g(theta)
  if (is.null(cluster)) {
    B <- crossprod(G) / n
  } else {
    # Stata's gmm vce(cluster) applies NO small-sample factor (verified
    # against fixtures: an M/(M-1) correction overstates SEs by exactly
    # sqrt(M/(M-1))).
    Gc <- rowsum(G, group = cluster)
    B <- crossprod(Gc) / n
  }
  V <- Ainv %*% B %*% t(Ainv) / n
  dimnames(V) <- list(names(theta), names(theta))
  V
}
