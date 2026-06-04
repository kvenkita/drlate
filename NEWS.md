# drlate 0.1.0

Extensions beyond the Stata original:

* **Diagnostics:** `plot()` methods for instrument propensity score overlap,
  covariate balance (love plot), and implied-weight distributions;
  `balance()` returns standardized mean differences; `print()`/`summary()`
  report first-stage strength and flag weak instruments.
* **Bootstrap inference:** `drlate(vcov = "bootstrap")` provides
  nonparametric bootstrap standard errors and percentile confidence
  intervals (cluster bootstrap when `cluster` is supplied), with optional
  parallelism.
* **Weak-instrument-robust inference:** `confint(method = "fieller")`
  inverts the joint test of the numerator and denominator, returning
  bounded, complement, or whole-line confidence sets as appropriate.
* **DR Hausman test:** `dr_hausman()` implements the doubly robust test of
  unconfoundedness from Słoczyński, Uysal & Wooldridge (2022, Section 5) —
  proposed in the paper but not available in the Stata package — with an
  analytic standard error from one jointly stacked moment system.
* **Estimator comparison:** `drlate_compare()` runs IPWRA/IPW/AIPW/RA in
  one call with a comparison table and dot-whisker plot.
* **Documentation site:** pkgdown website with the primer vignette
  featured, deployed via GitHub Actions.

# drlate 0.0.0.9000

* Initial R port of the Stata package `drlate` v1.0.0 (SSC S459708).
* LATE and LATT estimation via IPWRA, IPW, AIPW, and RA.
* Linear, logit, and Poisson outcome/treatment models.
* Logit (MLE), CBPS, and IPT instrument propensity score models.
* Joint sandwich inference across all estimation stages, with robust and
  cluster-robust variants and sampling weights.
* Overlap diagnostics (`pstolerance`, `osample`).
* Golden-fixture test harness for numerical equivalence with Stata:
  estimates match to ~1e-9 and standard errors to ~1e-6 across all 33
  validation scenarios.
