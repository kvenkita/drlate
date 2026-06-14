# drlate (development version)

New outcome and treatment model families, completing parity with the Stata
`lateffects` `omodel`/`tmodel` options:

* `omodel` gains `"probit"` (binary outcome, probit link) and `"flogit"` /
  `"fprobit"` for **fractional** outcomes in `[0, 1]` (e.g. proportions or
  rates). `tmodel` gains `"probit"`. The fractional families share all
  estimation with their binary counterparts and only relax the response to the
  unit interval. These families are internally validated (the probit
  quasi-likelihood score is the one already used for the probit instrument
  propensity score); direct numerical parity with `lateffects` is planned for a
  later release.

Postestimation diagnostics mirroring the Stata `lateffects` suite (StataNow):

* `complier_means()` reports population versus complier covariate means, the
  complier averages computed with the normalized Abadie-kappa weights (Stata's
  `estat compliers`). `kappa_weights()` returns those weights (the `genkappa`
  object) for use in other complier summaries.
* `balance_test()` implements the Imai and Ratkovic (2014) overidentification
  test for whether the instrument propensity score balances the covariates
  (Stata's `latebalance overid`); cluster-robust when the fit is.
* `balance(detail = TRUE)` adds IPW-weighted arm means and unweighted and
  weighted variance ratios alongside the standardized mean differences
  (Stata's `latebalance summarize`).
* `plot()` gains `type = "balance_density"` (covariate kernel densities by
  instrument arm, raw versus weighted; Stata's `latebalance density`) and a
  `geom = "density"` option for `type = "overlap"` (Stata's `lateoverlap`).

These diagnostics are internally validated; direct numerical parity with the
`lateffects` postestimation commands is planned for a later release.

# drlate 0.2.0

* New kappa-weighting estimators of the LATE from Słoczyński, Uysal, and
  Wooldridge (2025, *JBES* 43(1), 164-177): `method = "kappa"` (tau_a),
  `"kappa0"` (tau_a,0), and `"kappa10"` (tau_a,10), validated against the
  Stata `kappalate` command. Cluster-robust SEs, sampling weights, the
  bootstrap, and (for `"kappa"`/`"kappa0"`) Fieller confidence sets carry
  over from the existing machinery.
* Printed output now shows the kappalate names for the IPW estimators:
  normalized IPW is `tau_u`, unnormalized IPW is `tau_a,1`.
* New `ivmodel = "probit"` (kappalate's `zmodel(probit)`) for the
  weighting estimators (`"ipw"` and the kappa methods), completing
  coverage of the kappalate command's options.
* `drlate_compare()` now reports each kappa estimator's own normalization
  in the `normalized` column, and `?drlate` documents that the kappa
  denominators are kappa-weight means — estimating the same complier share
  as the IPW first-stage contrast, but as a different sample statistic.

# drlate 0.1.0 (patch)

Changes from an internal econometric audit (Monte Carlo evidence in
`data-raw/mc-review.R` and `data-raw/mc-weak2.R`):

* **Correctness (deliberate divergence from Stata 1.0.0):** the
  unnormalized LATT-AIPW estimator now computes the treated share with
  sampling weights. Stata uses an unweighted mean there, which leaves its
  `w1` moment condition nonzero under pweights and invalidates the joint
  variance; with uniform weights the two coincide exactly (all validated
  configurations are unaffected).
* Fieller sets handle the degenerate-quadratic regimes explicitly
  (half-line when `denom^2 = q * V_dd`, single-point tangency) instead of
  collapsing them into "whole line", and the complement-set print states
  that the set is unbounded.
* The weak-instrument advisory now triggers at first-stage `F < 10`
  (`|z| < 3.16`) instead of `|z| < 2`, and the printout reports
  `z^2 ~ F`.
* Cluster bootstrap warns when there are fewer than 30 clusters;
  documentation notes that failed bootstrap draws concentrate where
  identification is weak (prefer the Fieller set in that regime) and that
  seeded results are reproducible per fixed number of cores.
* `drlate_compare()` documents that IPW/RA rows use reduced adjustment
  sets (estimator and specification change together) and de-duplicates
  rows after normalization auto-switching.
* `dr_hausman()` documents that the LATT and ATT halves adjust on the
  instrument- and treatment-equation covariates respectively.

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
