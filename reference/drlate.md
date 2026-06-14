# Doubly robust estimation of the LATE and LATT

Estimates the local average treatment effect (LATE) or the local average
treatment effect on the treated (LATT) with a binary instrument,
following Słoczyński, Uysal, and Wooldridge (2022). A faithful R port of
the Stata package `drlate` (SSC S459708): point estimates come from
sequential weighted regressions, and standard errors are computed
jointly for the instrument propensity score, the outcome regression, the
treatment regression, and the causal estimand by stacking all moment
conditions into a single M-estimation system.

## Usage

``` r
drlate(
  outcome,
  treatment,
  instrument,
  data,
  omodel = c("linear", "logit", "probit", "poisson", "flogit", "fprobit"),
  tmodel = c("logit", "probit", "linear", "poisson"),
  ivmodel = c("logit", "cbps", "ipt", "probit"),
  method = c("ipwra", "ipw", "aipw", "ra", "kappa", "kappa0", "kappa10"),
  estimand = c("late", "latt"),
  normalized = TRUE,
  weights = NULL,
  cluster = NULL,
  pstolerance = 1e-05,
  osample = FALSE,
  subset = NULL,
  keep_data = TRUE,
  vcov = c("analytic", "bootstrap"),
  boot_reps = 999L,
  boot_seed = NULL,
  cores = 1L
)
```

## Arguments

- outcome:

  A formula `y ~ covariates` for the outcome model. Use `y ~ 1` for no
  covariates (required when `method = "ipw"`).

- treatment:

  A formula `d ~ covariates` for the treatment model.

- instrument:

  A formula `z ~ covariates` for the instrument propensity score model;
  `z` must be binary 0/1. Use `z ~ 1` when `method = "ra"`.

- data:

  A data frame containing all variables.

- omodel:

  Outcome model family: `"linear"` (default; continuous), `"logit"` or
  `"probit"` (outcome must be 0/1), `"poisson"` (outcome must be
  non-negative), or `"flogit"` / `"fprobit"` (fractional outcome in
  `[0, 1]`, e.g. a proportion). The `f`-prefixed families share all
  estimation with `"logit"` / `"probit"` and only relax the response to
  the unit interval, matching the Stata `lateffects` `omodel` options.

- tmodel:

  Treatment model family: `"logit"` (default; treatment must be 0/1),
  `"probit"`, `"linear"`, or `"poisson"`.

- ivmodel:

  Instrument propensity score model: `"logit"` (maximum likelihood;
  default), `"cbps"` (covariate balancing, Imai and Ratkovic 2014; not
  available with `estimand = "latt"`), `"ipt"` (inverse probability
  tilting, Graham, Pinto, and Egel 2012), or `"probit"` (maximum
  likelihood; mirrors kappalate's `zmodel(probit)` and is available only
  for the weighting estimators that command covers — `"ipw"`, `"kappa"`,
  `"kappa0"`, `"kappa10"` — with `estimand = "late"`).

- method:

  Estimator: `"ipwra"` (inverse-probability-weighted regression
  adjustment; default), `"ipw"`, `"aipw"`, `"ra"`, or one of the
  kappa-weighting estimators of Słoczyński, Uysal, and Wooldridge
  (2025): `"kappa"` (unnormalized Abadie kappa; kappalate's `tau_a`),
  `"kappa0"` (untreated-arm kappa; `tau_a,0`), or `"kappa10"`
  (normalized kappa; `tau_a,10`). The kappa estimators require
  intercept-only outcome and treatment formulas, a binary treatment, and
  `estimand = "late"`; `ivmodel = "cbps"` is available for `"kappa"`
  only, and `"ipt"` for none of them. drlate's normalized and
  unnormalized `"ipw"` coincide with kappalate's `tau_u` and `tau_a,1`.

- estimand:

  `"late"` (default) or `"latt"`.

- normalized:

  Logical; use normalized moment conditions (default `TRUE`). Only
  relevant for `method = "ipw"` and `method = "aipw"`.

- weights:

  Optional sampling weights (a numeric vector, or a column name in
  `data` given as a string).

- cluster:

  Optional cluster identifier for clustered standard errors (a vector,
  or a column name in `data` given as a string).

- pstolerance:

  Overlap tolerance: estimation stops with an error if any estimated
  instrument propensity score is below `pstolerance` or above
  `1 - pstolerance`. Default `1e-5`.

- osample:

  Logical; if `TRUE`, overlap violations do not stop estimation with an
  error. Instead `drlate()` returns (invisibly) a logical vector marking
  the violating observations.

- subset:

  Optional logical or integer vector selecting rows of `data`.

- keep_data:

  Logical; retain the internal estimation context (model matrices,
  fitted propensity scores, weights) on the returned object (default
  `TRUE`). Required by
  [`plot.drlate()`](https://kvenkita.github.io/drlate/reference/plot.drlate.md),
  [`balance()`](https://kvenkita.github.io/drlate/reference/balance.md),
  and the bootstrap; set to `FALSE` for a leaner object.

- vcov:

  `"analytic"` (default) for the joint M-estimation sandwich replicating
  the Stata package, or `"bootstrap"` for nonparametric bootstrap
  standard errors and percentile confidence intervals (whole clusters
  are resampled when `cluster` is supplied). The analytic variance is
  always computed and stored either way. Draws that fail (degenerate
  resamples, non-convergence, overlap violations) are dropped and
  counted; because such failures concentrate where identification is
  weak, a non-trivial failure rate is itself a sign that percentile
  intervals are unreliable and the Fieller set
  (`confint(., method = "fieller")`) should be preferred.

- boot_reps:

  Number of bootstrap replications (default 999).

- boot_seed:

  Optional seed for reproducible bootstrap draws. Results are
  reproducible for a fixed number of `cores`; serial and parallel runs
  use different (both valid) random streams.

- cores:

  Number of CPU cores for the bootstrap (default 1). Values above 1 use
  a PSOCK cluster and require the package to be installed (not merely
  loaded with `devtools::load_all()`).

## Value

An object of class `"drlate"`, a list with components including
`coefficients` (the causal estimate, the numerator effect of Z on Y, and
the denominator effect of Z on D), `vcov3` (their variance matrix,
diagonal by construction, as in the Stata package), `vcov_full` (the
joint variance matrix of all stacked parameters), `theta` (all stacked
parameter estimates), `N`, `dmeanz1`, `dmeanz0`, and the call. For
`method = "kappa10"` only the causal estimate is reported (the estimator
is a difference of two ratios, so no single numerator/denominator pair
exists). For `"kappa"` and `"kappa0"` the third coefficient is the mean
of the corresponding kappa weight: under the LATE assumptions it
estimates the same complier share as the IPW first-stage contrast (the
population ATE of Z on D), but it is a different sample statistic and
the two can diverge under propensity score misspecification.

## References

Słoczyński, T., S. D. Uysal, and J. M. Wooldridge (2022). "Doubly Robust
Estimation of Local Average Treatment Effects Using Inverse Probability
Weighted Regression Adjustment."
[doi:10.48550/arXiv.2208.01300](https://doi.org/10.48550/arXiv.2208.01300)

Słoczyński, T., S. D. Uysal, and J. M. Wooldridge (2025). "Abadie's
Kappa and Weighting Estimators of the Local Average Treatment Effect."
*Journal of Business & Economic Statistics* 43(1), 164–177.
[doi:10.1080/07350015.2024.2332763](https://doi.org/10.1080/07350015.2024.2332763)

## Examples

``` r
data(drlate_sim)
fit <- drlate(lwage ~ age + educ, nvstat ~ age + educ,
              rsncode ~ age + educ, data = drlate_sim)
summary(fit)
#> 
#> Local average treatment effect
#> Number of obs    : 2,000
#> Estimator        : IPWRA
#> Outcome model    : linear
#> Treatment model  : logit
#> Instrument model : logit (MLE)
#> 
#>              Estimate Std. Error z value   Pr(>|z|) [95% conf. interval]
#> LATE: D on Y   0.4705    0.07915   5.944  2.786e-09     0.3153    0.6256
#> ATE: Z on Y    0.2845    0.05043   5.642  1.679e-08     0.1857    0.3834
#> ATE: Z on D    0.6048    0.01837  32.929 8.326e-238     0.5688    0.6408
#> 
#> First stage (Z on D): z = 32.93 (z^2 ~ first-stage F = 1084)
```
