# Compare drlate estimators in one call

Runs several estimators on the same specification and collects the
causal estimates with their confidence intervals — the sensitivity
comparison applied papers routinely report. Formula restrictions are
handled automatically: `method = "ipw"` drops the outcome/treatment
covariates and `method = "ra"` drops the instrument covariates (each
with a message), matching the requirements of those estimators.

## Usage

``` r
drlate_compare(
  outcome,
  treatment,
  instrument,
  data,
  methods = c("ipwra", "ipw", "aipw", "ra"),
  both_norms = FALSE,
  ...
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

- methods:

  Estimators to run (any of the `method` values accepted by
  [`drlate()`](https://kvenkita.github.io/drlate/reference/drlate.md)).

- both_norms:

  Logical; also run the unnormalized variants of `"ipw"` and `"aipw"`
  (default `FALSE`).

- ...:

  Passed on to
  [`drlate()`](https://kvenkita.github.io/drlate/reference/drlate.md)
  (e.g. `omodel`, `tmodel`, `ivmodel`, `estimand`, `weights`,
  `cluster`).

## Value

An object of class `"drlate_compare"`: a data frame with columns
`method`, `normalized`, `estimate`, `se`, `ci_lo`, `ci_hi`, with a
`print` method and a dot-whisker `plot` method.

## Details

Because IPW carries no outcome/treatment regressions and RA carries no
instrument propensity score, the automatic formula adjustment means the
rows do not share a single adjustment specification: differences between
the IPW or RA row and the doubly robust rows reflect both the estimator
*and* the reduced specification. Read the comparison as a robustness
display, not as a test that isolates estimator choice; the doubly robust
rows (IPWRA, AIPW) are the like-for-like pair.

## Examples

``` r
cmp <- drlate_compare(lwage ~ age + educ, nvstat ~ age + educ,
                      rsncode ~ age + educ, data = drlate_sim)
#> method = "ipw": dropping outcome/treatment covariates (weighted means only).
#> method = "ra": dropping instrument covariates (no propensity score).
cmp
#> Estimator comparison (LATE)
#> 
#>   estimator estimate     se           95% CI
#>       ipwra   0.4705 0.0792 [0.3153, 0.6256]
#>   ipw (nrm)   0.4741 0.0793 [0.3187, 0.6295]
#>  aipw (nrm)   0.4702 0.0792 [0.3150, 0.6254]
#>          ra   0.4597 0.0792 [0.3045, 0.6150]
```
