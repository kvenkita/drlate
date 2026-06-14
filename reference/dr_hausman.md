# Doubly robust Hausman test of unconfoundedness

Tests whether the treatment is unconfounded given the covariates, using
the comparison proposed by Słoczyński, Uysal, and Wooldridge (2022,
Section 5), building on Donald, Hsu, and Lieli (2014). Under **one-sided
noncompliance** (nobody takes the treatment without the instrument:
\\\Pr(D = 1 \mid Z = 0) = 0\\), the LATT identified through the
instrument equals the ATT identified through unconfoundedness of the
treatment — so a significant difference between the doubly robust LATT
estimate (which uses the instrument) and the doubly robust ATT estimate
(which does not) is evidence against unconfoundedness. Unlike the
textbook OLS-vs-IV Hausman test, this comparison is robust to treatment
effect heterogeneity.

## Usage

``` r
dr_hausman(
  outcome,
  treatment,
  instrument,
  data,
  omodel = c("linear", "logit", "poisson"),
  tmodel = c("logit", "linear", "poisson"),
  ivmodel = c("logit", "ipt"),
  weights = NULL,
  cluster = NULL,
  pstolerance = 1e-05,
  subset = NULL
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

  Instrument propensity score model for the LATT half: `"logit"`
  (default) or `"ipt"`.

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

- subset:

  Optional logical or integer vector selecting rows of `data`.

## Value

An object of class `"htest"` with the z statistic, p-value, and the DR
LATT, DR ATT, and difference estimates.

## Details

The DR ATT estimator follows the paper's equation (33): a treatment
propensity score \\\Pr(D = 1 \mid X)\\ is fitted by logit QMLE on the
treatment-equation covariates; the outcome model is fitted on the
untreated sample weighted by the odds \\\hat p/(1-\hat p)\\; and
\\\hat\tau\_{ATT}\\ is the treated-sample mean outcome minus the mean
imputed counterfactual. The standard error of the difference comes from
stacking the moment conditions of *both* estimators (and the difference)
into one M-estimation system, so the covariance between them is
accounted for analytically — the analytic option suggested in the paper.

Note that the two halves adjust on their respective formulas: the LATT
half's propensity score uses the *instrument*-equation covariates, while
the ATT half's uses the *treatment*-equation covariates (both share the
outcome model). Supply the same covariate set to all three formulas
unless you intend them to differ.

## References

Słoczyński, T., S. D. Uysal, and J. M. Wooldridge (2022). "Doubly Robust
Estimation of Local Average Treatment Effects Using Inverse Probability
Weighted Regression Adjustment."
[doi:10.48550/arXiv.2208.01300](https://doi.org/10.48550/arXiv.2208.01300)

Donald, S. G., Y.-C. Hsu, and R. P. Lieli (2014). "Testing the
Unconfoundedness Assumption via Inverse Probability Weighted Estimators
of (L)ATT." *Journal of Business & Economic Statistics* 32(3), 395-415.

## Examples

``` r
d <- drlate_sim
d$nvstat[d$rsncode == 0] <- 0L   # impose one-sided noncompliance
dr_hausman(lwage ~ age + educ, nvstat ~ age + educ,
           rsncode ~ age + educ, data = d)
#> 
#>  Doubly robust Hausman test of unconfoundedness
#>  (Sloczynski-Uysal-Wooldridge 2022, one-sided noncompliance)
#> 
#> data:  d
#> z = -5.7425, p-value = 9.331e-09
#> alternative hypothesis: two.sided
#> sample estimates:
#>    DR LATT     DR ATT difference 
#>  0.3760331  0.6323210 -0.2562878 
#> 
```
