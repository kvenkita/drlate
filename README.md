# drlate

Doubly robust estimation of the local average treatment effect (LATE) and
the local average treatment effect on the treated (LATT) in R, following
[Słoczyński, Uysal & Wooldridge (2022)](https://doi.org/10.48550/arXiv.2208.01300).

This is a faithful R port of the Stata package
[`drlate`](https://ideas.repec.org/c/boc/bocode/s459708.html) (SSC S459708)
by S. Derya Uysal, Tymon Słoczyński, and Jeffrey M. Wooldridge: point
estimates come from the same sequential weighted regressions, and standard
errors are computed jointly for the instrument propensity score, the outcome
regression, the treatment regression, and the causal estimand by stacking
all moment conditions into a single M-estimation sandwich — the R equivalent
of the original's `gmm, onestep iterate(0)` construction. The test suite
verifies numerical equivalence against Stata-generated fixtures.

## Features

| | |
|---|---|
| Estimands | LATE, LATT |
| Estimators (`method`) | IPWRA (default), IPW, AIPW, RA |
| Outcome / treatment models | linear, logit, Poisson |
| Instrument PS models (`ivmodel`) | logit MLE (default), CBPS, IPT |
| Weighting | normalized (default) or unnormalized moments |
| Inference | robust or cluster-robust joint sandwich SEs; sampling weights |

## Installation

```r
# install.packages("remotes")
remotes::install_github("kailasv/drlate")
```

## Usage

```r
library(drlate)
data(drlate_sim)

fit <- drlate(lwage ~ age + educ,    # outcome model
              nvstat ~ age + educ,   # treatment model
              rsncode ~ age + educ,  # instrument propensity score model
              data = drlate_sim)
summary(fit)
#> Local average treatment effect
#> Number of obs    : 2,000
#> Estimator        : IPWRA
#> Outcome model    : linear
#> Treatment model  : logit
#> Instrument model : logit (MLE)
#>
#>              Estimate Std. Error z value   Pr(>|z|) [95% conf. interval]
#> LATE: D on Y   0.4414    0.07632   5.783  7.343e-09     0.2918    0.5910
#> ATE: Z on Y    0.2669    0.04666   5.721  1.057e-08     0.1755    0.3584
#> ATE: Z on D    0.6048    0.01837  32.929 8.326e-238     0.5688    0.6408
```

See `vignette("drlate")` for the full methodology and more examples.

## Validating against Stata

`inst/stata/make-fixtures.do` runs the original Stata `drlate` over a grid
of scenarios on the public SIPP extract and exports estimates and variances
as CSV fixtures; the testthat suite then asserts equality (estimates to
1e-6, standard errors to 1e-4 relative). Run it from the package root in
Stata (`ssc install drlate` first), then `devtools::test()`.

## License

MIT. Portions derived from the Stata package `drlate`, © 2026 S. Derya
Uysal, Tymon Słoczyński, and Jeffrey M. Wooldridge (MIT licensed).
