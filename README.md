# drlate

`drlate` estimates the local average treatment effect (LATE) and the local
average treatment effect on the treated (LATT) from observational data with
a binary instrument, following
[Słoczyński, Uysal & Wooldridge (2022)](https://doi.org/10.48550/arXiv.2208.01300).
It is an R port of the Stata package
[`drlate`](https://ideas.repec.org/c/boc/bocode/s459708.html) (SSC S459708)
by S. Derya Uysal, Tymon Słoczyński, and Jeffrey M. Wooldridge: point
estimates come from the same sequential weighted regressions, standard
errors from the same jointly stacked M-estimation system and its sandwich
variance, and the test suite verifies numerical equivalence against
Stata-generated fixtures across 33 scenarios. On top of the port, the
package provides design diagnostics, weak-instrument-robust and bootstrap
inference, and the paper's doubly robust Hausman test of unconfoundedness.

## Installation

```r
# install.packages("remotes")
remotes::install_github("kvenkita/drlate")
```

## Example

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
#> LATE: D on Y   0.4705    0.07915   5.944  2.786e-09     0.3153    0.6256
#> ATE: Z on Y    0.2845    0.05043   5.642  1.679e-08     0.1857    0.3834
#> ATE: Z on D    0.6048    0.01837  32.929 8.326e-238     0.5688    0.6408
#>
#> First stage (Z on D): z = 32.93 (z^2 ~ first-stage F = 1084)
```

## Citation

If you use drlate in your research, please cite:

> Venkitasubramanian, K. (2026). drlate: Doubly Robust Estimation
> of the Local Average Treatment Effect in R. R package version 0.1.0.
> https://github.com/kvenkita/drlate

This package implements the doubly-robust LATE estimator introduced in:

> Uysal, D., Słoczyński, T., & Wooldridge, J. M. (2026). DRLATE: Stata
> module to perform doubly robust estimation of the local average
> treatment effect (LATE) and the local average treatment effect on the
> treated (LATT). Statistical Software Components S459708, Boston College
> Department of Economics.

(`citation("drlate")` prints both entries with BibTeX.)

## Features

| | |
|---|---|
| Estimands | LATE, LATT |
| Estimators (`method`) | IPWRA (default), IPW, AIPW, RA |
| Outcome / treatment models | linear, logit, Poisson |
| Instrument propensity score models (`ivmodel`) | logit MLE (default), CBPS, IPT |
| Weighting | normalized (default) or unnormalized moments; sampling weights |
| Standard errors | joint sandwich over all estimation stages; robust or cluster-robust |
| Diagnostics | `plot(fit)` for propensity-score overlap, covariate balance, and weight distributions; `balance()` tables; first-stage strength on every printout |
| Fieller confidence sets | `confint(fit, method = "fieller")` — weak-instrument-robust |
| Bootstrap | `drlate(..., vcov = "bootstrap")` — cluster-aware percentile intervals |
| DR Hausman test | `dr_hausman()` — test of unconfoundedness under one-sided noncompliance (paper, Section 5) |
| Estimator comparison | `drlate_compare()` with a dot-whisker plot |
| Overlap | `pstolerance` enforcement, `osample` violator flagging |

## Documentation

The [package website](https://kvenkita.github.io/drlate/) serves the
[primer](https://kvenkita.github.io/drlate/articles/drlate-primer.html),
the [package overview and Stata replication](https://kvenkita.github.io/drlate/articles/drlate.html),
and the function reference.

## License

MIT. Portions derived from the Stata package `drlate`, © 2026 S. Derya
Uysal, Tymon Słoczyński, and Jeffrey M. Wooldridge (MIT licensed).
