# drlate

Doubly robust estimation of the local average treatment effect (LATE) and
the local average treatment effect on the treated (LATT) in R, following
[Słoczyński, Uysal & Wooldridge (2022)](https://doi.org/10.48550/arXiv.2208.01300).

The estimation core is a faithful R port of the Stata package
[`drlate`](https://ideas.repec.org/c/boc/bocode/s459708.html) (SSC S459708)
by S. Derya Uysal, Tymon Słoczyński, and Jeffrey M. Wooldridge: point
estimates come from the same sequential weighted regressions, and standard
errors are computed jointly for the instrument propensity score, the outcome
regression, the treatment regression, and the causal estimand by stacking
all moment conditions into a single M-estimation system and computing its
sandwich variance. This is the R equivalent of the original's
`gmm, onestep iterate(0)` construction. The test suite
verifies numerical equivalence against Stata-generated fixtures (estimates
to ~1e-9, standard errors to ~1e-6, across 33 scenarios).

On top of the port, the package adds diagnostics and inference tools that
the Stata original does not provide, each validated by Monte Carlo (see the
[validation report](https://kvenkita.github.io/drlate/articles/validation-review.html)).

## Features

### The estimation core (Stata-equivalent)

| | |
|---|---|
| Estimands | LATE, LATT |
| Estimators (`method`) | IPWRA (default), IPW, AIPW, RA |
| Outcome / treatment models | linear, logit, Poisson |
| Instrument PS models (`ivmodel`) | logit MLE (default), CBPS, IPT |
| Weighting | normalized (default) or unnormalized moments |
| Inference | robust or cluster-robust joint sandwich SEs; sampling weights |
| Overlap | `pstolerance` enforcement, `osample` violator flagging |

### Beyond the Stata package

| | |
|---|---|
| Diagnostics | `plot(fit)`: propensity-score **overlap**, covariate **balance** (love plot), implied **weight** distributions; `balance()` returns weighted/unweighted standardized mean differences |
| First-stage strength | every printout reports the first-stage z and z² ≈ F, flagging weakness below F = 10 |
| Weak-IV-robust inference | `confint(fit, method = "fieller")` — Fieller/Anderson–Rubin confidence sets with honest unbounded regimes (Monte Carlo coverage 0.955 where the Wald interval degenerates) |
| Bootstrap | `drlate(..., vcov = "bootstrap")` — paper-recommended nonparametric bootstrap; cluster-aware, parallel, failure-counting |
| Specification test | `dr_hausman()` — the doubly robust Hausman test of unconfoundedness from the paper's Section 5 (**not implemented in the Stata package**), with an analytic SE from one jointly stacked moment system; Monte Carlo size 0.047 at nominal 5% |
| Sensitivity | `drlate_compare()` — IPWRA/IPW/AIPW/RA side by side with a dot-whisker plot |

## Installation

```r
# install.packages("remotes")
remotes::install_github("kvenkita/drlate")
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
#> LATE: D on Y   0.4705    0.07915   5.944  2.786e-09     0.3153    0.6256
#> ATE: Z on Y    0.2845    0.05043   5.642  1.679e-08     0.1857    0.3834
#> ATE: Z on D    0.6048    0.01837  32.929 8.326e-238     0.5688    0.6408
#>
#> First stage (Z on D): z = 32.93 (z^2 ~ first-stage F = 1084)
```

Diagnostics, robustness, and the unconfoundedness test:

```r
plot(fit, "overlap")                      # propensity-score overlap by arm
plot(fit, "balance")                      # love plot, unweighted vs IPW
plot(fit, "weights")                      # implied weight distributions

confint(fit, method = "fieller")          # weak-IV-robust confidence set
drlate(..., vcov = "bootstrap")           # bootstrap SEs / percentile CIs

cmp <- drlate_compare(lwage ~ age + educ, nvstat ~ age + educ,
                      rsncode ~ age + educ, data = drlate_sim)
plot(cmp)                                 # estimator sensitivity

# DR Hausman test of unconfoundedness (one-sided noncompliance)
d <- drlate_sim; d$nvstat[d$rsncode == 0] <- 0L
dr_hausman(lwage ~ age + educ, nvstat ~ age + educ,
           rsncode ~ age + educ, data = d)
```

## Documentation

- [**Primer**](https://kvenkita.github.io/drlate/articles/drlate-primer.html) —
  doubly robust LATE estimation from intuition to practice, with worked
  examples of every feature
- [**Package overview and Stata replication**](https://kvenkita.github.io/drlate/articles/drlate.html)
- [**Validation report and peer review**](https://kvenkita.github.io/drlate/articles/validation-review.html) —
  port equivalence, theory conformance, Monte Carlo evidence

## Validating against Stata

`inst/stata/make-fixtures.do` runs the original Stata `drlate` over a grid
of scenarios on the public SIPP extract and exports estimates and variances
as CSV fixtures; the testthat suite then asserts equality (estimates to
1e-6, standard errors to 1e-4 relative; observed agreement is two to three
orders of magnitude tighter). Run it from the package root in Stata
(`ssc install drlate` first), then `devtools::test()`.

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

## License

MIT. Portions derived from the Stata package `drlate`, © 2026 S. Derya
Uysal, Tymon Słoczyński, and Jeffrey M. Wooldridge (MIT licensed).
