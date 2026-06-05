# Simulated example data for drlate

A simulated dataset with a binary instrument, a binary treatment with
two-sided noncompliance, and continuous, positive, and binary outcome
variables, designed to exercise every model family supported by
[`drlate()`](https://kvenkita.github.io/drlate/reference/drlate.md). The
complier average treatment effect (LATE) used in the data-generating
process is 0.5. The treatment is genuinely endogenous (compliance type
shifts the baseline outcome, so naive OLS is biased upward) and the
instrument is only conditionally valid (its propensity depends on `age`
and `educ`, so the raw Wald ratio is biased too).

## Usage

``` r
drlate_sim
```

## Format

A data frame with 2,000 rows and 7 variables:

- lwage:

  continuous outcome

- kwage:

  positive outcome (for Poisson models), `exp(lwage / 2)`

- hijob:

  binary outcome (for logit models)

- nvstat:

  binary treatment

- rsncode:

  binary instrument

- age:

  continuous covariate

- educ:

  factor covariate with levels `hs`, `college`, `graduate`

## Source

Simulated; see `data-raw/drlate_sim.R` in the package sources.
