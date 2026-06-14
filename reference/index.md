# Package index

## Estimation

Fit doubly robust LATE/LATT estimators.

- [`drlate()`](https://kvenkita.github.io/drlate/reference/drlate.md) :
  Doubly robust estimation of the LATE and LATT
- [`drlate_compare()`](https://kvenkita.github.io/drlate/reference/drlate_compare.md)
  : Compare drlate estimators in one call

## Diagnostics

Check overlap, balance, weights, and instrument strength.

- [`plot(`*`<drlate>`*`)`](https://kvenkita.github.io/drlate/reference/plot.drlate.md)
  : Diagnostic plots for drlate fits
- [`balance()`](https://kvenkita.github.io/drlate/reference/balance.md)
  : Covariate balance across instrument arms
- [`balance_test()`](https://kvenkita.github.io/drlate/reference/balance_test.md)
  : Imai-Ratkovic covariate-balance test

## Complier profiling

Characterize the complier subpopulation via Abadie’s kappa.

- [`complier_means()`](https://kvenkita.github.io/drlate/reference/complier_means.md)
  : Complier covariate means
- [`kappa_weights()`](https://kvenkita.github.io/drlate/reference/kappa_weights.md)
  : Abadie's kappa weights

## Inference

Confidence intervals, including weak-instrument-robust sets.

- [`confint(`*`<drlate>`*`)`](https://kvenkita.github.io/drlate/reference/confint.drlate.md)
  : Confidence intervals for drlate fits

## Specification tests

- [`dr_hausman()`](https://kvenkita.github.io/drlate/reference/dr_hausman.md)
  : Doubly robust Hausman test of unconfoundedness

## Data

- [`drlate_sim`](https://kvenkita.github.io/drlate/reference/drlate_sim.md)
  : Simulated example data for drlate

## Package

- [`drlate-package`](https://kvenkita.github.io/drlate/reference/drlate-package.md)
  : drlate: Doubly Robust Estimation of Local Average Treatment Effects
