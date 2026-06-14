# Complier covariate means

Compares the average of each covariate in the full estimation sample
with its average in the complier subpopulation, the latter computed with
the normalized Abadie kappa weights of
[`kappa_weights()`](https://kvenkita.github.io/drlate/reference/kappa_weights.md).
Because the local average treatment effect is a causal effect for
compliers, knowing how compliers differ from the population aids
interpretation. This is the Stata `estat compliers` postestimation
feature.

## Usage

``` r
complier_means(object, vars = NULL)
```

## Arguments

- object:

  A fitted
  [`drlate()`](https://kvenkita.github.io/drlate/reference/drlate.md)
  object (with `keep_data = TRUE`) using an instrument propensity score
  (any `method` except `"ra"`).

- vars:

  Optional character vector selecting a subset of the model covariates.
  Defaults to all covariates across the three model formulas.

## Value

A data frame with one row per covariate and columns `variable`,
`population_mean`, `complier_mean`, and `difference`
(`complier_mean - population_mean`).

## Details

Covariate values are reported on their original scale.

## See also

[`kappa_weights()`](https://kvenkita.github.io/drlate/reference/kappa_weights.md)

## Examples

``` r
fit <- drlate(lwage ~ age + educ, nvstat ~ age + educ,
              rsncode ~ age + educ, data = drlate_sim)
complier_means(fit)
#>       variable population_mean complier_mean   difference
#> 1          age         34.5560    34.3303393 -0.225660695
#> 2  educcollege          0.3615     0.3590211 -0.002478943
#> 3 educgraduate          0.1395     0.1431700  0.003670025
```
