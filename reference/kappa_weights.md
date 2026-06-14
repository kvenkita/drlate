# Abadie's kappa weights

Returns the per-observation Abadie kappa weight implied by a fitted
[`drlate()`](https://kvenkita.github.io/drlate/reference/drlate.md)
object, \$\$\kappa = 1 - \frac{D(1 - Z)}{1 - p(X)} - \frac{(1 - D)
Z}{p(X)},\$\$ where \\p(X)\\ is the estimated instrument propensity
score. The kappa weights identify the complier subpopulation: for any
function \\g\\ of the data, \\E\[g \mid \mathrm{complier}\] = E\[\kappa
g\] / E\[\kappa\]\\ (Abadie 2003). They are the weights used by
[`complier_means()`](https://kvenkita.github.io/drlate/reference/complier_means.md)
and are the Stata `estat compliers, genkappa()` object.

## Usage

``` r
kappa_weights(object, normalize = TRUE)
```

## Arguments

- object:

  A fitted
  [`drlate()`](https://kvenkita.github.io/drlate/reference/drlate.md)
  object (with `keep_data = TRUE`) using an instrument propensity score
  (any `method` except `"ra"`).

- normalize:

  Logical. If `TRUE` (default), the returned weights are the
  sampling-weighted, normalized weights \\w\kappa / \sum w\kappa\\ that
  sum to one (the form used to compute complier averages). If `FALSE`,
  the raw kappa values are returned.

## Value

A numeric vector with one entry per estimation-sample observation.

## See also

[`complier_means()`](https://kvenkita.github.io/drlate/reference/complier_means.md)
