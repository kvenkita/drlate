# Confidence intervals for drlate fits

Confidence intervals for drlate fits

## Usage

``` r
# S3 method for class 'drlate'
confint(object, parm, level = 0.95, method = c("default", "fieller"), ...)
```

## Arguments

- object:

  A fitted
  [`drlate()`](https://kvenkita.github.io/drlate/reference/drlate.md)
  object.

- parm:

  Coefficients to include (names or indices); defaults to all three
  reported quantities.

- level:

  Confidence level.

- method:

  `"default"` gives Wald intervals from the joint sandwich (or bootstrap
  percentile intervals when the fit used `vcov = "bootstrap"`).
  `"fieller"` inverts the test of `num - t * denom = 0` using the joint
  covariance of the numerator and denominator, giving a confidence set
  for the LATE/LATT ratio that remains valid when the first stage is
  weak; the set may be an interval, the complement of an interval, or
  the whole line, and is returned as a `"drlate_fieller"` object with
  its own print method.

- ...:

  Currently unused.

## Value

For `method = "default"`, a numeric matrix with one row per requested
coefficient (`parm`) and two columns holding the lower and upper
confidence limits. The columns are labelled with the corresponding
percentiles (for the default 95% level, `"2.5 %"` and `"97.5 %"`). The
limits are Wald intervals from the joint sandwich covariance, or
percentile intervals from the resampling draws when the fit was computed
with `vcov = "bootstrap"`.

For `method = "fieller"`, an object of class `"drlate_fieller"`: a list
describing the weak-instrument-robust confidence set for the LATE/LATT
ratio (its endpoints and shape, the estimand name, and the confidence
level), with its own `print` method. Because a Fieller set need not be a
bounded interval, it is returned in this form rather than as a matrix of
endpoints.
