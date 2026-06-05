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
