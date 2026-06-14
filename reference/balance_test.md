# Imai-Ratkovic covariate-balance test

Tests whether the estimated instrument propensity score balances the
covariates, using the overidentification test of Imai and Ratkovic
(2014). The propensity-score MLE score equations identify the
coefficients; the covariate-balancing (CBPS) moments are the
overidentifying restrictions. A large statistic is evidence that the
propensity-score model does not balance the covariates — a
misspecification diagnostic. This is the Stata `latebalance overid`
postestimation feature.

## Usage

``` r
balance_test(object)
```

## Arguments

- object:

  A fitted
  [`drlate()`](https://kvenkita.github.io/drlate/reference/drlate.md)
  object (with `keep_data = TRUE`) using a logistic or probit instrument
  propensity score.

## Value

An object of class `drlate_balance_test`: a list with `statistic`
(Hansen's J), `df`, `p.value`, `ivmodel`, and `n`, with a `print`
method.

## References

Imai, K. and Ratkovic, M. (2014). Covariate Balancing Propensity Score.
*Journal of the Royal Statistical Society B* 76(1), 243–263.

## See also

[`balance()`](https://kvenkita.github.io/drlate/reference/balance.md)
for the standardized-mean-difference diagnostics.

## Examples

``` r
fit <- drlate(lwage ~ age + educ, nvstat ~ age + educ,
              rsncode ~ age + educ, data = drlate_sim)
balance_test(fit)
#> Imai-Ratkovic covariate-balance test (overidentification)
#> 
#>   Hansen J = 3.0473   df = 4   p-value = 0.5499
#>   Instrument propensity score: logit (n = 2000)
#> 
#>   H0: the propensity-score model balances the covariates.
```
