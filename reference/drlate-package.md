# drlate: Doubly Robust Estimation of Local Average Treatment Effects

Estimates the local average treatment effect (LATE) and the local
average treatment effect on the treated (LATT) using observational data
with a binary instrument, following Sloczynski, Uysal, and Wooldridge
(2022)
[doi:10.48550/arXiv.2208.01300](https://doi.org/10.48550/arXiv.2208.01300)
. Supports inverse probability weighted regression adjustment (IPWRA),
inverse probability weighting (IPW), augmented inverse probability
weighting (AIPW), and regression adjustment (RA), with linear, logistic,
or Poisson outcome and treatment models, and instrument propensity
scores estimated by maximum likelihood, covariate balancing (CBPS), or
inverse probability tilting (IPT). Standard errors are computed jointly
for all estimation stages by stacking the moment conditions of every
model into a single M-estimation system. A faithful R port of the Stata
package 'drlate' (Statistical Software Components S459708).

## See also

Useful links:

- <https://github.com/kvenkita/drlate>

- <https://kvenkita.github.io/drlate/>

- Report bugs at <https://github.com/kvenkita/drlate/issues>

## Author

**Maintainer**: Kailas Venkitasubramanian <kailasv@gmail.com>

Other contributors:

- S. Derya Uysal (Author of the original Stata package 'drlate')
  \[contributor, copyright holder\]

- Tymon Sloczynski (Author of the original Stata package 'drlate')
  \[contributor, copyright holder\]

- Jeffrey M. Wooldridge (Author of the original Stata package 'drlate')
  \[contributor, copyright holder\]
