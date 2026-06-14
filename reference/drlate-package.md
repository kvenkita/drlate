# drlate: Doubly Robust Estimation of Local Average Treatment Effects

Estimates the local average treatment effect (LATE) and the local
average treatment effect on the treated (LATT) using observational data
with a binary instrument, implementing the complete estimator suite of
Sloczynski, Uysal, and Wooldridge: the doubly robust estimators of
Sloczynski, Uysal, and Wooldridge (2022)
[doi:10.48550/arXiv.2208.01300](https://doi.org/10.48550/arXiv.2208.01300)
– inverse probability weighted regression adjustment (IPWRA), inverse
probability weighting (IPW), augmented inverse probability weighting
(AIPW), and regression adjustment (RA) – and the Abadie-kappa weighting
estimators of Sloczynski, Uysal, and Wooldridge (2025)
[doi:10.1080/07350015.2024.2332763](https://doi.org/10.1080/07350015.2024.2332763)
. Supports linear, logistic, probit, Poisson, and fractional
(fractional-logit and fractional-probit) outcome and treatment models,
and instrument propensity scores estimated by maximum likelihood,
covariate balancing (CBPS), or inverse probability tilting (IPT).
Standard errors are computed jointly for all estimation stages by
stacking the moment conditions of every model into a single M-estimation
system; weak-instrument-robust Fieller confidence sets, cluster-aware
bootstrap inference, design diagnostics, and a doubly robust
Hausman-type test of unconfoundedness are included. Estimates and
standard errors are validated against the authors' Stata commands
'drlate' (Statistical Software Components S459708) and 'kappalate'
(S459257).

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
