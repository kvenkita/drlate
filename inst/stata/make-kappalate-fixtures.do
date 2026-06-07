*! Golden-fixture generator for the drlate kappa-weighting estimators.
*!
*! Run from the R package root directory (the folder containing DESCRIPTION):
*!     cd <path-to-drlate-R-package>
*!     do inst/stata/make-kappalate-fixtures.do
*!
*! Requires the Stata kappalate package:  ssc install kappalate
*! Writes one CSV per scenario into tests/testthat/fixtures/.

version 17
clear all
set more off

cap mkdir "tests/testthat/fixtures"

* --------------------------------------------------------------------------
* Helper: write e(b)/e(V) of the last kappalate call to a CSV fixture.
* Column j of e(b) is written as bj/vj; the R tests map positions to
* estimators per scenario (which(all) + zmodel(logit):
* 1=tau_a, 2=tau_a,1, 3=tau_a,0, 4=tau_a,10, 5=tau_u;
* which(all) + zmodel(cbps): 1=tau_a, 2=tau_u).
* --------------------------------------------------------------------------
capture program drop kapout
program define kapout
    args id
    tempname f b V
    matrix `b' = e(b)
    matrix `V' = e(V)
    local k = colsof(`b')
    file open `f' using "tests/testthat/fixtures/`id'.csv", write replace
    file write `f' "name,value" _n
    file write `f' "N," %21.0g (e(N)) _n
    file write `f' "k," %21.0g (`k') _n
    forvalues j = 1/`k' {
        file write `f' "b`j'," %21.0g (el(`b',1,`j')) _n
        file write `f' "v`j'," %21.0g (el(`V',`j',`j')) _n
    }
    file close `f'
    di as txt "wrote fixture: `id'"
end

capture program drop runkap
program define runkap
    gettoken id 0 : 0
    cap noi `0'
    if _rc {
        di as err "SCENARIO FAILED (`id'): rc = " _rc
    }
    else {
        kapout `id'
    }
end

* --------------------------------------------------------------------------
* Data preparation: identical to inst/stata/make-fixtures.do
* --------------------------------------------------------------------------
capture confirm file "tests/testthat/fixtures/sipp.dta"
if !_rc {
    use "tests/testthat/fixtures/sipp.dta", clear
}
else {
    use "https://people.brandeis.edu/~tslocz/sipp.dta", clear
}
drop if kwage == . | educ == . | rsncode == 999
generate double lwage = ln(kwage)
generate cluvar = educ

* --------------------------------------------------------------------------
* Scenarios
* --------------------------------------------------------------------------
runkap kappalate_logit_all  kappalate lwage (nvstat = rsncode) age_5, zmodel(logit) which(all)
runkap kappalate_cbps_all   kappalate lwage (nvstat = rsncode) age_5, zmodel(cbps) which(all)
runkap kappalate_logit_clu  kappalate lwage (nvstat = rsncode) age_5, zmodel(logit) which(all) vce(cluster cluvar)

* One-sided noncompliance: force D = 0 whenever Z = 0
preserve
replace nvstat = 0 if rsncode == 0
runkap kappalate_logit_onesided kappalate lwage (nvstat = rsncode) age_5, zmodel(logit) which(all)
restore

di as txt "kappalate fixtures done."
