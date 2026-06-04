*! Golden-fixture generator for the drlate R package.
*!
*! Run from the R package root directory (the folder containing DESCRIPTION):
*!     cd <path-to-drlate-R-package>
*!     do inst/stata/make-fixtures.do
*!
*! Requires the Stata drlate package:  ssc install drlate
*! Writes one CSV per scenario into tests/testthat/fixtures/.

version 17
clear all
set more off

cap mkdir "tests/testthat/fixtures"

* --------------------------------------------------------------------------
* Helper: write e() results of the last drlate call to a CSV fixture
* --------------------------------------------------------------------------
capture program drop fixout
program define fixout
    args id
    tempname f
    file open `f' using "tests/testthat/fixtures/`id'.csv", write replace
    file write `f' "name,value" _n
    file write `f' "N," %21.0g (e(N)) _n
    file write `f' "dmeanz1," %21.0g (e(dmeanz1)) _n
    file write `f' "dmeanz0," %21.0g (e(dmeanz0)) _n
    file write `f' "b_late," %21.0g (el(e(b),1,1)) _n
    file write `f' "b_num," %21.0g (el(e(b),1,2)) _n
    file write `f' "b_denom," %21.0g (el(e(b),1,3)) _n
    file write `f' "v_late," %21.0g (el(e(V),1,1)) _n
    file write `f' "v_num," %21.0g (el(e(V),2,2)) _n
    file write `f' "v_denom," %21.0g (el(e(V),3,3)) _n
    if e(N_clust) != . {
        file write `f' "N_clust," %21.0g (e(N_clust)) _n
    }
    file close `f'
    di as txt "wrote fixture: `id'"
end

* --------------------------------------------------------------------------
* Helper: run a scenario, tolerating failures so the grid completes
* --------------------------------------------------------------------------
capture program drop runfix
program define runfix
    gettoken id 0 : 0
    cap noi `0'
    if _rc {
        di as err "SCENARIO FAILED (`id'): rc = " _rc
    }
    else {
        fixout `id'
    }
end

* --------------------------------------------------------------------------
* Data preparation (identical processing is replicated on the R side)
* --------------------------------------------------------------------------
use "https://people.brandeis.edu/~tslocz/sipp.dta", clear
drop if kwage == . | educ == . | rsncode == 999
generate double lwage = ln(kwage)

* Binary outcome for logit models: above-median log wage
quietly summarize lwage, detail
generate byte hiwage = lwage > r(p50)

* Deterministic, value-based sampling weight (independent of sort order)
generate double wpw = 1 + (kwage - floor(kwage))

* Cluster variable for vce(cluster) scenarios
* (educ has a modest number of distinct values)
generate cluvar = educ

* ==========================================================================
* LATE scenarios
* ==========================================================================

* --- IPWRA (default) ---
runfix late_ipwra_lin_logit_logit drlate (lwage age_5) (nvstat age_5) (rsncode age_5)
runfix late_ipwra_lin_lin_logit   drlate (lwage age_5) (nvstat age_5, linear) (rsncode age_5)
runfix late_ipwra_lin_pois_logit  drlate (lwage age_5) (nvstat age_5, poisson) (rsncode age_5)
runfix late_ipwra_logit_logit_logit drlate (hiwage age_5, logit) (nvstat age_5) (rsncode age_5)
runfix late_ipwra_pois_logit_logit  drlate (kwage age_5, poisson) (nvstat age_5) (rsncode age_5)

* --- IPW (no outcome/treatment covariates) ---
runfix late_ipw_nrm   drlate (lwage) (nvstat) (rsncode age_5), method(ipw)
runfix late_ipw_unnrm drlate (lwage) (nvstat) (rsncode age_5), method(ipw) unnrm

* --- AIPW ---
runfix late_aipw_nrm   drlate (lwage age_5) (nvstat age_5) (rsncode age_5), method(aipw)
runfix late_aipw_unnrm drlate (lwage age_5) (nvstat age_5) (rsncode age_5), method(aipw) unnrm
runfix late_aipw_nrm_logit_y drlate (hiwage age_5, logit) (nvstat age_5) (rsncode age_5), method(aipw)
runfix late_aipw_nrm_pois_y  drlate (kwage age_5, poisson) (nvstat age_5) (rsncode age_5), method(aipw)

* --- RA (no instrument covariates) ---
runfix late_ra          drlate (lwage age_5) (nvstat age_5) (rsncode), method(ra)
runfix late_ra_pois_y   drlate (kwage age_5, poisson) (nvstat age_5) (rsncode), method(ra)

* --- CBPS instrument model ---
runfix late_ipwra_cbps drlate (lwage age_5) (nvstat age_5) (rsncode age_5, cbps)
runfix late_ipw_cbps   drlate (lwage) (nvstat) (rsncode age_5, cbps), method(ipw)
runfix late_aipw_cbps  drlate (lwage age_5) (nvstat age_5) (rsncode age_5, cbps), method(aipw)

* --- IPT instrument model ---
runfix late_ipwra_ipt drlate (lwage age_5) (nvstat age_5) (rsncode age_5, ipt)
runfix late_ipw_ipt   drlate (lwage) (nvstat) (rsncode age_5, ipt), method(ipw)
runfix late_aipw_ipt  drlate (lwage age_5) (nvstat age_5) (rsncode age_5, ipt), method(aipw)

* --- Weights and clustering ---
runfix late_ipwra_pw      drlate (lwage age_5) (nvstat age_5) (rsncode age_5) [pw=wpw]
runfix late_ipwra_cluster drlate (lwage age_5) (nvstat age_5) (rsncode age_5), vce(cluster cluvar)
runfix late_ipwra_pw_cluster drlate (lwage age_5) (nvstat age_5) (rsncode age_5) [pw=wpw], vce(cluster cluvar)
runfix late_aipw_pw       drlate (lwage age_5) (nvstat age_5) (rsncode age_5) [pw=wpw], method(aipw)

* --- Multiple covariates ---
runfix late_ipwra_multix drlate (lwage age_5 educ) (nvstat age_5 educ) (rsncode age_5 educ)

* ==========================================================================
* LATT scenarios
* ==========================================================================

runfix latt_ipwra      drlate (lwage age_5) (nvstat age_5) (rsncode age_5), latt
runfix latt_ipw_nrm    drlate (lwage) (nvstat) (rsncode age_5), latt method(ipw)
runfix latt_ipw_unnrm  drlate (lwage) (nvstat) (rsncode age_5), latt method(ipw) unnrm
runfix latt_aipw_nrm   drlate (lwage age_5) (nvstat age_5) (rsncode age_5), latt method(aipw)
runfix latt_aipw_unnrm drlate (lwage age_5) (nvstat age_5) (rsncode age_5), latt method(aipw) unnrm
runfix latt_ra         drlate (lwage age_5) (nvstat age_5) (rsncode), latt method(ra)
runfix latt_ipwra_ipt  drlate (lwage age_5) (nvstat age_5) (rsncode age_5, ipt), latt
runfix latt_ipwra_pw   drlate (lwage age_5) (nvstat age_5) (rsncode age_5) [pw=wpw], latt
runfix latt_ipwra_cluster drlate (lwage age_5) (nvstat age_5) (rsncode age_5), latt vce(cluster cluvar)

di as txt _n "Done. Fixtures written to tests/testthat/fixtures/."
