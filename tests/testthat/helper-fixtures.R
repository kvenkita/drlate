# Helpers shared by the golden-fixture tests.

#' Path to a Stata-generated fixture; NULL if absent.
fixture_path <- function(id) {
  p <- test_path("fixtures", paste0(id, ".csv"))
  if (file.exists(p)) p else NULL
}

#' Read a fixture CSV (name,value rows) into a named numeric vector.
read_fixture <- function(id) {
  p <- fixture_path(id)
  x <- utils::read.csv(p, stringsAsFactors = FALSE)
  stats::setNames(as.numeric(x$value), x$name)
}

#' Download (once per session) and process the SIPP extract exactly as
#' inst/stata/make-fixtures.do does. Returns NULL when offline/unavailable.
sipp_data <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    rds <- test_path("fixtures", "sipp.rds")
    if (file.exists(rds)) {
      cache <<- readRDS(rds)
      return(cache)
    }
    if (!requireNamespace("haven", quietly = TRUE)) return(NULL)
    # Prefer the local copy used by make-fixtures.do; fall back to the URL
    local_dta <- test_path("fixtures", "sipp.dta")
    src <- if (file.exists(local_dta)) local_dta
           else "https://people.brandeis.edu/~tslocz/sipp.dta"
    dta <- tryCatch(haven::read_dta(src), error = function(e) NULL)
    if (is.null(dta)) return(NULL)
    d <- as.data.frame(dta)
    d <- d[!is.na(d$kwage) & !is.na(d$educ) & d$rsncode != 999, ]
    d$lwage <- log(d$kwage)
    # Stata: summarize, detail; gen hiwage = lwage > r(p50)
    d$hiwage <- as.integer(d$lwage > stats::median(d$lwage))
    # Deterministic value-based weight, as in the do-file
    d$wpw <- 1 + (d$kwage - floor(d$kwage))
    d$cluvar <- d$educ
    saveRDS(d, rds)
    cache <<- d
    d
  }
})

#' Compare an R drlate fit against a Stata fixture.
expect_matches_fixture <- function(fit, id,
                                   tol_b = 1e-6, tol_se = 1e-4) {
  fx <- read_fixture(id)
  b <- unname(coef(fit))
  se <- unname(sqrt(diag(fit$vcov3)))
  expect_equal(fit$N, unname(fx["N"]))
  expect_equal(unname(fit$dmeanz1), unname(fx["dmeanz1"]), tolerance = 1e-6)
  expect_equal(unname(fit$dmeanz0), unname(fx["dmeanz0"]), tolerance = 1e-6)
  expect_equal(b[1], unname(fx["b_late"]),  tolerance = tol_b)
  expect_equal(b[2], unname(fx["b_num"]),   tolerance = tol_b)
  expect_equal(b[3], unname(fx["b_denom"]), tolerance = tol_b)
  expect_equal(se[1], unname(sqrt(fx["v_late"])),  tolerance = tol_se)
  expect_equal(se[2], unname(sqrt(fx["v_num"])),   tolerance = tol_se)
  expect_equal(se[3], unname(sqrt(fx["v_denom"])), tolerance = tol_se)
  if ("N_clust" %in% names(fx)) {
    expect_equal(fit$N_clust, unname(fx["N_clust"]))
  }
}

#' Skip unless both the fixture and the SIPP data are available.
skip_if_no_fixture <- function(id) {
  skip_on_cran()
  if (is.null(fixture_path(id))) {
    skip(paste0("fixture '", id, "' not generated (run ",
                "inst/stata/make-fixtures.do in Stata)"))
  }
  if (is.null(sipp_data())) skip("SIPP data unavailable (offline?)")
}
