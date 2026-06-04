# Verify the numerical claims made in the primer vignette
setwd("C:/Users/kyle/Documents/Projects/Personal/drlate")
devtools::load_all(quiet = TRUE)
data(drlate_sim)

cat("naive OLS:",
    coef(lm(lwage ~ nvstat + age + educ, drlate_sim))["nvstat"], "\n")
cat("raw Wald:",
    with(drlate_sim,
         (mean(lwage[rsncode == 1]) - mean(lwage[rsncode == 0])) /
         (mean(nvstat[rsncode == 1]) - mean(nvstat[rsncode == 0]))), "\n")

f <- drlate(lwage ~ age + educ, nvstat ~ age + educ, rsncode ~ age + educ,
            data = drlate_sim)
cat("ipwra:", coef(f)[1], " CI:", confint(f)[1, ], "\n")

cat("PS wrong:",
    coef(drlate(lwage ~ age + educ, nvstat ~ age + educ, rsncode ~ 1,
                data = drlate_sim))[1], "\n")
cat("regs wrong:",
    coef(drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age + educ,
                data = drlate_sim))[1], "\n")
cat("ipt:",
    coef(drlate(lwage ~ age + educ, nvstat ~ age + educ,
                rsncode ~ age + educ, data = drlate_sim,
                ivmodel = "ipt"))[1], "\n")
cat("latt:",
    coef(drlate(lwage ~ age + educ, nvstat ~ age + educ,
                rsncode ~ age + educ, data = drlate_sim,
                estimand = "latt"))[1], "\n")
cat("poisson kwage:",
    coef(drlate(kwage ~ age + educ, nvstat ~ age + educ,
                rsncode ~ age + educ, data = drlate_sim,
                omodel = "poisson"))[1], "\n")
