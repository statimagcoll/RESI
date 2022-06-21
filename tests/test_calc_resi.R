# Test on `calc_resi.R` and `resi.R`
# devtools::install_github("statimagcoll/RESI", ref = "longitudinal", force = TRUE)
# LONGITUDINAL
library(RESI)
library(ggplot2)
library(mvtnorm)
library(gee)
library(geepack)
library(nlme)
library(lme4)


set.seed(1213)
sim = sim_data_cont(N = 200, S = c(1, 1, 1), pi = 0.5, ni_range = c(3, 10), rho.G = 0, sigma0 = 1, sigma.t = 0, sigma.e = 1, rho.e = 0, fixed.design = TRUE)
test_data = sim$data

# testing calc_resi.lme, gee, ...
fit_lme = lme(fixed = y ~ time + trt, data = test_data,
              random = reStruct(~ 1 | id, pdClass= "pdSymm", REML = FALSE))
summary(fit_lme)
# fit_gee = gee(y ~ time + trt, data = test_data, id = id, corstr = "exchangeable")
fit_gee = geeglm(y ~ time + trt, data = test_data, id = id, corstr = "exchangeable")

fit_lmer = lmer(y ~  time  + trt + (1|id) , data = test_data)

calc_resi(fit_gee)
calc_resi(fit_lme)
calc_resi(fit_lmer)
calc_resi(fit_lmer, robust.var = FALSE)

fit_lmer = lmer(y ~ -1 + time  + trt + (1|id) , data = test_data)
calc_resi(fit_lmer, robust.var = TRUE)
calc_resi(fit_lmer, robust.var = FALSE)


# testing resi.geeglm, resi.lmerMod
resi(fit_gee, nboot = 1000)
resi(fit_lmer, nboot = 10)
