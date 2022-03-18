# test
# Testing `sim_data_cont`
library(RESI)
library(ggplot2)
library(mvtnorm)
library(gee)
library(geepack)
library(nlme)

set.seed(1213)
sim = sim_data_cont(N = 200, beta = c(-1, 1, 2), pi = 0.5, ni_range = c(3, 10), rho.G = 0, sigma0 = 1, sigma.t = 0, sigma.e = 1, rho.e = 0, fixed.design = TRUE)
test_data = sim$data
# ggplot(test_data, aes(x = time, y = y, col = trt, group = id)) + geom_line()
# a random intercept only model
test_data$subj_id = test_data$id

lme()

mod.gee <- geeglm(y ~ time + trt, data = test_data, id = subj_id, corstr = "exchangeable")
A = summary(mod.gee)

mod1 <- gee(y ~ time + trt, data = test_data, id = subj_id, corstr = "exchangeable")
summary(mod1)
mod.gee2 <- geeglm(y ~ time + trt, data = test_data, id = subj_id, corstr = "independence")
summary(mod.gee2)

mod.lme <- lme(fixed = y ~ time + trt, method = "ML", data = test_data,
            random = reStruct(~ time | subj_id, pdClass= "pdSymm", REML = FALSE))


library(RESI)
anoes.lme(mod.lme)
