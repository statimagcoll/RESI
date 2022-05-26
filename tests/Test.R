# test
# Testing `sim_data_cont`
library(RESI)
library(ggplot2)
library(mvtnorm)
library(gee)
library(geepack)
library(nlme)
# devtools::install_github("statimagcoll/RESI", ref = "longitudinal", force = TRUE)


set.seed(1213)

sim = sim_data_cont(N = 200, beta = c(-1, 1, 2), pi = 0.5, ni_range = c(3, 10), rho.G = 0, sigma0 = 1, sigma.t = 0, sigma.e = 1, rho.e = 0, fixed.design = TRUE)
test_data = sim$data
# ggplot(test_data, aes(x = time, y = y, col = trt, group = id)) + geom_line()
# a random intercept only model
# test_data$subj_id = test_data$id

mod1 <- gee(y ~ time + trt, data = test_data, id = id, corstr = "exchangeable")
summary(mod1)
mod.gee2 <- geeglm(y ~ time + trt, data = test_data, id = subj_id, corstr = "independence")
summary(mod.gee2)



# ----------------------
sim = sim_data_cont(N = 200, beta = c(-1, 1, 2), pi = 0.5, ni_range = c(3, 10), rho.G = 0, sigma0 = 1, sigma.t = 0, sigma.e = 1, rho.e = 0, fixed.design = TRUE)
test_data = sim$data

mod.gee <- geeglm(y ~ time + trt, data = test_data, id = id, corstr = "exchangeable")
A = summary(mod.gee)

mod.lme <- lme(fixed = y ~ time + trt, method = "ML", data = test_data,
            random = reStruct(~ time | id, pdClass= "pdSymm", REML = FALSE))

ess(mod.lme)
ess(mod.gee)

# ------------
# Testing simulated data
N = 20
pi = 0.5
ni_range = c(5, 5)
rho.G = 0.2
sigma0 = 1
sigma.t = 0
sigma.e = 1
rho.e = 0
fixed.design = TRUE
S = rep(0.5, 3)

sim = sim_data_cont(N = N, S = S, pi = pi, ni_range = ni_range, rho.G = 0, sigma0 = 1, sigma.t = 0, sigma.e = 1, rho.e = 0, fixed.design = TRUE)

test_data = sim$data
fit <- lme(fixed = y ~ time + trt, method = "ML", data = test_data,
    random = reStruct(~ 1 | id, pdClass= "pdSymm"))

RESI::resi(fit)

summary(fit)

sim$true_beta
sim$true_sd


ggplot(test_data, aes(x = time, y = y, col = trt, group = id)) + geom_line()

fit <- lme(fixed = y ~ time + trt, data = test_data,
           random = reStruct(~ 1 | id, pdClass= "pdSymm"))

summary(fit)

RESI::resi(fit)

# Testing Bayesian Bootstrap sampling
y <- factor(c(0, 1, 1, 1, 1 ,0, 1 ,1 ,1 ,0 ,0 ,0 ,0 ,1, 1, 0, 1 ,1 ,1, 0 ,1 ,1, 0 ,0 ,0, 1 ,1))
x <- factor(c(rep(1, 8), rep(2, 6), rep(3, 13)))
x2 <- c(29 ,29 ,14 ,20 ,33 ,20 ,26, 17 ,17 ,24, 28 ,23 ,19, 16 ,42 ,13, 18 ,60, 37, 16 ,36 ,30 ,26 ,17, 64, 32 ,24)
df.resi <- data.frame(y,x,x2)

data = bayes.samp(df.resi)

glm(y ~ x + x2, data = data, family = "binomial", weights = g)



