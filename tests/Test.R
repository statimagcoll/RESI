# test
# Testing `sim_data_cont`
library(RESI)
library(ggplot2)
library(mvtnorm)
library(gee)
library(geepack)
library(nlme)
library(lme4)
library(dplyr)
# devtools::install_github("statimagcoll/RESI", ref = "longitudinal", force = TRUE)


set.seed(1213)

sim = sim_data_cont(N = 200, S = c(1, 1, 1), pi = 0.5, ni_range = c(3, 10), rho.G = 0, sigma0 = 1, sigma.t = 0, sigma.e = 1, rho.e = 0, fixed.design = TRUE)
test_data = sim$data

# testing calc_resi.lme, gee, ...
fit_lme = lme(fixed = y ~ time + trt, data = test_data,
              random = reStruct(~ 1 | id, pdClass= "pdSymm", REML = FALSE))
summary(fit_lme)
fit_gee = gee(y ~ time + trt, data = test_data, id = id, corstr = "exchangeable")

fit_lmer = lmer(y ~ -1 + time  + trt + (1|id) , data = test_data)


calc_resi(fit_gee)
calc_resi(fit_lme)
calc_resi(fit_lmer, robust.var = FALSE)



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
S = rep(0, 3)

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



# testing `ess`
N = 20
pi = 0.5
ni_range = c(5, 5)
rho.G = 0.2
sigma0 = 1
sigma.t = 0
true_sigma.e = 1
rho.e = 0
fixed.design = TRUE
S = rep(1, 3)

sim = sim_data_cont(N = N, S = S, pi = pi, ni_range = ni_range, rho.G = 0, sigma0 = 1, sigma.t = 0, sigma.e = 1, rho.e = 0, fixed.design = TRUE)

test_data = sim$data

fit = geepack::geeglm(y ~ time + trt, data = test_data, id = id, corstr = "exchangeable")

obj = fit
object = fit
xTab = output

# -------------------
library(magrittr)
library(Matrix)
library(splines)

# Testing data simulation function
N = 200
pi = 0.5
ni_range = c(5, 5)
ni = 5
rho.G = 0.2
sigma0 = 1
sigma.t = 0
sigma.e = 1
rho.e = 0
fixed.trt = fixed.time = FALSE
S = rep(1, 3)

corstr.e = "exchangeable" # the true structure
cor_spec = "exchangeable" # the working cor structrue


AR1 = AR1_str(0.5, 5)
comp = comp_str(1, 0.5, 5)

true_sigma.e = comp_str(1, 0.5, 5)
work_sigma.e = comp_to_AR1(true_sigma.e)

time_points = (0:(ni_range[2]-1))
# The true cov of \hat{beta} given X and working covariance structure
# Using numeric methods to find the variance of each parameter estimator
# (assuming the model is correctly specified)
# The covariance matrix of \hat{beta}
sum_A = sum_B = sum_ind = 0
temp_A = temp_B = 0 # for the var(\hat{\beta}) for `trt` when time is random
rep = 1e2
temp_N = 100
for (k in 1:rep){
  if (fixed.trt) {
    trt_temp = sample(rep(0:1, times = c(ceiling(temp_N*pi), temp_N - ceiling(temp_N*pi)) ), replace = FALSE)
  } else {trt_temp = rbinom(temp_N, 1, pi)}
  # duplicate it for longitudinal data
  trt_temp = rep(trt_temp, each = ni_range[1])

  if (fixed.time){
    time_temp = rep(0:(ni-1), times = temp_N)
  } else {
    time_list = sapply(1:temp_N, function(x) runif(ni, 0, 5) %>% sort)
    time_temp = c(unlist(time_list))
  }

  X = cbind(1, time = time_temp, trt = trt_temp)

  # for the Cov(beta_hat) under the correlation specification
  big_work_cov = lapply(1:temp_N, function(X) work_sigma.e) %>% bdiag %>% as.matrix
  big_true_cov = lapply(1:temp_N, function(X) true_sigma.e) %>% bdiag %>% as.matrix
  sum_A = sum_A + t(X) %*% solve(big_work_cov) %*% X
  sum_B = sum_B + t(X) %*% solve(big_work_cov) %*% big_true_cov %*% solve(big_work_cov) %*% X

  # for the Cov(beta_hat) under independence specification
  sum_ind = sum_ind + t(X) %*% X

  # Note: the sum's here = sum of (X^T Sigma X) or its corresponding sandwich form
}
COV_beta = solve(sum_A) %*% sum_B %*% solve(sum_A) * rep * temp_N  # = N * (the asymptotic cov of beta hat - beta_0)
COV_beta_ind = solve(sum_ind / (rep*temp_N * ni_range[1])) # =  N * (the asymptotic cov of beta hat under independence )



sim = sim_data_cont(N = N, S = S, pi = pi, ni_range = ni_range, true_sigma.e = comp, work_sigma.e = comp,
                    COV_beta = COV_beta, COV_beta_ind = COV_beta_ind, fixed.trt = TRUE, fixed.time = TRUE)

sim_dat = sim$data


fit_gee = geepack::geeglm(y ~ ns(time, df = 2) + trt, data = sim_dat, id = id, corstr = "exchangeable")

resi_pe.geeglm(fit_gee, anova = TRUE)

cs_dat = sim_dat %>% group_by(id) %>% slice_sample(n = 1)

fit_lm = lm(y ~ ns(time, df = 2) + trt, data = cs_dat)
a = resi(fit_lm, data = cs_dat)


fit_glm = glm(y ~ time + trt, data = sim_dat)
resi_pe.glm(fit_glm, data = sim_dat)





