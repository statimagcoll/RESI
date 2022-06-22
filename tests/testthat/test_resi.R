library(car)
library(splines)
library(sandwich)
library(haven)
library(lme4)
library(dplyr)
library(survival)

data = read.csv(test_path("testdata", "insurance.csv"))
mod = glm(charges ~ region * age + bmi + sex, data = data)
mod.s = glm(charges ~ region * ns(age, df=3) + bmi + sex, data = data)
mod.r = glm(charges ~ bmi +sex, data = data)
mod.s.r = glm(charges ~ region * ns(age, df=3), data = data)
mod.lm = lm(charges ~ region * age + bmi + sex, data = data)
mod.lm.s = lm(charges ~ region*ns(age, df = 3) + bmi + sex, data = data)
mod.lm.s.r = lm(charges ~ region * ns(age, df=3), data = data)
mod.one.predictor = glm(charges ~ bmi, data = data)
mod.one.pred.multi = glm(charges ~ region, data = data)
mod.one.pred.lm = lm(charges ~ region, data = data)
data[,"smoker"] = recode(data[,"smoker"], yes = 1, no = 0)
mod.log = glm(smoker ~ age + region, data = data, family = "binomial")

# missing data models
data.nas = data
data.nas[12:17, 'bmi'] = NA
mod.na = lm(charges ~ region * age + bmi + sex, data = data.nas)
mod.log.na = glm(smoker ~ age + region + bmi, data = data.nas, family = "binomial")

#nls
x <- seq(0, 100, 1)
y<-((runif(1, 10, 20)*x)/(runif(1, 0, 10) + x)) +
  rnorm(101, 0, 1)
data.nls = data.frame(X = x, Y = y)
mod.nls <-nls(Y~a * X/(b + X), data = data.nls, start = list(a = 1, b = 2))

#survival
data.surv = survival::lung
data.surv.test = data.surv[203:nrow(data.surv),]
rownames(data.surv.test) = letters[1:26]
mod.surv = survreg(Surv(time, status) ~ age + sex + ph.karno, data=data.surv,
                   dist='weibull', robust = TRUE)
# check for missing value
mod.surv.red = survreg(Surv(time, status) ~ age, data=data.surv,
                   dist='weibull', robust = TRUE)
mod.surv.red.test = survreg(Surv(time, status) ~ age + sex + ph.karno, data=data.surv[which(1:nrow(data.surv)!=mod.surv$na.action),],dist='weibull', robust = TRUE)

mod.surv.nr = survreg(Surv(time, status) ~ age + sex + ph.karno, data=data.surv,
                      dist='weibull', robust = FALSE)
mod.surv.nr.red.test = survreg(Surv(time, status) ~ age + sex + ph.karno, data=data.surv[which(rownames(data.surv)!=mod.surv$na.action),],dist='weibull', robust = FALSE)
mod.coxph =  coxph(Surv(time, status) ~ age + sex + wt.loss, data=lung, robust = TRUE)
mod.coxph.red =  coxph(Surv(time, status) ~ age, data=lung, robust = TRUE)
mod.coxph.nr = coxph(Surv(time, status) ~ age + sex + wt.loss, data=lung, robust = FALSE)
mod.coxph.na = coxph(Surv(time, status) ~ age + ph.karno + wt.loss, data=lung, robust = TRUE)

# hurdle
data.hurdle = pscl::bioChemists
data.hurdle.na = data.hurdle
data.hurdle.na[12:17, 'kid5'] = NA
mod.hurdle = pscl::hurdle(art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment, data = data.hurdle)
mod.hurd.r = pscl::hurdle(art ~ kid5 + phd | kid5 + phd, data = data.hurdle)
mod.hurdle.na = pscl::hurdle(art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment, data = data.hurdle.na)
mod.hurdle.r.na = pscl::hurdle(art ~ mar + phd | mar + phd, data = data.hurdle.na)

# zeroinfl
mod.zinf = pscl::zeroinfl(art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment, data = data.hurdle)
mod.zinf.r = pscl::zeroinfl(art ~ kid5 + phd | kid5 + phd, data = data.hurdle)

# gee
data.gee = read.csv(test_path("testdata", "depression.csv"))
mod.gee = gee::gee(depression ~ diagnose + drug*time, data = data.gee, id = id, family = binomial, corstr = "independence")
mod.geeglm = geepack::geeglm(depression ~ diagnose + drug*time, data = data.gee, id = id, family = binomial, corstr = "independence")
mod.lme = nlme::lme(distance ~ age + Sex, data = nlme::Orthodont, random = ~ 1)


## tests
test_that("Specifying non-allowed vcov produces warning",{
  expect_warning(resi(mod.nls, data = data.nls, nboot = 10, vcov = sandwich::vcovHC), "Sandwich vcov function not applicable for nls model type, vcovfunc set to regtools::nlshc")
  expect_warning(resi(mod.surv, data = data.surv, nboot = 10, vcov = sandwich::vcovHC), "vcovfunc argument ignored for survreg objects")
  expect_warning(resi(mod.coxph, data = data.surv, nboot = 10, vcov = sandwich::vcovHC), "vcovfunc argument ignored for coxph objects")
})
