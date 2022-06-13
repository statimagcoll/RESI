library(car)
library(splines)
library(sandwich)
library(haven)
library(lme4)
library(dplyr)
library(survival)

data = read.csv("./tests/testdata/insurance.csv")
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
data[,"smoker"] = dplyr::recode(data[,"smoker"], yes = 1, no = 0)
mod.log = glm(smoker ~ age + region, data = data, family = "binomial")

#nls
x <- seq(0, 100, 1)
y<-((runif(1, 10, 20)*x)/(runif(1, 0, 10) + x)) +
  rnorm(101, 0, 1)
data.nls = data.frame(X = x, Y = y)
mod.nls <-nls(Y~a * X/(b + X), data = data.nls, start = list(a = 1, b = 2))

#survival
data.surv = survival::lung
mod.surv = survreg(Surv(time, status) ~ age + sex + ph.karno, data=data.surv,
                   dist='weibull', robust = TRUE)
mod.surv.red = survreg(Surv(time, status) ~ age, data=data.surv,
                   dist='weibull', robust = TRUE)
mod.surv.nr = survreg(Surv(time, status) ~ age + sex + ph.karno, data=data.surv,
                      dist='weibull', robust = FALSE)
mod.coxph =  coxph(Surv(time, status) ~ age + sex + wt.loss, data=lung, robust = TRUE)
mod.coxph.red =  coxph(Surv(time, status) ~ age, data=lung, robust = TRUE)
mod.coxph.nr = coxph(Surv(time, status) ~ age + sex + wt.loss, data=lung, robust = FALSE)

test_that("anoes default produces correct anova table", {
  anoes = anoes(mod, nboot = 10)
  anoes.s = anoes(mod.s, nboot = 10)
  anoes.anova = data.frame(anoes$anova[1:(nrow(anoes$anova))-1,1:3])
  anoes.s.anova = data.frame(anoes.s$anova[1:(nrow(anoes$anova))-1,1:3])
  Anova = data.frame(Anova(mod, test.statistic = 'Wald', vcov. = sandwich::vcovHC))
  Anova.s = data.frame(Anova(mod.s, test.statistic = 'Wald', vcov. = sandwich::vcovHC))
  expect_equal(anoes.anova, Anova)
  expect_equal(anoes.s.anova, Anova.s)
})

test_that("can specify additional arguments of anova", {
  anoes.type3 = anoes(mod, type = 3, nboot = 10)
  anoes.anova.type3 = data.frame(anoes.type3$anova[1:(nrow(anoes.type3$anova))-1,1:3])
  anova.type3 = data.frame(Anova(mod, test.statistic = 'Wald', vcov. = sandwich::vcovHC, type = 3))
  expect_equal(anoes.anova.type3, anova.type3)
})

test_that("the RESI estimates are in between the confidence intervals", {
  anoes = anoes(mod, nboot = 10)
  summary.cols = anoes$summary[, c('RESI', 'LL', 'UL')]
  anova.cols = anoes$anova[,c('RESI', 'LL', 'UL')]
  expect_true(all(summary.cols[,1] >= summary.cols[,2]) & all(summary.cols[,1] <= summary.cols[,3]))
  expect_true(all(anova.cols[,1] >= anova.cols[,2]) & all(anova.cols[,1] <= anova.cols[,3]))
})

test_that("wald.test and waldtest still consistent",{
  expect_equal(unname(lmtest::waldtest(mod.lm, vcov = vcovHC, test = "Chisq")$Chisq[2]),
               unname(aod::wald.test(vcovHC(mod.lm), coef(mod.lm),
                                     Terms = 2:nrow(vcovHC(mod.lm)))$result$chi2["chi2"]))
})

# can use other robust.var besides vcovHC

# binomial family works

# other family works

# an lm works? if changed to the update function instead of glm

# expect warning for specifying F/F for anova/summary

# summary and anova output agree on overall

# summary and anova output agree on variables they should

# summary tables match with/without using a reduced model

# all bootstrap type works
