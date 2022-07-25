library(car)
library(splines)
library(sandwich)
library(survival)
library(nlme)
library(lme4)

data = RESI::insurance
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
mod.log = glm(smoker ~ age + region, data = data, family = "binomial")
mod.pois = glm(age ~ region + charges + sex*smoker, data = data, family = "poisson")

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
mod.surv.nr = survreg(Surv(time, status) ~ age + sex + ph.karno, data=data.surv,
                      dist='weibull', robust = FALSE)
mod.survexp = survreg(Surv(time, status) ~ age + sex + ph.karno, data=data.surv,
                      dist='exp', robust = FALSE)
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
data.gee = RESI::depression
mod.gee = gee::gee(depression ~ diagnose + drug*time, data = data.gee, id = id, family = binomial, corstr = "independence")
mod.geeglm = geepack::geeglm(depression ~ diagnose + drug*time, data = data.gee, id = id, family = binomial, corstr = "independence")
mod.lme = nlme::lme(distance ~ age + Sex, data = nlme::Orthodont, random = ~ 1)
mod.lmerMod = lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
data.lmer = nlme::Orthodont
data.lmer$nsex <- as.numeric(data.lmer$Sex=="Male")
data.lmer$nsexage <- with(data.lmer, nsex*age)
mod.lmerMod2 = lmer(distance ~ age + (age|Subject) + (0+nsex|Subject) +
                      (0 + nsexage|Subject), data=data.lmer)

## tests
test_that("Specifying non-allowed vcov produces warning",{
  expect_warning(resi(mod.nls, data = data.nls, nboot = 10, vcovfunc = sandwich::vcovHC), "Sandwich vcov function not applicable for nls model type, vcovfunc set to regtools::nlshc")
  expect_warning(resi(mod.surv, data = data.surv, nboot = 10, vcovfunc = sandwich::vcovHC), "vcovfunc argument ignored for survreg objects")
  expect_warning(resi(mod.coxph, data = data.surv, nboot = 10, vcovfunc = sandwich::vcovHC), "vcovfunc argument ignored for coxph objects")
})

test_that("Bayesian bootstrap is not allowed for survival models",{
  expect_warning(resi(mod.surv, data = data.surv, nboot = 10, boot.method = "bayes"), "Bayesian bootstrap not currently supported for survreg models, using non-parametric bootstrap")
  expect_warning(resi(mod.coxph, data = data.surv, nboot = 10, boot.method = "bayes"), "Bayesian bootstrap not currently supported for coxph models, using non-parametric bootstrap")
})

test_that("resi produces the correct estimates", {
  expect_equal(unname(resi(mod, nboot = 1)$estimates), c(0.35982590, -0.06733537, -0.02670248, -0.03341748, -0.00246893, 0.14980504, 0.15238719,
                                                         0.05839381, 0.01667149, 0.03610640, -0.01497171, 0.03655364, 0.29533571, 0.14991488,
                                                         0.05159896, 0.01617303), tolerance = 1e-07)
  expect_equal(unname(resi(mod.lm, nboot = 1)$estimates), c(0.359825898, -0.067297337, -0.026687397, -0.033398602, -0.002467535, 0.149720414, 0.152301107,
                                                            0.058360820, 0.016662071, 0.036086004, -0.014963251, 0.036616942, 0.296220354, 0.150361134,
                                                            0.051742893, 0.016116372))
  expect_equal(unname(resi(mod.s, nboot = 1, data = data)$estimates[c(1, 6, 11, 24)]), c(0.35634034, 0.06143486, 0.00748577, 0.05339006), tolerance = 1e-07)
  # expect_equal(unname(resi(mod.nls, nboot = 1, data = data.nls)$estimates), c(11.763649, 8.979576, 1.352517), tolerance = 1e-07)
  # need to figure out what is happening with nls
  expect_equal(unname(resi(mod.surv, nboot = 1, data = data.surv)$estimates), c(0.19357703, 0.51465342, -0.08062289, 0.20010885, 0.10512266, -0.27578597, 0.04597918,
                                                                                0.18883217, 0.08163828), tolerance = 1e-07)
  expect_equal(unname(resi(mod.coxph, nboot = 1, data = data.surv)$estimates), c(0.224354226, 0.136138740, -0.213293162, 0.008641209, 0.117732151, 0.202042262, 0.000000000))
  expect_equal(unname(resi(mod.hurdle, nboot = 1)$estimates[c(1, 2, 5, 8, 10)]), c(0.28063439, 0.12770489, -0.06898884, 0.02686305, 0.06039791), tolerance = 1e-07)
  expect_equal(unname(resi(mod.zinf, nboot = 1)$estimates[c(1, 2, 5, 8, 10)]), c(0.23725614, 0.11892936, -0.07014934, -0.03481673, -0.03461235), tolerance = 1e-07)
  expect_equal(unname(resi(mod.gee, nboot = 10, data = data.gee)$coefficients[,'RESI']), c(0.0000, 0.4850899, 0.0000, 0.5591547, 0.2889370), tolerance = 1e-07)
  expect_equal(unname(resi(mod.geeglm, nboot = 10, data = data.gee)$coefficients[,'RESI']), c(0.0000, 0.4850899, 0.000, 0.5591547, 0.2889370), tolerance = 1e-07)
  expect_equal(unname(resi(mod.lme, nboot = 10, data = data.gee)$coefficients[,'RESI']), c(3.659090, 1.739166, 0.512371), tolerance = 1e-07)
})

test_that("RESI estimates are in between the confidence limits", {
  resi.obj = resi(mod, nboot = 500)
  expect_true(all(resi.obj$coefficients$RESI >= resi.obj$coefficients$`2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`97.5%`))
  expect_true(all(resi.obj$anova$RESI >= resi.obj$anova$`2.5%`) & all(resi.obj$anova$RESI <= resi.obj$anova$`97.5%`))
  # bayes option for glm not seeming right currently
  # resi.obj = resi(mod, boot.method = "bayes")
  # expect_true(all(resi.obj$coefficients$RESI >= resi.obj$coefficients$`2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`97.5%`))
  # expect_true(all(resi.obj$anova$RESI >= resi.obj$anova$`2.5%`) & all(resi.obj$anova$RESI <= resi.obj$anova$`97.5%`))
  resi.obj = resi(mod.lm, nboot = 500)
  expect_true(all(resi.obj$coefficients$RESI >= resi.obj$coefficients$`2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`97.5%`))
  expect_true(all(resi.obj$anova[1:5, "RESI"] >= resi.obj$anova[1:5, "2.5%"]) & all(resi.obj$anova[1:5, "RESI"] <= resi.obj$anova[1:5, "RESI"]))
  resi.obj = resi(mod.nls, data = data.nls)
  expect_true(all(resi.obj$coefficients$RESI >= resi.obj$coefficients$`2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`97.5%`))
  resi.obj = resi(mod.surv, data = data.surv, nboot = 500)
  expect_true(all(resi.obj$coefficients$RESI >= resi.obj$coefficients$`2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`97.5%`))
  expect_true(all(resi.obj$anova$RESI >= resi.obj$anova$`2.5%`) & all(resi.obj$anova$RESI <= resi.obj$anova$`97.5%`))
  resi.obj = resi(mod.coxph, data = data.surv, nboot = 500)
  expect_true(all(resi.obj$coefficients$RESI >= resi.obj$coefficients$`2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`97.5%`))
  expect_true(all(resi.obj$anova$RESI >= resi.obj$anova$`2.5%`) & all(resi.obj$anova$RESI <= resi.obj$anova$`97.5%`))
  resi.obj = resi(mod.hurdle, nboot = 500)
  expect_true(all(resi.obj$coefficients$RESI >= resi.obj$coefficients$`2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`97.5%`))
  # resi.obj = resi(mod.hurdle, boot.method = "bayes", nboot = 500)
  # expect_true(all(resi.obj$coefficients$RESI >= resi.obj$coefficients$`2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`97.5%`))
  resi.obj = resi(mod.zinf, nboot = 500)
  expect_true(all(resi.obj$coefficients$RESI >= resi.obj$coefficients$`2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`97.5%`))
  # resi.obj = resi(mod.gee, data = data.gee, nboot = 500)
  # expect_true(all(resi.obj$coefficients$RESI >= resi.obj$coefficients$`2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`97.5%`))
  resi.obj = resi(mod.geeglm, data = data.gee, nboot = 500)
  expect_true(all(resi.obj$coefficients$RESI >= resi.obj$coefficients$`2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`97.5%`))
  resi.obj = resi(mod.lme, nboot = 500)
  expect_true(all(resi.obj$coefficients$RESI >= resi.obj$coefficients$`2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`97.5%`))
})

test_that("z2S.form = 2 returns same abs. RESI as Chi-sq",{
  resi.obj = resi(mod, nboot = 1, z2S.form = 2)
  expect_equal(abs(resi.obj$coefficients[6:7, 'RESI']), resi.obj$anova[3:4, 'RESI'], tolerance = 1e-07)
  resi.obj = resi(mod.surv, data = data.surv, nboot = 1, z2S.form = 2)
  expect_equal(abs(resi.obj$coefficients[2:4, 'RESI']), resi.obj$anova[1:3, 'RESI'], tolerance = 1e-07)
  resi.obj = resi(mod.coxph, data = data.surv, nboot = 1, z2S.form = 2)
  expect_equal(abs(resi.obj$coefficients[1:3, 'RESI']), resi.obj$anova[1:3, 'RESI'], tolerance = 1e-07)
})

test_that("t2S.form = 2 returns same abs. RESI as F", {
  resi.obj = resi(mod.lm, nboot = 1, t2S.form = 2)
  expect_equal(abs(resi.obj$coefficients[6:7,'RESI']), resi.obj$anova[3:4,'RESI'], tolerance = 1e-07)
})

test_that("wald.test and waldtest still consistent",{
  expect_equal(unname(lmtest::waldtest(mod.lm, vcov = vcovHC, test = "Chisq")$Chisq[2]),
               unname(aod::wald.test(vcovHC(mod.lm), coef(mod.lm),
                                     Terms = 2:nrow(vcovHC(mod.lm)))$result$chi2["chi2"]))
})

