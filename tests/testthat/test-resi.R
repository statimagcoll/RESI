library(car)
#library(splines)
library(sandwich)
#library(survival)
library(nlme)
#library(lme4)
library(tibble)
#library(geepack)

data = RESI::insurance
mod = glm(charges ~ region * age + bmi + sex, data = data)
mod.r = glm(charges ~ bmi +sex, data = data)
mod.lm = lm(charges ~ region * age + bmi + sex, data = data)
mod.log = glm(smoker ~ age + region, data = data, family = "binomial")
mod.glm.int = glm(smoker ~ 1, data = data, family = "binomial")
mod.lm.int = lm(charges ~ 1, data = data)
# splines
if(requireNamespace("splines")){
  mod.s = glm(charges ~ region * splines::ns(age, df=3) + bmi + sex, data = data)
  mod.s.r = glm(charges ~ region * splines::ns(age, df=3), data = data)
  mod.lm.s = lm(charges ~ region*splines::ns(age, df = 3) + bmi + sex, data = data)
  mod.lm.s.r = lm(charges ~ region * splines::ns(age, df=3), data = data)
}

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
mod.nls = nls(Y~a * X/(b + X), data = data.nls, start = list(a = 1, b = 2))

#survival
if(requireNamespace("survival")){
  data.surv = survival::lung
  mod.surv = survival::survreg(survival::Surv(time, status) ~ age + sex + ph.karno, data=data.surv,
                     dist='weibull', robust = TRUE)
  mod.surv.red = survival::survreg(survival::Surv(time, status) ~ age, data=data.surv[which(!(1:nrow(data.surv)%in% mod.surv$na.action)),],
                         dist='weibull', robust = TRUE)
  mod.surv.nr = survival::survreg(survival::Surv(time, status) ~ age + sex + ph.karno, data=data.surv,
                        dist='weibull', robust = FALSE)
  mod.surv.int = survival::survreg(survival::Surv(time, status) ~ 1, data=data.surv,
                                   dist='weibull', robust = TRUE)
  mod.coxph =  survival::coxph(survival::Surv(time, status) ~ age + sex + wt.loss, data=survival::lung, robust = TRUE)
  mod.coxph.red =  survival::coxph(survival::Surv(time, status) ~ age, data=survival::lung, robust = TRUE)
}


# hurdle
if(requireNamespace("pscl")){
  data.hurdle = pscl::bioChemists
  mod.hurdle = pscl::hurdle(art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment, data = data.hurdle)
  mod.hurdle.r = pscl::hurdle(art ~ kid5 + phd | kid5 + phd, data = data.hurdle)
  mod.hurdle.int = pscl::hurdle(art ~ 1, data = data.hurdle)

  # zeroinfl
  mod.zinf = pscl::zeroinfl(art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment, data = data.hurdle)
  mod.zinf.r = pscl::zeroinfl(art ~ kid5 + phd | kid5 + phd, data = data.hurdle)
  mod.zinf.int = pscl::zeroinfl(art ~ 1, data = data.hurdle)
}

# gee and lme
if(requireNamespace("gee")){
  data.gee = RESI::depression
  mod.gee = gee::gee(depression ~ diagnose + drug*time, data = data.gee,
                     id = id, family = binomial, corstr = "independence")
  mod.gee.int = gee::gee(depression ~ 1, data = data.gee,
                         id = id, family = binomial, corstr = "independence")
}

if(requireNamespace("geepack")){
  library(geepack)
  mod.geeglm = geepack::geeglm(depression ~ diagnose + drug*time, data = data.gee,
                               id = id, family = binomial, corstr = "independence")
  mod.geeglm.r = geepack::geeglm(depression ~ diagnose, data = data.gee, id = id,
                                 family = binomial, corstr = "independence")
  mod.geeglm.int = geepack::geeglm(depression ~ 1, data = data.gee, id = id,
                                 family = binomial, corstr = "independence")
}

mod.lme = nlme::lme(distance ~ age + Sex, data = nlme::Orthodont, random = ~ 1)
mod.lme.int = nlme::lme(distance ~ 1, data = nlme::Orthodont, random = ~ 1)

if(requireNamespace("lme4")){
  mod.lmerMod = lme4::lmer(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
  data.lmer = nlme::Orthodont
  data.lmer$nsex <- as.numeric(data.lmer$Sex=="Male")
  data.lmer$nsexage <- with(data.lmer, nsex*age)
  mod.lmerMod2 = lme4::lmer(distance ~ age + (age|Subject) + (0+nsex|Subject) +
                        (0 + nsexage|Subject), data=data.lmer)
  mod.lmerMod.int = lme4::lmer(Reaction ~ 1 + (1 | Subject), lme4::sleepstudy)
}

# tibbles
if(requireNamespace("tibble")){
  tib.dat <- tibble::as_tibble(data)
  mod.glm.tib <- glm(charges ~ region * age + bmi + sex, data = tib.dat)
  mod.lm.tib <- lm(charges ~ region * age + bmi + sex, data = tib.dat)
  tib.nls <- tibble::as_tibble(data.nls)
  mod.nls.tib <- nls(Y~a * X/(b + X), data = tib.nls, start = list(a = 1, b = 2))
  if(requireNamespace("survival")){
    tib.surv <- tibble::as_tibble(data.surv)
    mod.surv.tib <- survival::survreg(survival::Surv(time, status) ~ age + sex + ph.karno, data=tib.surv,
                            dist='weibull', robust = TRUE)
    mod.coxph.tib <-  survival::coxph(survival::Surv(time, status) ~ age + sex + wt.loss,
                            data=tibble::as_tibble(survival::lung), robust = TRUE)
  }
  if(requireNamespace("pscl")){
    tib.hurdle <- tibble::as_tibble(data.hurdle)
    mod.hurdle.tib <- pscl::hurdle(art ~ fem + mar + kid5 + phd + ment | fem + mar +
                                     kid5 + phd + ment, data = tib.hurdle)
    mod.zinf.tib <- pscl::zeroinfl(art ~ fem + mar + kid5 + phd + ment | fem + mar +
                                     kid5 + phd + ment, data = tib.hurdle)
  }
  if(requireNamespace("gee")){
    tib.gee <- tibble::as_tibble(data.gee)
    mod.gee.tib <- gee::gee(depression ~ diagnose + drug*time, data = tib.gee,
                            id = id, family = binomial, corstr = "independence")
  }
  if(requireNamespace("geepack")){
    mod.geeglm.tib <- geepack::geeglm(depression ~ diagnose + drug*time,
                                      data = tib.gee, id = id, family = binomial,
                                      corstr = "independence")
  }
  if(requireNamespace("lme4")){
    lmerModtibdat = tibble::as_tibble(lme4::sleepstudy)
    mod.lmerMod.tib <- lme4::lmer(Reaction ~ Days + (Days | Subject),
                                  lmerModtibdat)
  }}


## tests
test_that("Specifying non-allowed vcov produces warning",{
  expect_warning(resi(mod.nls, data = data.nls, nboot = 1, vcovfunc = sandwich::vcovHC),
                 "Sandwich vcov function not applicable for nls model type, vcovfunc set to regtools::nlshc")
  if(requireNamespace("survival")){
  expect_warning(resi(mod.surv, data = data.surv, nboot = 1, vcovfunc = sandwich::vcovHC),
                 "vcovfunc argument ignored for survreg objects")
  expect_warning(resi(mod.coxph, data = data.surv, nboot = 1, vcovfunc = sandwich::vcovHC),
                 "vcovfunc argument ignored for coxph objects")}
})


test_that("boot.method = 'bayes' only works for lm and nls", {
  expect_true(resi(mod.lm, nboot = 1, boot.method = "bayes")$boot.method == "bayes")
  expect_true(resi(mod.nls, data = data.nls, nboot = 1, boot.method = "bayes")$boot.method == "bayes")
  expect_error(resi(mod, nboot = 1, boot.method = "bayes"))
  if(requireNamespace("survival")){
    expect_error(resi(mod.surv, data = data.surv, nboot = 1, boot.method = "bayes"))
    expect_error(resi(mod.coxph, data = data.surv, nboot = 1, boot.method = "bayes"))}
  if(requireNamespace("pscl")){
    expect_error(resi(mod.hurdle, nboot = 1, boot.method = "bayes"))
    expect_error(resi(mod.zinf, nboot = 1, boot.method = "bayes"))}
})

test_that("data is needed for certain models", {
  expect_error(resi(mod.s, nboot = 1))
  expect_error(resi(mod.lm.s, nboot = 1))
  expect_error(resi(mod.nls, nboot = 1))
  if(requireNamespace("survival")){
  expect_error(resi(mod.surv, nboot = 1))
  expect_error(resi(mod.coxph, nboot = 1))}
  if(requireNamespace("gee")){
  expect_error(resi(mod.gee))}
})

# using different boot.results storage now
# test_that("boot.results stores correctly",{
#   expect_equal(colnames(resi(mod.log, nboot = 2, store.boot = TRUE)$boot.results$t), c("Overall", "(Intercept)", "age", "regionnorthwest", "regionsoutheast", "regionsouthwest", "age", "region"))
# })

test_that("resi produces the correct estimates", {
  expect_equal(unname(resi(mod, nboot = 1)$estimates), c(0.35982590, -0.06733537, -0.02670248, -0.03341748, -0.00246893, 0.14980504, 0.15238719,
                                                         0.05839381, 0.01667149, 0.03610640, -0.01497171, 0.03655364, 0.29533571, 0.14991488,
                                                         0.05159896, 0.01617303), tolerance = 1e-07)
  expect_equal(unname(resi(mod.lm, nboot = 1)$estimates), c(0.3595407548, -0.0672973369,
                                                            -0.0266873974, -0.0333986021,
                                                            -0.0024675352, 0.1497204142,
                                                            0.1523011075, 0.0583608197,
                                                            0.0166620711, 0.0360860037,
                                                            -0.0149632506, 0.0364798507,
                                                            0.2951113263, 0.1497981921,
                                                            0.0515491712, 0.0160560332), tolerance = 1e-07)
  expect_equal(unname(resi(mod.s, nboot = 1, data = data)$estimates[c(1, 6, 11, 24)]),
               c(0.35634034, 0.06143486, 0.00748577, 0.05339006), tolerance = 1e-07)
  if(requireNamespace("survival")){
  expect_equal(unname(resi(mod.surv, nboot = 1, data = data.surv)$estimates),
               c(0.194002944, 0.515785774, -0.080800275, 0.200549139, 0.105353954, -0.276392759,
                 0.046080344, 0.189247644, 0.081817903), tolerance = 1e-07)
  expect_equal(unname(resi(mod.coxph, nboot = 1, data = data.surv)$estimates),
               c(0.224354226, 0.136138740, -0.213293162, 0.008641209, 0.117732151,
                 0.202042262, 0.000000000))
  }
  if(requireNamespace("pscl")){
  expect_equal(unname(resi(mod.hurdle, nboot = 1)$estimates[c(1, 2, 5, 8, 10)]),
               c(0.28063439, 0.12770489, -0.06898884, 0.02686305, 0.06039791),
               tolerance = 1e-07)
  expect_equal(unname(resi(mod.zinf, nboot = 1)$estimates[c(1, 2, 5, 8, 10)]),
               c(0.23725614, 0.11892936, -0.07014934, -0.03481673, -0.03461235),
               tolerance = 1e-05)
  }
  if(requireNamespace("gee")){
  expect_equal(unname(resi(mod.gee, nboot = 10, data = data.gee)$coefficients[,'L-RESI']),
               c(-0.02585617, -0.48811210, 0.01414410, 0.56177861, -0.29398258),
               tolerance = 1e-07)
  expect_equal(unname(resi(mod.gee, nboot = 10, data = data.gee)$coefficients[,'CS-RESI']),
               c(-0.007006305,-0.121182453,0.003769394,0.119402133,-0.069223594),
               tolerance = 1e-07)
  }
  if(requireNamespace("geepack")){
  expect_equal(unname(resi(mod.geeglm, nboot = 10, data = data.gee)$coefficients[,'L-RESI']),
               c(-0.02585617, -0.48811210, 0.01414410, 0.56177861, -0.29398258),
               tolerance = 1e-07)
  expect_equal(unname(resi(mod.geeglm, nboot = 10, data = data.gee)$coefficients[,'CS-RESI']),
               c(-0.007006305,-0.121182453,0.003769394,0.119402133,-0.069223594),
               tolerance = 1e-07)
  }
  expect_equal(unname(suppressWarnings(resi(mod.lme, nboot = 10)$coefficients[,'RESI'])),
               c(3.659090, 1.739166, 0.512371), tolerance = 1e-07)
  if(requireNamespace("lme4")){
  expect_equal(unname(suppressWarnings(resi(mod.lmerMod, nboot = 10)$coefficients[,'RESI'])),
               c(8.434942, 1.533073), tolerance = 1e-07)}
})

test_that("RESI estimates are in between the confidence limits", {
  resi.obj = resi(mod, nboot = 500, store.boot = TRUE)
  expect_true(all(resi.obj$coefficients$RESI >= resi.obj$coefficients$`2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`97.5%`))
  expect_true(all(resi.obj$anova$RESI >= resi.obj$anova$`2.5%`) & all(resi.obj$anova$RESI <= resi.obj$anova$`97.5%`))
  An.obj = car::Anova(resi.obj, alpha = 0.01)
  expect_true(all(An.obj$RESI >= An.obj$`0.05%`) & all(An.obj$RESI <= An.obj$`99.5%`))
  resi.obj = resi(mod.lm, nboot = 500, store.boot = TRUE)
  expect_true(all(resi.obj$coefficients$RESI >= resi.obj$coefficients$`2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`97.5%`))
  expect_true(all(resi.obj$anova[1:5, "RESI"] >= resi.obj$anova[1:5, "2.5%"]) & all(resi.obj$anova[1:5, "RESI"] <= resi.obj$anova[1:5, "RESI"]))
  An.obj = car::Anova(resi.obj, alpha = 0.01)
  expect_true(all(An.obj$RESI[1:5] >= An.obj$`0.05%`[1:5]) & all(An.obj$RESI[1:5] <= An.obj$`99.5%`[1:5]))
  resi.obj = resi(mod.nls, data = data.nls)
  expect_true(all(resi.obj$coefficients$RESI >= resi.obj$coefficients$`2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`97.5%`))
  if(requireNamespace("survival")){
  resi.obj = resi(mod.coxph, data = data.surv, nboot = 500)
  expect_true(all(resi.obj$coefficients$RESI >= resi.obj$coefficients$`2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`97.5%`))
  expect_true(all(resi.obj$anova$RESI >= resi.obj$anova$`2.5%`) & all(resi.obj$anova$RESI <= resi.obj$anova$`97.5%`))}
  # commented for time
  # resi.obj = resi(mod.hurdle, nboot = 500)
  # expect_true(all(resi.obj$coefficients$RESI >= resi.obj$coefficients$`2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`97.5%`))
  if(requireNamespace("geepack")){
  resi.obj = resi(mod.geeglm, nboot = 500)
  expect_true(all(resi.obj$coefficients$`L-RESI` >= resi.obj$coefficients$`L 2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`L 97.5%`))
  expect_true(all(resi.obj$coefficients$`CS-RESI` >= resi.obj$coefficients$`CS 2.5%`) & all(resi.obj$coefficients$RESI <= resi.obj$coefficients$`CS 97.5%`))
  expect_true(all(resi.obj$anova$`L-RESI` >= resi.obj$anova$`L 2.5%`) & all(resi.obj$anova$RESI <= resi.obj$anova$`L 97.5%`))
  expect_true(all(resi.obj$anova$`CS-RESI` >= resi.obj$anova$`CS 2.5%`) & all(resi.obj$anova$RESI <= resi.obj$anova$`CS 97.5%`))}
})

if(requireNamespace("gee") & requireNamespace("geepack")){
test_that("vcov methods same for gee and geeglm",{
  expect_equal(unname(mod.gee$robust.variance), mod.geeglm$geese$vbeta)
  modg1 = glm(formula = formula(mod.gee), family = mod.gee$family, data = data.gee,
              contrasts = mod.gee$contrasts)
  modg2 = glm(formula = formula(mod.geeglm), family = mod.geeglm$family, data = data.gee,
              contrasts = mod.geeglm$contrasts)
  expect_equal(sandwich::vcovHC(modg1, type = "HC0"), sandwich::vcovHC(modg2, type = "HC0"))
})

test_that("boot.results same (approx) for gee and geeglm",{
  set.seed(123)
  resi.obj = resi(mod.geeglm, nboot = 5, store.boot = T, anova = F)
  set.seed(123)
  resi.obj2 = resi(mod.gee, data = data.gee, nboot = 5, store.boot = T)
  expect_equal(resi.obj$boot.results$t[,-1], resi.obj2$boot.results$t, tolerance = 1e-09)
})}

test_that("unbiased = FALSE returns same abs. RESI as Chi-sq/F",{
  resi.obj = resi(mod, nboot = 1, unbiased = FALSE)
  expect_equal(abs(resi.obj$coefficients[6:7, 'RESI']), resi.obj$anova[3:4, 'RESI'], tolerance = 1e-03)
  if(requireNamespace("survival")){
  resi.obj = resi(mod.surv, data = data.surv, nboot = 1, unbiased = FALSE)
  expect_equal(abs(resi.obj$coefficients[2:4, 'RESI']), resi.obj$anova[1:3, 'RESI'], tolerance = 1e-07)
  resi.obj = resi(mod.coxph, data = data.surv, nboot = 1, unbiased = FALSE)
  expect_equal(abs(resi.obj$coefficients[1:3, 'RESI']), resi.obj$anova[1:3, 'RESI'], tolerance = 1e-07)}
  resi.obj = resi(mod.lm, nboot = 1, unbiased = FALSE)
  expect_equal(abs(resi.obj$coefficients[6:7,'RESI']), resi.obj$anova[3:4,'RESI'], tolerance = 1e-07)
})

test_that("Using a matrix for data works", {
  expect_true(class(resi(mod.nls, nboot = 10, data = as.matrix(data.nls))) == "resi")
})

test_that("Specifying a reduced model only changes overall output", {
  resi.obj = resi(mod, nboot = 1)
  resi.objr = resi(mod, model.reduced = mod.r, nboot = 1)
  expect_equal(resi.obj$coefficients$RESI, resi.objr$coefficients$RESI)
  expect_equal(resi.obj$anova$RESI, resi.objr$anova$RESI)
  expect_false(resi.obj$overall$RESI[2] == resi.objr$overall$RESI[2])
  resi.obj = resi(mod.lm.s, data = data, nboot = 1)
  resi.objr = resi(mod.lm.s, model.reduced = mod.lm.s.r, data = data, nboot = 1)
  expect_equal(resi.obj$coefficients$RESI, resi.objr$coefficients$RESI)
  expect_equal(resi.obj$anova$RESI, resi.objr$anova$RESI)
  expect_false(resi.obj$overall$RESI[2] == resi.objr$overall$RESI[2])
  if(requireNamespace("survival")){
  resi.obj = resi(mod.surv, data = data.surv, nboot = 1)
  resi.objr = resi(mod.surv, model.reduced = mod.surv.red, data = data.surv[which(!(1:nrow(data.surv)%in% mod.surv$na.action)),], nboot = 10)
  expect_false(resi.obj$overall$RESI[2] == resi.objr$overall$RESI[2])}
  if(requireNamespace("pscl")){
  resi.obj = resi(mod.hurdle, nboot = 1)
  resi.objr = resi(mod.hurdle, model.reduced = mod.hurdle.r, nboot = 1)
  expect_equal(resi.obj$coefficients$RESI, resi.objr$coefficients$RESI)
  expect_equal(resi.obj$anova$RESI, resi.objr$anova$RESI)
  expect_false(resi.obj$overall$RESI[2] == resi.objr$overall$RESI[2])}
})

test_that("specifying additional vcov args works",{
  expect_false(any(resi(mod, nboot = 1, vcov.args = list(type = "HC0"))$coefficients[,'RESI'] == resi(mod, nboot = 1)$coefficients[,'RESI']))
})

test_that("specifying additional Anova args works",{
  expect_false(all(resi(mod, nboot = 1, Anova.args = list(type = "3"))$anova[2:6,'RESI'] == resi(mod, nboot = 1)$anova[,'RESI']))
  if(requireNamespace("survival")){
  expect_false(length(resi(mod.surv, data = data.surv, nboot = 1, Anova.args = list(type = "3"))$anova$RESI) == length(resi(mod.surv, data = data.surv, nboot = 1)$anova$RESI))}
})

test_that("vcovfunc = vcov changes naive.var to TRUE",{
  expect_true(resi(mod, vcovfunc = vcov, nboot = 10)$naive.var == TRUE)
})

if(requireNamespace("tibble")){
test_that("tibbles work", {
  expect_equal(unname(resi(mod.glm.tib, nboot = 1)$estimates), c(0.35982590, -0.06733537, -0.02670248, -0.03341748, -0.00246893, 0.14980504, 0.15238719,
                                                                 0.05839381, 0.01667149, 0.03610640, -0.01497171, 0.03655364, 0.29533571, 0.14991488,
                                                                 0.05159896, 0.01617303), tolerance = 1e-07)
  expect_equal(unname(resi(mod.lm.tib, nboot = 1)$estimates), c(0.3595407548, -0.0672973369,
                                                                -0.0266873974, -0.0333986021,
                                                                -0.0024675352, 0.1497204142,
                                                                0.1523011075, 0.0583608197,
                                                                0.0166620711, 0.0360860037,
                                                                -0.0149632506, 0.0364798507,
                                                                0.2951113263, 0.1497981921,
                                                                0.0515491712, 0.0160560332))
  if(requireNamespace("survival")){
  expect_equal(unname(resi(mod.surv.tib, nboot = 1, data = data.surv)$estimates), c(0.194002944, 0.515785774, -0.080800275, 0.200549139, 0.105353954, -0.276392759,
                                                                                    0.046080344, 0.189247644, 0.081817903), tolerance = 1e-07)
  expect_equal(unname(resi(mod.coxph.tib, nboot = 1, data = data.surv)$estimates), c(0.224354226, 0.136138740, -0.213293162, 0.008641209, 0.117732151, 0.202042262, 0.000000000))
  }
  if(requireNamespace("pscl")){
  expect_equal(unname(resi(mod.hurdle.tib, nboot = 1)$estimates[c(1, 2, 5, 8, 10)]), c(0.28063439, 0.12770489, -0.06898884, 0.02686305, 0.06039791), tolerance = 1e-07)
  expect_equal(unname(resi(mod.zinf.tib, nboot = 1)$estimates[c(1, 2, 5, 8, 10)]), c(0.23725614, 0.11892936, -0.07014934, -0.03481673, -0.03461235), tolerance = 1e-05)
  }
  if(requireNamespace("gee")){
  expect_equal(unname(resi(mod.gee.tib, nboot = 10, data = data.gee)$coefficients[,'L-RESI']), c(-0.02585617, -0.48811210, 0.01414410, 0.56177861, -0.29398258), tolerance = 1e-07)
  }
  if(requireNamespace("geepack")){
  expect_equal(unname(resi(mod.geeglm.tib, nboot = 10, data = data.gee)$coefficients[,'L-RESI']), c(-0.02585617, -0.48811210, 0.01414410, 0.56177861, -0.29398258), tolerance = 1e-07)
  }
  if(requireNamespace("lme4")){
  expect_equal(unname(suppressWarnings(resi(mod.lmerMod.tib, nboot = 10))$coefficients[,'RESI']),c(8.434942, 1.533073), tolerance = 1e-07)
  }
})}

test_that("intercept-only full model works and doesn't have overall element", {
  expect_true(is.null(resi(mod.glm.int, nboot = 1, Anova.args = list(type = 3))$overall))
  expect_true(is.null(resi(mod.lm.int, nboot = 1, Anova.args = list(type = 3))$overall))
  if(requireNamespace("survival")){
    expect_true(is.null(resi(mod.surv.int, data = data.surv, nboot = 1, Anova.args = list(type = 3))$overall))
  }
  if(requireNamespace("pscl")){
    expect_true(is.null(resi(mod.hurdle.int, nboot = 1, Anova.args = list(type = 3))$overall))
    expect_true(is.null(resi(mod.zinf.int, nboot = 1, Anova.args = list(type = 3))$overall))
  }
  # to be implemented
  # if(requireNamespace("geepack")){
  #   expect_true(is.null(resi(mod.geeglm.int, nboot = 1)$overall))
  # }

})

test_that("same model.full and model.reduced produces error",{
  expect_error(resi(model.full = mod, model.reduced = mod, nboot = 1))
})

# important because some resi_pe methods use waldtest and some use wald.test
test_that("wald.test and waldtest still consistent",{
  expect_equal(unname(lmtest::waldtest(mod.lm, vcov = vcovHC, test = "Chisq")$Chisq[2]),
               unname(aod::wald.test(vcovHC(mod.lm), coef(mod.lm),
                                     Terms = 2:nrow(vcovHC(mod.lm)))$result$chi2["chi2"]))
})

# issue from github
test_that("overall = FALSE doesn't produce an error", {
  expect_equal(unname(resi(mod.lm, nboot = 1, overall = F, coefficients = F)$estimates),
               c(0.03647985, 0.29511133, 0.14979819, 0.05154917, 0.01605603), tolerance = 1e-05)
})


