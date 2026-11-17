egsingle_ci_data <- function() {
  egsingle <- local({
    env <- new.env()
    utils::data("egsingle", package = "mlmRev", envir = env)
    env$egsingle
  })
  egsingle$size_100 <- egsingle$size / 100
  egsingle$lowinc_10 <- egsingle$lowinc / 10
  egsingle
}
form_ci <- math ~ year + female + black + hispanic + size_100 + lowinc_10 +
  (1 | schoolid/childid)

test_that("fast CR2 equals clubSandwich CR2 and GLS pieces reproduce lme4", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("mlmRev")
  skip_if_not_installed("clubSandwich")
  eg <- egsingle_ci_data()
  sub <- droplevels(subset(eg, schoolid %in% levels(schoolid)[1:12]))
  fit <- lme4::lmer(form_ci, data = sub, REML = TRUE)
  expect_equal(unname(.mixedCiCR2(fit, "schoolid")),
               unname(as.matrix(clubSandwich::vcovCR(fit, cluster = sub$schoolid, type = "CR2"))),
               tolerance = 1e-8)

  base <- .mixedCiBase(fit, "schoolid")
  gls <- .mixedCiGls(base, base$theta, base$sigma)
  expect_equal(unname(.mixedCiDeviance(base, c(base$theta, base$sigma))), unname(lme4::REMLcrit(fit)))
  expect_equal(unname(solve(gls$info)), unname(as.matrix(stats::vcov(fit))), tolerance = 1e-6)
  pieces <- .mixedCiClusterPieces(base, gls)
  expect_equal(unname(matrix(colSums(pieces$I_s), ncol(base$X))), unname(gls$info))
  expect_equal(unname(colSums(pieces$U_s)), rep(0, ncol(base$X)), tolerance = 1e-6)
})

test_that("fast CR2 matches clubSandwich for random slopes and singular fits", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("clubSandwich")
  ss <- lme4::sleepstudy
  slope <- lme4::lmer(Reaction ~ Days + (Days | Subject), data = ss)
  expect_equal(unname(.mixedCiCR2(slope, "Subject")),
               unname(as.matrix(clubSandwich::vcovCR(slope, cluster = ss$Subject, type = "CR2"))),
               tolerance = 1e-8)
  set.seed(1)
  d <- data.frame(g = factor(rep(1:30, each = 4)), x = rnorm(120))
  d$y <- d$x + rnorm(120)
  sing <- suppressMessages(lme4::lmer(y ~ x + (1 | g), data = d))
  skip_if_not(lme4::isSingular(sing))
  expect_equal(unname(.mixedCiCR2(sing, "g")),
               unname(as.matrix(clubSandwich::vcovCR(sing, cluster = d$g, type = "CR2"))),
               tolerance = 1e-8)
})

test_that("design-variance linearization matches resampling of school profiles", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("mlmRev")
  eg <- egsingle_ci_data()
  fit <- lme4::lmer(form_ci, data = eg, REML = FALSE)
  base <- .mixedCiBase(fit, "schoolid")
  gls <- .mixedCiGls(base, base$theta, base$sigma)
  pieces <- .mixedCiClusterPieces(base, gls)
  p <- ncol(base$X)
  counts <- tapply(lme4::getME(fit, "flist")[["childid:schoolid"]], base$cluster,
                   function(x) length(unique(x)))
  L <- matrix(colnames(base$X) == "year", 1L)
  Vhat <- solve(gls$info)
  design_sd <- function(n_k, cnt) {
    with_design <- .mixedCiContrast(L, n_k, gls$beta, Vhat, gls$info, pieces, FALSE, TRUE, cnt)
    without <- .mixedCiContrast(L, n_k, gls$beta, Vhat, gls$info, pieces, FALSE, FALSE, cnt)
    sqrt(drop(with_design$var_R - without$var_R))
  }
  set.seed(1)
  g <- nrow(pieces$I_s)
  k <- which(L[1, ])
  boot <- replicate(2000, {
    s <- sample.int(g, g, replace = TRUE)
    Istar <- matrix(colSums(pieces$I_s[s, , drop = FALSE]), p)
    gls$beta[k] / sqrt(sum(counts[s]) * solve(Istar)[k, k])
  })
  expect_equal(design_sd(sum(counts), counts), stats::sd(boot), tolerance = 0.1)
  # Omitting the unit-count perturbation grossly overstates the design SD.
  expect_gt(design_sd(sum(counts), NULL), 2 * stats::sd(boot))
})

test_that("resi() on lmerMod returns ordered RESI intervals for each method", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("mlmRev")
  skip_if_not_installed("clubSandwich")
  skip_if_not_installed("numDeriv")
  eg <- egsingle_ci_data()
  sub <- droplevels(subset(eg, schoolid %in% levels(schoolid)[1:20]))
  fit <- lme4::lmer(form_ci, data = sub, REML = TRUE)
  check <- function(out) {
    for (tab in list(out$coefficients, out$anova)) {
      expect_true(all(c("2.5%", "97.5%") %in% names(tab)))
      expect_true(all(tab[["2.5%"]] <= tab$RESI & tab$RESI <= tab[["97.5%"]]))
    }
  }
  for (method in c("normal", "qf", "cf")) {
    check(resi(fit, vcovfunc = stats::vcov, ci.method = method, level = 2))
  }
  robust <- resi(fit, vcov.args = list(type = "CR2"), cluster = "schoolid",
                 ci.method = "normal", level = 2)
  check(robust)
  # CR2 intervals use the t(g - 1) reference: signed normal widths scale by
  # qt / qnorm relative to the normal reference.
  Vcr2 <- .mixedCiCR2(fit, "schoolid")
  n_coef <- setNames(robust$coefficients$n, rownames(robust$coefficients))
  ci_z <- .resi_ci_mixed(fit, Vcr2, TRUE, n_coef, ci.method = "normal",
                         cluster = "schoolid", ref.df = Inf)$normal$coefficients
  expect_equal(robust$coefficients[["97.5%"]] - robust$coefficients[["2.5%"]],
               unname(ci_z[, 2] - ci_z[, 1]) * qt(0.975, 19) / qnorm(0.975),
               tolerance = 1e-6)
  expect_equal(robust$coefficients$Robust.SE,
               unname(sqrt(diag(.mixedCiCR2(fit, "schoolid")))))
  expect_warning(none <- resi(fit, vcovfunc = stats::vcov, ci.method = "none", level = 2),
                 "point estimates only")
  expect_false("2.5%" %in% names(none$coefficients))
})

test_that("resi() on lme reproduces the lmerMod intervals and falls back when unsupported", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("mlmRev")
  skip_if_not_installed("clubSandwich")
  skip_if_not_installed("numDeriv")
  eg <- egsingle_ci_data()
  sub <- droplevels(subset(eg, schoolid %in% levels(schoolid)[1:20]))
  sub$math[c(5, 50, 500)] <- NA
  fit_lmer <- lme4::lmer(form_ci, data = sub, REML = TRUE)
  fit_lme <- nlme::lme(math ~ year + female + black + hispanic + size_100 + lowinc_10,
                       random = ~ 1 | schoolid/childid, data = sub, na.action = na.omit)
  base <- .mixedCiBase(fit_lme, "schoolid")
  gls <- .mixedCiGls(base, base$theta, base$sigma)
  expect_equal(unname(gls$beta), unname(nlme::fixef(fit_lme)), tolerance = 1e-8)
  expect_equal(.mixedCiDeviance(base, c(base$theta, base$sigma)),
               -2 * as.numeric(stats::logLik(fit_lme)), tolerance = 1e-8)
  cols <- c("RESI", "2.5%", "97.5%")
  for (args in list(list(vcovfunc = stats::vcov),
                    list(vcov.args = list(type = "CR2"), cluster = "schoolid"))) {
    a <- suppressMessages(do.call(resi, c(list(fit_lme, ci.method = "qf", level = 2), args)))
    b <- suppressMessages(do.call(resi, c(list(fit_lmer, ci.method = "qf", level = 2), args)))
    expect_equal(unname(as.matrix(a$coefficients[cols])), unname(as.matrix(b$coefficients[cols])),
                 tolerance = 1e-4)
    expect_equal(unname(as.matrix(as.data.frame(a$anova)[cols])),
                 unname(as.matrix(as.data.frame(b$anova)[cols])), tolerance = 1e-4)
  }
  slope <- nlme::lme(Reaction ~ Days, random = ~ Days | Subject, data = lme4::sleepstudy)
  expect_equal(unname(.mixedCiCR2(slope, "Subject")),
               unname(as.matrix(clubSandwich::vcovCR(slope, cluster = lme4::sleepstudy$Subject,
                                                     type = "CR2"))),
               tolerance = 1e-8)
  hetero <- nlme::lme(Reaction ~ Days, random = ~ 1 | Subject, data = lme4::sleepstudy,
                      weights = nlme::varIdent(form = ~ 1 | Days > 4))
  expect_warning(out <- suppressMessages(resi(hetero, vcovfunc = stats::vcov)),
                 "not available for this lme fit")
  expect_false("2.5%" %in% names(out$coefficients))
})
