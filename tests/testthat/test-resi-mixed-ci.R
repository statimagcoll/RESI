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
  expect_equal(robust$coefficients$Robust.SE,
               unname(sqrt(diag(.mixedCiCR2(fit, "schoolid")))))
  expect_warning(none <- resi(fit, vcovfunc = stats::vcov, ci.method = "none", level = 2),
                 "point estimates only")
  expect_false("2.5%" %in% names(none$coefficients))
})
