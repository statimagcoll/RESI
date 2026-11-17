test_that("mixed-model RESI normalizes between-school terms by g and within-school terms by N", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("mlmRev")

  egsingle <- local({
    env <- new.env()
    utils::data("egsingle", package = "mlmRev", envir = env)
    env$egsingle
  })
  egsingle$schoolid <- factor(egsingle$schoolid)
  egsingle$childid <- factor(egsingle$childid)
  egsingle$size_100 <- egsingle$size / 100
  egsingle$lowinc_10 <- egsingle$lowinc / 10

  fit <- lme4::lmer(
    math ~ year + female + black + hispanic + size_100 + lowinc_10 +
      (1 | schoolid/childid),
    data = egsingle,
    REML = FALSE
  )
  pe <- resi_pe(fit, anova = TRUE, vcovfunc = stats::vcov, level = 2)

  expect_true(all(c("RESI", "Term", "n", "Variation level", "n level") %in% names(pe$coefficients)))
  expect_false(any(c("L-RESI", "CS-RESI") %in% names(pe$coefficients)))
  between <- c("(Intercept)", "size_100", "lowinc_10")
  expected_n <- ifelse(rownames(pe$coefficients) %in% between, 60, 1721)
  expect_equal(pe$coefficients$n, expected_n)
  expect_equal(pe$coefficients$`n level`,
               ifelse(rownames(pe$coefficients) %in% between, "schoolid", "childid:schoolid"))
  expect_equal(pe$coefficients["year", "Variation level"], "observation within childid:schoolid")
  expect_equal(pe$coefficients["lowinc_10", "Variation level"], "schoolid")
  expect_true(all(c("RESI", "n", "Variation level", "n level") %in% names(pe$anova)))
  expect_equal(pe$anova$n, ifelse(rownames(pe$anova) %in% between, 60, 1721))

  prop_cols <- c("Var prop schoolid", "Var prop childid:schoolid")
  expect_true(all(prop_cols %in% names(pe$coefficients)))
  expect_true(all(prop_cols %in% names(pe$anova)))
  expect_true(all(is.na(pe$coefficients["(Intercept)", prop_cols])))
  expect_equal(unname(unlist(pe$coefficients["lowinc_10", prop_cols])), c(1, 0))
  year_props <- unname(unlist(pe$coefficients["year", prop_cols]))
  expect_true(all(year_props > 0) && sum(year_props) < 1)
  expect_equal(pe$anova["black", "Var prop schoolid"],
               pe$coefficients["black1", "Var prop schoolid"])

  pe_default <- expect_message(
    resi_pe(fit, anova = FALSE, vcovfunc = stats::vcov),
    "level 1")
  expect_true(all(pe_default$coefficients$n == 60))
  expect_true(all(pe_default$coefficients$`n level` == "schoolid"))
  pe_global <- suppressMessages(resi_pe(fit, anova = FALSE, n = 100,
                                        vcovfunc = stats::vcov))
  expect_true(all(pe_global$coefficients$n == 100))

  pe_no_level <- suppressWarnings(resi_pe(fit, anova = FALSE, include_level = FALSE))
  expect_false("Variation level" %in% names(pe_no_level$coefficients))
  expect_false("n level" %in% names(pe_no_level$coefficients))
  expect_true("n" %in% names(pe_no_level$coefficients))
})

test_that("nlme mixed-model RESI uses the same normalizer API", {
  skip_if_not_installed("nlme")

  fit <- nlme::lme(distance ~ age + Sex, data = nlme::Orthodont,
                   random = ~1 | Subject, method = "ML")
  pe <- suppressMessages(resi_pe(fit, anova = TRUE, vcovfunc = stats::vcov))

  expect_true(all(c("RESI", "n", "Variation level", "n level") %in% names(pe$coefficients)))
  expect_false(any(c("L-RESI", "CS-RESI") %in% names(pe$coefficients)))
  expect_true(all(pe$coefficients$n == 27))
  expect_equal(pe$coefficients["age", "Variation level"], "observation within Subject")
  expect_true(all(pe$coefficients$`n level` == "Subject"))
  expect_true(all(c("RESI", "n", "Variation level") %in% names(pe$anova)))
  expect_equal(pe$coefficients["age", "Var prop Subject"], 0)
  expect_equal(pe$anova["Sex", "Var prop Subject"], 1)
})

test_that("a random slope at a coarser level sets the term's n to that level's count", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("nlme")
  skip_if_not_installed("mlmRev")
  egsingle <- local({
    env <- new.env()
    utils::data("egsingle", package = "mlmRev", envir = env)
    env$egsingle
  })
  slope_school <- suppressWarnings(lme4::lmer(
    math ~ year + female + size + (year | schoolid) + (1 | childid:schoolid), data = egsingle))
  pe <- resi_pe(slope_school, vcovfunc = stats::vcov, level = 2)
  expect_equal(pe$coefficients[c("year", "femaleMale"), "n"], c(60, 1721))
  expect_equal(pe$coefficients["year", "n level"], "schoolid")
  expect_equal(pe$anova["year", "n"], 60)
  expect_equal(pe$coefficients["year", "Variation level"], "observation within childid:schoolid")

  slope_child <- lme4::lmer(
    math ~ year + female + size + (1 | schoolid) + (year | childid:schoolid), data = egsingle)
  expect_equal(resi_pe(slope_child, vcovfunc = stats::vcov, level = 2)$coefficients["year", "n"], 1721)

  slope_lme <- nlme::lme(math ~ year + female + size,
                         random = list(schoolid = ~ year, childid = ~ 1), data = egsingle)
  expect_equal(resi_pe(slope_lme, vcovfunc = stats::vcov, level = 2)$coefficients["year", "n"], 60)
})
