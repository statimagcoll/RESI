test_that("resi_pe.gam converts mgcv Wald statistics", {
  skip_if_not_installed("mgcv")
  ins <- RESI::insurance
  ins$region <- factor(ins$region)
  b <- mgcv::gam(log10(charges) ~ smoker + region + s(age) + s(bmi),
                 data = ins, method = "REML")
  s <- mgcv::summary.gam(b)
  r <- resi_pe(b)
  n <- nrow(ins)
  expect_equal(r$coefficients["smoker", "RESI"],
               t2S(s$p.table["smoker", "t value"], b$df.residual, n))
  expect_equal(r$anova["s(age)", "RESI"],
               f2S(s$s.table["s(age)", "F"], s$s.table["s(age)", "Ref.df"],
                   b$df.residual, n))
  expect_equal(r$anova["region", "RESI"],
               f2S(s$pTerms.table["region", "F"], 3, b$df.residual, n))

  bb <- mgcv::gam(I(charges > 15000) ~ smoker + s(age), family = binomial,
                  data = ins)
  sb <- mgcv::summary.gam(bb)
  rb <- resi_pe(bb)
  expect_equal(rb$anova["s(age)", "RESI"],
               chisq2S(sb$s.table["s(age)", "Chi.sq"],
                       sb$s.table["s(age)", "Ref.df"], n))
  expect_true("z value" %in% colnames(rb$coefficients))
})

test_that("resi_pe.gam returns NA for smooths without a null space", {
  skip_if_not_installed("mgcv")
  ins <- RESI::insurance
  ins$region <- factor(ins$region)
  b <- mgcv::bam(log10(charges) ~ smoker + s(age) + s(region, bs = "re"),
                 data = ins)
  expect_warning(r <- resi_pe(b), "null space")
  expect_true(is.na(r$anova["s(region)", "RESI"]))
  expect_false(is.na(r$anova["s(age)", "RESI"]))
  expect_error(resi(b), "resi_pe")
})
