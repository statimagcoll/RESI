## Tests for asymptotic CI functions (resi_pe_asymptotic, resi ci.method)
## These are fast (no bootstrap) tests checking mathematical properties.

library(sandwich)
library(car)

data <- RESI::insurance

# Fit a few standard models used throughout
mod.lm  <- lm(charges ~ region * age + bmi + sex, data = data)
mod.glm <- glm(smoker ~ age + region, data = data, family = "binomial")


# ---------------------------------------------------------------------------
# Helper: run both ci.method values and return the output
# ---------------------------------------------------------------------------
asym_n <- resi_pe_asymptotic(mod.lm,  ci.method = "normal")
asym_q <- resi_pe_asymptotic(mod.lm,  ci.method = "qf")
asym_g <- resi_pe_asymptotic(mod.glm, ci.method = "normal")


# ===========================================================================
# 1. Output structure
# ===========================================================================
test_that("resi_pe_asymptotic returns expected structure", {
  expect_s3_class(asym_n, "resi")

  # Both tables present
  expect_true(!is.null(asym_n$coefficients))
  expect_true(!is.null(asym_n$anova))

  # CI columns named correctly for default alpha = 0.05
  ci_cols <- c("2.5%", "97.5%")
  expect_true(all(ci_cols %in% colnames(asym_n$coefficients)))
  expect_true(all(ci_cols %in% colnames(asym_n$anova)))
  expect_true(all(ci_cols %in% colnames(asym_q$anova)))

  # RESI column present
  expect_true("RESI" %in% colnames(asym_n$coefficients))
  expect_true("RESI" %in% colnames(asym_n$anova))

  # custom alpha produces correct column names
  out_01 <- resi_pe_asymptotic(mod.lm, alpha = 0.01, ci.method = "normal")
  expect_true(all(c("0.5%", "99.5%") %in% colnames(out_01$coefficients)))

  # GLM also works
  expect_s3_class(asym_g, "resi")
  expect_true(!is.null(asym_g$anova))
})


# ===========================================================================
# 2. Point estimates match resi_pe()
# ===========================================================================
test_that("resi_pe_asymptotic RESI point estimates match resi_pe()", {
  pe_lm  <- resi_pe(mod.lm)
  pe_glm <- resi_pe(mod.glm)

  # lm: coefficients
  expect_equal(asym_n$coefficients[, "RESI"],
               pe_lm$coefficients[, "RESI"],
               tolerance = 1e-6,
               label = "lm coefficients RESI")

  # lm: anova (compare on non-Residuals rows)
  anova_rows <- rownames(asym_n$anova)
  expect_equal(asym_n$anova[anova_rows, "RESI"],
               pe_lm$anova[anova_rows, "RESI"],
               tolerance = 1e-6,
               label = "lm anova RESI")

  # glm: anova
  anova_rows_g <- rownames(asym_g$anova)
  expect_equal(asym_g$anova[anova_rows_g, "RESI"],
               pe_glm$anova[anova_rows_g, "RESI"],
               tolerance = 1e-6,
               label = "glm anova RESI")
})


# ===========================================================================
# 3. RESI is contained within its own CI (unsigned / anova)
# ===========================================================================
test_that("RESI point estimate is within asymptotic CIs (anova, unsigned)", {
  # normal method
  resi_col  <- asym_n$anova[, "RESI"]
  lo_col    <- asym_n$anova[, "2.5%"]
  hi_col    <- asym_n$anova[, "97.5%"]
  expect_true(all(lo_col <= resi_col + 1e-10),
              label = "normal: LCI <= RESI (anova)")
  expect_true(all(resi_col <= hi_col + 1e-10),
              label = "normal: RESI <= UCI (anova)")

  # qf method
  resi_col_q <- asym_q$anova[, "RESI"]
  lo_col_q   <- asym_q$anova[, "2.5%"]
  hi_col_q   <- asym_q$anova[, "97.5%"]
  expect_true(all(lo_col_q <= resi_col_q + 1e-10),
              label = "qf: LCI <= RESI (anova)")
  expect_true(all(resi_col_q <= hi_col_q + 1e-10),
              label = "qf: RESI <= UCI (anova)")
})


# ===========================================================================
# 4. Non-negativity of unsigned CIs (anova table)
# ===========================================================================
test_that("Asymptotic anova CIs are non-negative", {
  expect_true(all(asym_n$anova[, "2.5%"]  >= -1e-10))
  expect_true(all(asym_n$anova[, "97.5%"] >  0))
  expect_true(all(asym_q$anova[, "2.5%"]  >= -1e-10))
  expect_true(all(asym_q$anova[, "97.5%"] >  0))
})


# ===========================================================================
# 5. Monotonicity: wider alpha gives narrower CI
# ===========================================================================
test_that("Wider alpha gives narrower asymptotic CI", {
  out_05 <- resi_pe_asymptotic(mod.lm, alpha = 0.05, ci.method = "normal")
  out_10 <- resi_pe_asymptotic(mod.lm, alpha = 0.10, ci.method = "normal")

  width_05_coef <- out_05$coefficients[, "97.5%"] - out_05$coefficients[, "2.5%"]
  width_10_coef <- out_10$coefficients[, "95%"]   - out_10$coefficients[, "5%"]
  expect_true(all(width_10_coef <= width_05_coef + 1e-10),
              label = "coefficients: alpha=0.10 width <= alpha=0.05 width")

  width_05_an <- out_05$anova[, "97.5%"] - out_05$anova[, "2.5%"]
  width_10_an <- out_10$anova[, "95%"]   - out_10$anova[, "5%"]
  expect_true(all(width_10_an <= width_05_an + 1e-10),
              label = "anova: alpha=0.10 width <= alpha=0.05 width")
})


# ===========================================================================
# 6. Normal and QF CIs have similar widths for lm (large n = 1338)
# ===========================================================================
test_that("Normal and QF anova CIs are close for large n (lm)", {
  # With n=1338, both methods should agree within ~20% of CI width
  width_n <- asym_n$anova[, "97.5%"] - asym_n$anova[, "2.5%"]
  width_q <- asym_q$anova[, "97.5%"] - asym_q$anova[, "2.5%"]

  # Relative difference in width should be < 30% for each term
  rel_diff <- abs(width_q - width_n) / pmax(width_n, 1e-6)
  expect_true(all(rel_diff < 0.30),
              info  = paste("Relative width differences:", round(rel_diff, 3)),
              label = "normal vs QF CI widths agree within 30%")
})


# ===========================================================================
# 7. CI bounds agree on coefficients: normal and QF both use signed normal
# ===========================================================================
test_that("QF and normal coefficient CIs are identical (both use signed normal)", {
  out_qf <- resi_pe_asymptotic(mod.lm, ci.method = "qf")
  # Coefficients CIs are always from the signed normal path regardless of ci.method
  expect_equal(asym_n$coefficients[, "2.5%"],
               out_qf$coefficients[, "2.5%"],
               tolerance = 1e-8,
               label = "LCI: normal == qf for coefficients")
  expect_equal(asym_n$coefficients[, "97.5%"],
               out_qf$coefficients[, "97.5%"],
               tolerance = 1e-8,
               label = "UCI: normal == qf for coefficients")
})


# ===========================================================================
# 8. resi() with ci.method="normal" / "qf" skips bootstrap and returns CIs
# ===========================================================================
test_that("resi() ci.method='normal' returns CIs without bootstrap", {
  # Should complete quickly (no boot) and have CI columns
  out_r <- resi(mod.lm, ci.method = "normal")
  expect_s3_class(out_r, "resi")
  expect_true("ci.method" %in% names(out_r))
  expect_equal(out_r$ci.method, "normal")
  expect_true(all(c("2.5%", "97.5%") %in% colnames(out_r$coefficients)))
  expect_true(all(c("2.5%", "97.5%") %in% colnames(out_r$anova)))

  out_rq <- resi(mod.lm, ci.method = "qf")
  expect_s3_class(out_rq, "resi")
  expect_equal(out_rq$ci.method, "qf")
})


# ===========================================================================
# 9. Correctly handles anova = FALSE or coefficients = FALSE
# ===========================================================================
test_that("resi_pe_asymptotic respects anova/coefficients flags", {
  out_no_an  <- resi_pe_asymptotic(mod.lm, anova = FALSE)
  out_no_co  <- resi_pe_asymptotic(mod.lm, coefficients = FALSE)

  expect_null(out_no_an$anova)
  expect_true(!is.null(out_no_an$coefficients))

  expect_null(out_no_co$coefficients)
  expect_true(!is.null(out_no_co$anova))
})


# ===========================================================================
# 10. Under a known-null predictor, QF lower bound should be 0
# ===========================================================================
test_that("QF LCI = 0 for near-null predictor", {
  # region:age interaction is near-null in lm (RESI ~ 0)
  out_q <- resi_pe_asymptotic(mod.lm, ci.method = "qf")
  interaction_lci <- out_q$anova["region:age", "2.5%"]
  expect_equal(interaction_lci, 0, tolerance = 1e-6,
               label = "region:age QF lower CI should be 0")
})


# ===========================================================================
# 11. Signed coefficient CIs: negative RESI can have negative LCI
# ===========================================================================
test_that("Signed coefficient CIs allow negative values", {
  resi_signed <- asym_n$coefficients[, "RESI"]
  lci_signed  <- asym_n$coefficients[, "2.5%"]

  # If any RESI is negative, its LCI should also be negative
  neg_resi <- resi_signed[resi_signed < 0]
  if (length(neg_resi) > 0) {
    neg_lci <- lci_signed[resi_signed < 0]
    expect_true(all(neg_lci < 0),
                label = "Negative RESI should have negative LCI")
  }
})
