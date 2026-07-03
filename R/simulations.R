# [AI-generated: Claude Sonnet 4.6, 2026-06-15]
# [Prompt summary: Plasmode simulation pipeline for Goals 1 and 2 - RESI bootstrap CI evaluation]

# ============================================================
#  Internal helpers
# ============================================================

# Extract aligned RESI point estimates and CI bounds from one resi object
# for a specified set of row names (terms).
.simExtractRows <- function(resi_obj, table_name, terms, ci_lo, ci_hi) {
  tab <- resi_obj[[table_name]]
  if (is.null(tab)) {
    na_vec <- setNames(rep(NA_real_, length(terms)), terms)
    return(list(resi = na_vec, lo = na_vec, hi = na_vec))
  }
  idx <- match(terms, rownames(tab))
  list(
    resi = setNames(tab[idx, "RESI"],  terms),
    lo   = setNames(if (ci_lo %in% colnames(tab)) tab[idx, ci_lo] else
                      rep(NA_real_, length(terms)), terms),
    hi   = setNames(if (ci_hi %in% colnames(tab)) tab[idx, ci_hi] else
                      rep(NA_real_, length(terms)), terms)
  )
}

# Overwrite RESI column in a resi_pe object with direct (population) formulas:
#   T-stat  -> S = t / sqrt(n)          (signed)
#   Z-stat  -> S = z / sqrt(n)          (signed)
#   F-stat  -> S = sqrt(F * df / n)     (no rdf or df correction)
#   Chisq   -> S = sqrt(chisq / n)      (no df subtraction)
# This treats the full-dataset test statistics as the true population parameters
# rather than using small-sample-corrected estimators from resi_pe.
.simDirectRESI <- function(pe_obj, n, rdf = NULL) {
  if (!is.null(pe_obj$coefficients)) {
    tab <- pe_obj$coefficients
    if ("t value" %in% colnames(tab)) {
      # For parametric lm, rdf = n-p is passed so that S_true = t/sqrt(rdf),
      # which targets beta/sqrt(diag(n*(X'X)^{-1}*mean(resid^2))) -- the
      # population RESI using sigma^2 = mean(resid^2) with no df correction.
      denom <- if (!is.null(rdf)) sqrt(rdf) else sqrt(n)
      tab[, "RESI"] <- tab[, "t value"] / denom
    } else if ("z value" %in% colnames(tab)) {
      tab[, "RESI"] <- tab[, "z value"] / sqrt(n)
    }
    pe_obj$coefficients <- tab
  }
  if (!is.null(pe_obj$anova)) {
    tab <- pe_obj$anova
    if ("F" %in% colnames(tab)) {
      # Same logic: S_true = sqrt(F * df / rdf) when rdf is supplied.
      denom <- if (!is.null(rdf)) rdf else n
      tab[, "RESI"] <- sqrt(tab[, "F"] * tab[, "Df"] / denom)
    } else if ("Chisq" %in% colnames(tab)) {
      tab[, "RESI"] <- sqrt(tab[, "Chisq"] / n)
    }
    pe_obj$anova <- tab
  }
  pe_obj
}

# Compute per-term simulation metrics from a list of resi objects.
# Returns a data.frame with rows = terms and columns = bias, mse, coverage, width.
.simComputeMetrics <- function(reps, table_name, true_resi, ci_lo, ci_hi,
                                exclude_rows = NULL) {
  terms <- rownames(true_resi)
  if (!is.null(exclude_rows)) terms <- setdiff(terms, exclude_rows)
  if (length(terms) == 0L) return(NULL)

  extracted <- lapply(reps, .simExtractRows,
                      table_name = table_name, terms = terms,
                      ci_lo = ci_lo, ci_hi = ci_hi)

  resi_mat <- do.call(rbind, lapply(extracted, `[[`, "resi"))
  lo_mat   <- do.call(rbind, lapply(extracted, `[[`, "lo"))
  hi_mat   <- do.call(rbind, lapply(extracted, `[[`, "hi"))

  tv       <- true_resi[terms, "RESI"]
  diff_mat <- sweep(resi_mat, 2L, tv, "-")
  in_ci    <- sweep(lo_mat,   2L, tv, "<=") & sweep(hi_mat, 2L, tv, ">=")

  # upper_coverage: P(hi >= tv) — CI upper bound is at or above the true value
  # lower_coverage: P(lo <= tv) — CI lower bound is at or below the true value
  upper_cov_mat <- sweep(hi_mat, 2L, tv, ">=")
  lower_cov_mat <- sweep(lo_mat, 2L, tv, "<=")

  data.frame(
    bias           = colMeans(diff_mat,           na.rm = TRUE),
    mse            = colMeans(diff_mat ^ 2L,      na.rm = TRUE),
    empirical_sd   = apply(resi_mat, 2L, sd,      na.rm = TRUE),
    coverage       = colMeans(in_ci,              na.rm = TRUE),
    upper_coverage = colMeans(upper_cov_mat,      na.rm = TRUE),
    lower_coverage = colMeans(lower_cov_mat,      na.rm = TRUE),
    width          = colMeans(hi_mat - lo_mat,    na.rm = TRUE),
    row.names = terms,
    stringsAsFactors = FALSE
  )
}


# ============================================================
#  insurancePlasmodeSim
# ============================================================

#' Insurance Plasmode Simulation for RESI Evaluation
#'
#' Runs a plasmode simulation study using the \code{\link{insurance}} dataset to
#' evaluate RESI confidence interval performance. In each replicate, \code{n}
#' observations are resampled with replacement from the full insurance dataset
#' (\emph{N} = 1338). The RESI point estimates from \code{\link{resi_pe}} applied
#' to the full dataset are treated as the true parameter values for computing bias,
#' MSE, CI coverage, and CI width.
#'
#' Two models are evaluated:
#' \itemize{
#'   \item \strong{lm}: \code{log10(charges) ~ ns(age, df=3) * sex + bmi + smoker + region}
#'   \item \strong{glm}: \code{I(charges > 10000) ~ ns(age, df=3) * sex + bmi + smoker + region}
#'     with \code{family = binomial()}
#' }
#' Each model is evaluated under both parametric (\code{vcovfunc = stats::vcov}) and
#' robust (\code{vcovfunc = sandwich::vcovHC}) variance settings, yielding four
#' simulation conditions.
#'
#' Parallelization is via \code{\link[parallel]{mclapply}}, which uses forking and is
#' not supported on Windows (falls back to sequential evaluation on Windows).
#'
#' @param nsim Integer, number of simulation replicates per (setting, \code{n}) cell.
#'   Default 1000. Use 10 for initial testing.
#' @param n.vec Integer vector of sample sizes. Default
#'   \code{c(50, 100, 200, 500, 1000, 2000, 5000)}.
#' @param nboot Integer, bootstrap replicates per internal \code{\link{resi}} call.
#'   Default 500. Use 10 for initial testing. Ignored when \code{ci.method != "boot"}.
#' @param alpha Numeric, CI significance level. Default 0.05.
#' @param ci.method Character, CI method passed to \code{\link{resi}}. One of
#'   \code{"boot"} (bootstrap, default), \code{"normal"} (asymptotic truncated-normal),
#'   or \code{"qf"} (asymptotic quadratic-form / Imhof). When \code{ci.method != "boot"},
#'   \code{nboot} is ignored.
#' @param output.dir Character, path to the directory where all results are saved.
#'   Created if it does not exist. Defaults to \code{"resiBootSim"},
#'   \code{"resiAsympNormalSim"}, or \code{"resiAsympQFSim"} based on
#'   \code{ci.method} when \code{NULL}.
#' @param fixed.knots Logical. If \code{TRUE}, spline knots are fixed at the
#'   empirical tertiles of \code{age} in the full insurance dataset rather than
#'   re-selected by \code{df = 3} in each bootstrap sample. Default \code{FALSE}.
#' @param mc.cores.settings Integer, cores for the outer
#'   \code{mclapply} over (setting \eqn{\times} sample size) combinations. Default 1.
#' @param mc.cores.reps Integer, cores for the inner \code{mclapply} over simulation
#'   replicates within each (setting, \code{n}) cell. Default 1.
#'
#' @return Invisibly returns the summary metrics \code{data.frame}. Side effects:
#' \itemize{
#'   \item \code{output.dir/sim_raw/<setting>_n<n>.rds}: list of per-replicate
#'     \code{anova} and \code{coefficients} tables.
#'   \item \code{output.dir/summary_table.rds}: combined metrics table with columns
#'     \code{model, vcov, n, n_success, table, term, bias, mse, coverage, width}.
#' }
#' @seealso \code{\link{simFigures}}, \code{\link{simCompareMethodsFigures}},
#'   \code{\link{resi}}, \code{\link{resi_pe}}
#' @importFrom parallel mclapply
#' @importFrom splines ns
#' @importFrom sandwich vcovHC
#' @importFrom stats lm glm vcov binomial
#' @export
insurancePlasmodeSim <- function(nsim              = 1000L,
                                  n.vec             = c(50, 100, 200, 500, 1000, 2000, 5000),
                                  nboot             = 500L,
                                  alpha             = 0.05,
                                  ci.method         = c("boot", "normal", "qf", "cf"),
                                  output.dir        = NULL,
                                  fixed.knots       = FALSE,
                                  mc.cores.settings = 1L,
                                  mc.cores.reps     = 1L) {

  ci.method <- match.arg(ci.method)
  if (is.null(output.dir)) {
    output.dir <- switch(ci.method,
      boot   = "resiBootSim",
      normal = "resiAsympNormalSim",
      qf     = "resiAsympQFSim",
      cf     = "resiAsympCFSim"
    )
  }

  insurance <- RESI::insurance

  # Pre-compute expected factor levels for plasmode resampling check
  .fvars   <- c("sex", "smoker", "region")
  .flevels <- lapply(.fvars, function(v) unique(insurance[[v]]))
  names(.flevels) <- .fvars

  ci_lo <- paste0(alpha / 2 * 100, "%")
  ci_hi <- paste0((1 - alpha / 2) * 100, "%")

  if (fixed.knots) {
    .age_knots <- quantile(insurance$age, c(1/3, 2/3))
    .age_bk    <- range(insurance$age)
    lm_formula  <- eval(bquote(log10(charges) ~ splines::ns(age, knots = .(.age_knots), Boundary.knots = .(.age_bk)) * sex + bmi + smoker + region))
    glm_formula <- eval(bquote(I(charges > 15000) ~ splines::ns(age, knots = .(.age_knots), Boundary.knots = .(.age_bk)) * sex + bmi + smoker + region))
  } else {
    lm_formula  <- log10(charges) ~ splines::ns(age, df = 3) * sex + bmi + smoker + region
    glm_formula <- I(charges > 15000) ~ splines::ns(age, df = 3) * sex + bmi + smoker + region
  }

  model_settings <- list(
    list(type = "lm",  label = "lm_parametric",
         formula = lm_formula,  family = NULL,
         vcov_name = "parametric", vcovfunc = stats::vcov),
    list(type = "lm",  label = "lm_robust",
         formula = lm_formula,  family = NULL,
         vcov_name = "robust",     vcovfunc = sandwich::vcovHC),
    list(type = "glm", label = "glm_parametric",
         formula = glm_formula, family = stats::binomial(),
         vcov_name = "parametric", vcovfunc = stats::vcov),
    list(type = "glm", label = "glm_robust",
         formula = glm_formula, family = stats::binomial(),
         vcov_name = "robust",     vcovfunc = sandwich::vcovHC)
  )

  # --- True RESI values from the full dataset ------------------------------------
  message("Computing true RESI values from full insurance dataset (n = ",
          nrow(insurance), ")...")

  true_vals <- lapply(model_settings, function(s) {
    .formula <- s$formula
    .family  <- s$family
    full_mod <- if (s$type == "lm") {
      m <- lm(.formula, data = insurance)
      m$call[["formula"]] <- .formula
      m
    } else {
      m <- glm(.formula, data = insurance, family = .family)
      m$call[["formula"]] <- .formula
      m$call[["family"]]  <- .family
      m
    }
    # Robust true values: use HC0 (no hat-value correction) so that the
    # small-sample inflation from HC3 hat weights does not enter S_true.
    # Parametric lm true values: pass rdf = n-p so that S_true = t/sqrt(rdf),
    # targeting sigma^2 = mean(resid^2) (population value, no df correction).
    true_vcovfunc <- if (s$vcov_name == "robust") {
      function(x) sandwich::vcovHC(x, type = "HC0")
    } else {
      s$vcovfunc
    }
    rdf_true <- if (s$type == "lm" && s$vcov_name == "parametric") {
      full_mod$df.residual
    } else {
      NULL
    }
    .simDirectRESI(resi_pe(full_mod, data = insurance, vcovfunc = true_vcovfunc),
                   n = nrow(insurance), rdf = rdf_true)
  })
  names(true_vals) <- sapply(model_settings, `[[`, "label")

  # --- Output directories --------------------------------------------------------
  sim_raw_dir <- file.path(output.dir, "sim_raw")
  dir.create(sim_raw_dir,                      recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output.dir, "figures"), recursive = TRUE, showWarnings = FALSE)

  # --- (setting x n) grid -------------------------------------------------------
  grid <- expand.grid(
    setting_idx = seq_along(model_settings),
    n           = n.vec,
    stringsAsFactors = FALSE
  )

  message("Running ", nsim, " replicates for each of ", nrow(grid),
          " (setting x sample size) combinations ...")

  # --- Outer mclapply over (setting x n) ----------------------------------------
  all_metrics <- parallel::mclapply(seq_len(nrow(grid)), function(g) {

    s_idx <- grid$setting_idx[g]
    n     <- grid$n[g]
    s     <- model_settings[[s_idx]]
    tv    <- true_vals[[s_idx]]

    setting_label <- paste0(s$label, "_n", n)

    # --- Inner mclapply over replicates -----------------------------------------
    reps_raw <- parallel::mclapply(seq_len(nsim), function(i) {
      # Resample with replacement until the plasmode sample is well-conditioned:
      # (a) all factor levels present; (b) for GLM, both binary outcome values
      #     must appear within each smoker stratum (prevents complete separation).
      repeat {
        sample_data <- insurance[sample(nrow(insurance), n, replace = TRUE), ]
        ok <- all(vapply(.fvars, function(v)
          all(.flevels[[v]] %in% sample_data[[v]]), logical(1L)))
        if (ok && s$type == "glm") {
          biny <- as.integer(sample_data$charges > 15000)
          ok   <- length(unique(biny)) > 1L &&
            all(tapply(biny, sample_data$smoker,
                       function(x) length(unique(x)) > 1L))
        }
        if (ok) break
      }

      # Embed the evaluated formula (and family for glm) directly into mod$call
      # so that resi()'s internal update() calls can resolve them without needing
      # 's' to be in scope (it is not accessible inside boot::boot / resi_stat).
      .formula <- s$formula
      .family  <- s$family
      mod <- tryCatch({
        m <- if (s$type == "lm") {
          lm(.formula, data = sample_data)
        } else {
          glm(.formula, data = sample_data, family = .family)
        }
        m$call[["formula"]] <- .formula
        if (s$type == "glm") m$call[["family"]] <- .family
        m
      }, error = function(e) NULL)
      if (is.null(mod)) return(NULL)

      tryCatch(
        resi(mod, data = sample_data, nboot = nboot, vcovfunc = s$vcovfunc,
             alpha = alpha, store.boot = FALSE, ci.method = ci.method),
        error = function(e) NULL
      )
    }, mc.cores = mc.cores.reps)

    reps <- Filter(Negate(is.null), reps_raw)
    n_success <- length(reps)

    # Save minimal per-replicate output (anova + coefficients tables only)
    reps_to_save <- lapply(reps, function(r) {
      list(anova = r$anova, coefficients = r$coefficients)
    })
    # Retry loop: on macOS, SIGCHLD from mclapply child processes can interrupt
    # gzfile() with "Interrupted system call"; a brief pause + retry resolves it.
    for (.attempt in seq_len(5L)) {
      saved <- tryCatch({
        Sys.sleep(0.05 * .attempt)
        saveRDS(reps_to_save,
                file = file.path(sim_raw_dir, paste0(setting_label, ".rds")))
        TRUE
      }, error = function(e) FALSE)
      if (saved) break
      if (.attempt == 5L)
        warning("saveRDS failed after 5 attempts for ", setting_label)
    }

    if (n_success == 0L) {
      warning("All replicates failed for setting: ", setting_label)
      return(NULL)
    }

    # --- Compute metrics --------------------------------------------------------
    anova_metrics <- .simComputeMetrics(reps, "anova",
                                         tv$anova, ci_lo, ci_hi,
                                         exclude_rows = "Residuals")
    coef_metrics  <- .simComputeMetrics(reps, "coefficients",
                                         tv$coefficients, ci_lo, ci_hi,
                                         exclude_rows = "(Intercept)")

    rows <- list()
    if (!is.null(anova_metrics)) {
      rows[["anova"]] <- data.frame(
        model     = s$type,
        vcov      = s$vcov_name,
        n         = n,
        n_success = n_success,
        table     = "anova",
        term      = rownames(anova_metrics),
        anova_metrics,
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    }
    if (!is.null(coef_metrics)) {
      rows[["coef"]] <- data.frame(
        model     = s$type,
        vcov      = s$vcov_name,
        n         = n,
        n_success = n_success,
        table     = "coefficients",
        term      = rownames(coef_metrics),
        coef_metrics,
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)

  }, mc.cores = mc.cores.settings)

  # --- Combine and save summary table -------------------------------------------
  summary_table <- do.call(rbind, Filter(Negate(is.null), all_metrics))
  saveRDS(summary_table, file = file.path(output.dir, "summary_table.rds"))

  message("Done. Results saved to: ", output.dir)
  invisible(summary_table)
}


# ============================================================
#  simRecomputeSummary
# ============================================================

#' Recompute Summary Table from Raw Simulation Output
#'
#' Reads the per-replicate raw \code{.rds} files saved by
#' \code{\link{insurancePlasmodeSim}} and recomputes the summary metrics table
#' (bias, MSE, coverage, upper_coverage, lower_coverage, width).  Use this
#' whenever \code{.simComputeMetrics} has been updated (e.g. new metrics added)
#' without needing to rerun the full simulation.
#'
#' True RESI values are re-estimated from the full \code{\link{insurance}} dataset
#' using the same formulae and variance functions as the original simulation.
#'
#' @param output.dir Character, directory containing simulation output.
#'   Default \code{"resiBootSim"}.
#' @param alpha Numeric, CI significance level used in the original simulation.
#'   Default 0.05.
#' @param fixed.knots Logical. Must match the value used in the original
#'   \code{\link{insurancePlasmodeSim}} call so that the true RESI values are
#'   computed from the same model formula. Default \code{FALSE}.
#'
#' @return Invisibly returns the updated summary \code{data.frame}.
#'   Overwrites \code{output.dir/summary_table.rds}.
#' @seealso \code{\link{insurancePlasmodeSim}}, \code{\link{simFigures}}
#' @importFrom splines ns
#' @importFrom sandwich vcovHC
#' @importFrom stats lm glm vcov binomial
#' @export
simRecomputeSummary <- function(output.dir  = "resiBootSim",
                                 alpha       = 0.05,
                                 fixed.knots = FALSE) {

  insurance <- RESI::insurance

  ci_lo <- paste0(alpha / 2 * 100, "%")
  ci_hi <- paste0((1 - alpha / 2) * 100, "%")

  if (fixed.knots) {
    .age_knots <- quantile(insurance$age, c(1/3, 2/3))
    .age_bk    <- range(insurance$age)
    lm_formula  <- eval(bquote(log10(charges) ~ splines::ns(age, knots = .(.age_knots), Boundary.knots = .(.age_bk)) * sex + bmi + smoker + region))
    glm_formula <- eval(bquote(I(charges > 15000) ~ splines::ns(age, knots = .(.age_knots), Boundary.knots = .(.age_bk)) * sex + bmi + smoker + region))
  } else {
    lm_formula  <- log10(charges) ~ splines::ns(age, df = 3) * sex + bmi + smoker + region
    glm_formula <- I(charges > 15000) ~ splines::ns(age, df = 3) * sex + bmi + smoker + region
  }

  model_settings <- list(
    list(type = "lm",  label = "lm_parametric",
         formula = lm_formula,  family = NULL,
         vcov_name = "parametric", vcovfunc = stats::vcov),
    list(type = "lm",  label = "lm_robust",
         formula = lm_formula,  family = NULL,
         vcov_name = "robust",     vcovfunc = sandwich::vcovHC),
    list(type = "glm", label = "glm_parametric",
         formula = glm_formula, family = stats::binomial(),
         vcov_name = "parametric", vcovfunc = stats::vcov),
    list(type = "glm", label = "glm_robust",
         formula = glm_formula, family = stats::binomial(),
         vcov_name = "robust",     vcovfunc = sandwich::vcovHC)
  )

  message("Recomputing true RESI values from full insurance dataset...")
  true_vals <- lapply(model_settings, function(s) {
    .formula <- s$formula
    .family  <- s$family
    full_mod <- if (s$type == "lm") {
      m <- lm(.formula, data = insurance)
      m$call[["formula"]] <- .formula
      m
    } else {
      m <- glm(.formula, data = insurance, family = .family)
      m$call[["formula"]] <- .formula
      m$call[["family"]]  <- .family
      m
    }
    true_vcovfunc <- if (s$vcov_name == "robust") {
      function(x) sandwich::vcovHC(x, type = "HC0")
    } else {
      s$vcovfunc
    }
    rdf_true <- if (s$type == "lm" && s$vcov_name == "parametric") {
      full_mod$df.residual
    } else {
      NULL
    }
    .simDirectRESI(resi_pe(full_mod, data = insurance, vcovfunc = true_vcovfunc),
                   n = nrow(insurance), rdf = rdf_true)
  })
  names(true_vals) <- sapply(model_settings, `[[`, "label")

  sim_raw_dir <- file.path(output.dir, "sim_raw")
  raw_files   <- list.files(sim_raw_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(raw_files) == 0L)
    stop("No raw .rds files found in: ", sim_raw_dir)

  setting_labels <- sapply(model_settings, `[[`, "label")

  all_metrics <- lapply(raw_files, function(f) {
    fname <- sub("\\.rds$", "", basename(f))
    n     <- suppressWarnings(as.integer(sub(".*_n", "", fname)))
    label <- sub("_n[0-9]+$", "", fname)

    s_idx <- which(setting_labels == label)
    if (length(s_idx) == 0L || is.na(n)) {
      warning("Could not parse setting/n from filename: ", basename(f))
      return(NULL)
    }
    s  <- model_settings[[s_idx]]
    tv <- true_vals[[s_idx]]

    reps_raw  <- readRDS(f)
    n_success <- length(reps_raw)
    if (n_success == 0L) return(NULL)

    anova_metrics <- .simComputeMetrics(reps_raw, "anova",
                                         tv$anova, ci_lo, ci_hi,
                                         exclude_rows = "Residuals")
    coef_metrics  <- .simComputeMetrics(reps_raw, "coefficients",
                                         tv$coefficients, ci_lo, ci_hi,
                                         exclude_rows = "(Intercept)")

    rows <- list()
    if (!is.null(anova_metrics)) {
      rows[["anova"]] <- data.frame(
        model     = s$type,
        vcov      = s$vcov_name,
        n         = n,
        n_success = n_success,
        table     = "anova",
        term      = rownames(anova_metrics),
        anova_metrics,
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    }
    if (!is.null(coef_metrics)) {
      rows[["coef"]] <- data.frame(
        model     = s$type,
        vcov      = s$vcov_name,
        n         = n,
        n_success = n_success,
        table     = "coefficients",
        term      = rownames(coef_metrics),
        coef_metrics,
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  })

  summary_table <- do.call(rbind, Filter(Negate(is.null), all_metrics))
  saveRDS(summary_table, file = file.path(output.dir, "summary_table.rds"))
  message("Done. Updated summary_table.rds saved to: ", output.dir)
  invisible(summary_table)
}


# ============================================================
#  simEstimatorFigures
# ============================================================

#' Estimator Comparison Figures: t2S/f2S vs z2S/chisq2S
#'
#' Reads the per-replicate raw \code{.rds} files from one simulation directory and
#' produces Bias + MSE comparison figures contrasting the RESI estimators actually
#' used for \code{lm} models against the simpler z-statistic / chi-squared alternatives:
#' \itemize{
#'   \item \strong{Coefficients table}: \code{\link{t2S}} (used) vs \code{\link{z2S}} (alternative)
#'   \item \strong{Anova table}: \code{\link{f2S}} (used) vs \code{\link{chisq2S}} (alternative)
#' }
#' Only \code{lm} models are included; GLM models are skipped because z2S and
#' chisq2S are the natural estimators there.
#'
#' @param sim.dir Character, directory containing simulation output with a
#'   \code{sim_raw/} sub-directory.  Default \code{"resiBootSim"}.
#' @param figures.dir Character, output directory for PDFs.
#'   Default \code{file.path(sim.dir, "figures", "estimator_compare")}.
#' @param alpha Numeric, nominal level (used only in figure titles). Default 0.05.
#' @param fixed.knots Logical. Must match the value used in the original
#'   \code{\link{insurancePlasmodeSim}} call. Default \code{FALSE}.
#'
#' @return Invisibly returns the combined estimator-comparison \code{data.frame}.
#' @seealso \code{\link{simCompareMethodsFigures}}, \code{\link{simRecomputeSummary}}
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot lines points abline legend par axis plot.new layout
#' @importFrom stats lm vcov
#' @importFrom sandwich vcovHC
#' @importFrom splines ns
#' @export
simEstimatorFigures <- function(
    sim.dir     = "resiBootSim",
    figures.dir = NULL,
    alpha       = 0.05,
    fixed.knots = FALSE
) {
  if (is.null(figures.dir))
    figures.dir <- file.path(sim.dir, "figures", "estimator_compare")
  dir.create(figures.dir, recursive = TRUE, showWarnings = FALSE)

  insurance <- RESI::insurance

  # ---- Formula ---------------------------------------------------------------
  if (fixed.knots) {
    .age_knots <- quantile(insurance$age, c(1/3, 2/3))
    .age_bk    <- range(insurance$age)
    lm_formula <- eval(bquote(
      log10(charges) ~ splines::ns(age, knots = .(.age_knots),
                                   Boundary.knots = .(.age_bk)) * sex + bmi + smoker + region))
  } else {
    lm_formula <- log10(charges) ~ splines::ns(age, df = 3) * sex + bmi + smoker + region
  }

  # ---- lm settings only ------------------------------------------------------
  model_settings <- list(
    list(label = "lm_parametric", vcov_name = "parametric",
         vcovfunc = stats::vcov),
    list(label = "lm_robust",     vcov_name = "robust",
         vcovfunc = sandwich::vcovHC)
  )

  # ---- True S from full insurance dataset ------------------------------------
  message("Computing true RESI values from full insurance dataset...")
  true_vals <- lapply(model_settings, function(s) {
    m <- lm(lm_formula, data = insurance)
    m$call[["formula"]] <- lm_formula
    tv_vcov <- if (s$vcov_name == "robust") {
      function(x) sandwich::vcovHC(x, type = "HC0")
    } else {
      s$vcovfunc
    }
    rdf_tv <- if (s$vcov_name == "parametric") m$df.residual else NULL
    .simDirectRESI(resi_pe(m, data = insurance, vcovfunc = tv_vcov),
                   n = nrow(insurance), rdf = rdf_tv)
  })
  names(true_vals) <- sapply(model_settings, `[[`, "label")

  # ---- Colours / layout constants --------------------------------------------
  est_cols <- c("t2S / f2S" = "#1B6CA8", "z2S / chisq2S" = "#CC5500")
  row_h    <- 3.2
  leg_h    <- 0.55
  raw_dir  <- file.path(sim.dir, "sim_raw")

  # ---- Helper: extract one column from anova/coef table per replicate --------
  .extract_col <- function(reps, tbl_name, terms, col) {
    m <- matrix(NA_real_, nrow = length(reps), ncol = length(terms),
                dimnames = list(NULL, terms))
    for (ri in seq_along(reps)) {
      tab <- reps[[ri]][[tbl_name]]
      idx <- match(terms, rownames(tab))
      ok  <- !is.na(idx)
      if (any(ok)) m[ri, ok] <- tab[idx[ok], col]
    }
    m
  }

  # ---- Loop over lm settings -------------------------------------------------
  all_combined <- list()

  for (s in model_settings) {
    lbl   <- s$label
    vtype <- s$vcov_name
    tv    <- true_vals[[lbl]]

    raw_files <- list.files(raw_dir,
                            pattern = paste0("^", lbl, "_n[0-9]+\\.rds$"),
                            full.names = TRUE)
    if (length(raw_files) == 0L) {
      message("No raw files found for: ", lbl, " — skipping")
      next
    }

    # -- Per-n summaries -------------------------------------------------------
    per_n <- lapply(raw_files, function(f) {
      n_sim <- as.integer(sub(".*_n([0-9]+)\\.rds$", "\\1", basename(f)))
      reps  <- readRDS(f)
      if (length(reps) == 0L) return(NULL)

      rows <- list()

      # Anova: f2S (RESI col) vs chisq2S(F*Df, Df, n)
      a_terms <- setdiff(rownames(tv$anova), "Residuals")
      if (length(a_terms) > 0L) {
        resi_f <- .extract_col(reps, "anova", a_terms, "RESI")
        resi_c <- matrix(NA_real_, nrow = length(reps), ncol = length(a_terms),
                         dimnames = list(NULL, a_terms))
        for (ri in seq_along(reps)) {
          tab <- reps[[ri]]$anova
          idx <- match(a_terms, rownames(tab))
          ok  <- !is.na(idx)
          if (any(ok)) {
            F_i  <- tab[idx[ok], "F"]
            Df_i <- tab[idx[ok], "Df"]
            resi_c[ri, ok] <- chisq2S(F_i * Df_i, Df_i, n_sim)
          }
        }
        tv_a <- tv$anova[a_terms, "RESI"]
        diff_f <- sweep(resi_f, 2L, tv_a, "-")
        diff_c <- sweep(resi_c, 2L, tv_a, "-")
        rows$anova <- rbind(
          data.frame(model = "lm", vcov = vtype, n = n_sim, table = "anova",
                     term = a_terms, estimator = "t2S / f2S",
                     bias = colMeans(diff_f, na.rm = TRUE),
                     mse  = colMeans(diff_f^2, na.rm = TRUE),
                     row.names = NULL, stringsAsFactors = FALSE),
          data.frame(model = "lm", vcov = vtype, n = n_sim, table = "anova",
                     term = a_terms, estimator = "z2S / chisq2S",
                     bias = colMeans(diff_c, na.rm = TRUE),
                     mse  = colMeans(diff_c^2, na.rm = TRUE),
                     row.names = NULL, stringsAsFactors = FALSE)
        )
      }

      # Coefficients: t2S (RESI col) vs z2S(t, n)
      c_terms <- setdiff(rownames(tv$coefficients), "(Intercept)")
      if (length(c_terms) > 0L) {
        coef_col <- "t value"
        resi_t <- .extract_col(reps, "coefficients", c_terms, "RESI")
        t_mat  <- .extract_col(reps, "coefficients", c_terms, coef_col)
        resi_z <- matrix(z2S(t_mat, n_sim), nrow = nrow(t_mat), ncol = ncol(t_mat),
                         dimnames = dimnames(t_mat))
        tv_c   <- tv$coefficients[c_terms, "RESI"]
        diff_t <- sweep(resi_t, 2L, tv_c, "-")
        diff_z <- sweep(resi_z, 2L, tv_c, "-")
        rows$coef <- rbind(
          data.frame(model = "lm", vcov = vtype, n = n_sim, table = "coefficients",
                     term = c_terms, estimator = "t2S / f2S",
                     bias = colMeans(diff_t, na.rm = TRUE),
                     mse  = colMeans(diff_t^2, na.rm = TRUE),
                     row.names = NULL, stringsAsFactors = FALSE),
          data.frame(model = "lm", vcov = vtype, n = n_sim, table = "coefficients",
                     term = c_terms, estimator = "z2S / chisq2S",
                     bias = colMeans(diff_z, na.rm = TRUE),
                     mse  = colMeans(diff_z^2, na.rm = TRUE),
                     row.names = NULL, stringsAsFactors = FALSE)
        )
      }

      do.call(rbind, rows)
    })

    combined <- do.call(rbind, Filter(Negate(is.null), per_n))
    if (is.null(combined) || nrow(combined) == 0L) next
    all_combined[[lbl]] <- combined

    n_vals <- sort(unique(combined$n))

    # -- Figures ---------------------------------------------------------------
    for (tbl_name in c("anova", "coefficients")) {
      sub <- combined[combined$table == tbl_name, ]
      if (nrow(sub) == 0L) next

      terms   <- unique(sub$term)
      n_terms <- length(terms)
      tv_tbl  <- if (tbl_name == "anova") tv$anova else tv$coefficients

      fig_path <- file.path(figures.dir,
                            paste0("estimator_compare_", lbl, "_", tbl_name, ".pdf"))
      grDevices::pdf(fig_path, width = 7, height = row_h * n_terms + leg_h)

      # Layout: n_terms rows × 2 cols + 1 legend row
      n_panels <- 2L * n_terms
      lay_mat  <- matrix(c(seq_len(n_panels),
                           rep(n_panels + 1L, 2L)),
                         nrow = n_terms + 1L, ncol = 2L, byrow = TRUE)
      graphics::layout(lay_mat,
                       widths  = c(3.2, 3.2),
                       heights = c(rep(row_h, n_terms), leg_h))

      for (ti in seq_along(terms)) {
        term      <- terms[[ti]]
        tdata     <- sub[sub$term == term, ]
        true_s    <- if (term %in% rownames(tv_tbl)) tv_tbl[term, "RESI"] else NA_real_
        term_title <- if (is.finite(true_s)) sprintf("%s: S=%.3f", term, true_s) else term

        for (metric in c("bias", "mse")) {
          ylab <- if (metric == "bias") "Bias" else "MSE"
          ylim <- range(tdata[[metric]], na.rm = TRUE)
          if (metric == "bias") ylim <- range(c(ylim, 0), na.rm = TRUE)
          ylim <- ylim + diff(ylim) * c(-0.05, 0.05)

          graphics::par(mar = c(3.2, 3.0, 2.2, 0.5), mgp = c(1.8, 0.5, 0))
          graphics::plot(NULL, xlim = range(n_vals), ylim = ylim,
                         xlab = "Sample Size", ylab = ylab,
                         main = if (metric == "bias") term_title else "",
                         log = "x", xaxt = "n", bty = "l")
          graphics::axis(1, at = n_vals, labels = n_vals, las = 2L,
                         mgp = c(1.7, 0.35, 0))
          if (metric == "bias") graphics::abline(h = 0, lty = 2L, col = "gray40")

          for (est in names(est_cols)) {
            ed <- tdata[tdata$estimator == est, ]
            ed <- ed[order(ed$n), ]
            if (nrow(ed) == 0L) next
            graphics::lines(ed$n, ed[[metric]], col = est_cols[[est]], lwd = 2L)
            graphics::points(ed$n, ed[[metric]], col = est_cols[[est]], pch = 16L, cex = 0.9)
          }
        }
      }

      # Bottom legend
      est_labels <- if (tbl_name == "anova") c("f2S (used)", "chisq2S") else c("t2S (used)", "z2S")
      graphics::par(mar = c(0.2, 0.2, 0.2, 0.2))
      graphics::plot.new()
      graphics::legend("center",
                       legend = paste0(est_labels, "    "),
                       col    = unname(est_cols),
                       pch    = 16L, lwd = 2L, lty = 1L,
                       bty    = "n", cex = 0.9, ncol = 2L,
                       title  = "Estimator", title.font = 2L)
      grDevices::dev.off()
      message("Saved: ", fig_path)
    }
  }

  message("Estimator comparison figures saved to: ", figures.dir)
  invisible(do.call(rbind, all_combined))
}


# ============================================================
#  simFigures
# ============================================================

#' Simulation Performance Figures for RESI Evaluation
#'
#' Produces performance figures from the output of \code{\link{insurancePlasmodeSim}}.
#' Creates one PDF figure per (model type) \eqn{\times} (variance estimator) combination
#' (4 figures total by default). Each figure is a 2 \eqn{\times} 4 panel grid:
#' \itemize{
#'   \item \strong{Top row}: Anova-table RESI metrics (Bias, MSE, CI Coverage, CI Width)
#'   \item \strong{Bottom row}: Coefficients-table RESI metrics (same metrics; intercept excluded)
#'   \item \strong{Colors}: one colored line per model term (high-contrast matte palette)
#'   \item \strong{x-axis}: sample size on a log scale
#' }
#' Dashed reference lines are drawn at zero for Bias and at \eqn{1 - \alpha} for CI Coverage.
#'
#' @param output.dir Character, directory containing simulation output from
#'   \code{\link{insurancePlasmodeSim}}. Default \code{"resiBootSim"}.
#' @param alpha Numeric, nominal CI level used for the coverage reference line.
#'   Default 0.05.
#' @param ci.label Character, label for the CI method used in figure titles and file
#'   names. Default \code{NULL}, which auto-detects from \code{output.dir}:
#'   \code{"boot"} if the directory name contains "Boot"/"boot",
#'   \code{"normal"} if it contains "Normal"/"normal",
#'   \code{"qf"} if it contains "QF"/"qf", otherwise \code{basename(output.dir)}.
#'
#' @return Invisibly returns the summary metrics \code{data.frame}. Saves PDF figures
#'   to \code{file.path(output.dir, "figures")}, named
#'   \code{sim_<model>_<vcov>_<ci.label>.pdf}.
#' @seealso \code{\link{insurancePlasmodeSim}}
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot lines points abline legend par axis plot.new layout
#' @export
simFigures <- function(output.dir = "resiBootSim",
                        alpha      = 0.05,
                        ci.label   = NULL) {

  # Auto-detect ci.label from output.dir when not supplied
  if (is.null(ci.label)) {
    ci.label <- if (grepl("[Bb]oot",   output.dir)) "boot" else
                if (grepl("[Nn]ormal", output.dir)) "normal" else
                if (grepl("[Qq][Ff]",  output.dir)) "qf" else
                basename(output.dir)
  }

  # High-contrast matte palette (matplotlib tab10 + extensions)
  .sim_pal <- c(
    "#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD",
    "#8C564B", "#E377C2", "#17BECF", "#BCBD22", "#7F7F7F",
    "#AEC7E8", "#98DF8A"
  )

  summary_table <- readRDS(file.path(output.dir, "summary_table.rds"))
  if (is.null(summary_table) || nrow(summary_table) == 0L)
    stop("summary_table.rds in '", output.dir, "' is empty or NULL. ",
         "Re-run insurancePlasmodeSim() (all simulation replicates may have ",
         "failed -- check for 'All replicates failed' warnings).")

  figures_dir   <- file.path(output.dir, "figures")
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  n_vals        <- sort(unique(summary_table$n))
  metrics       <- c("bias", "mse", "coverage", "width")
  metric_labels <- c("Bias", "MSE", "CI Coverage", "CI Width")

  for (mtype in c("lm", "glm")) {
    for (vtype in c("parametric", "robust")) {

      sub <- summary_table[summary_table$model == mtype &
                             summary_table$vcov  == vtype, ]
      if (is.null(sub) || nrow(sub) == 0L) next

      fig_path <- file.path(
        figures_dir,
        paste0("sim_", mtype, "_", vtype, "_", ci.label, ".pdf")
      )

      anova_terms <- unique(sub$term[sub$table == "anova"])
      coef_terms  <- unique(sub$term[sub$table == "coefficients"])

      anova_cols <- setNames(
        .sim_pal[seq_along(anova_terms)], anova_terms)
      coef_cols  <- setNames(
        .sim_pal[seq_along(coef_terms)],  coef_terms)

      row_info <- list(
        list(table       = "anova",
             terms       = anova_terms,
             cols        = anova_cols,
             main_prefix = "Anova"),
        list(table       = "coefficients",
             terms       = coef_terms,
             cols        = coef_cols,
             main_prefix = "Coef")
      )

      # Layout: 4 rows × 5 cols; coverage column (col 3) is split into two
      # half-height panels (upper coverage on top, lower coverage on bottom).
      # Row 1:  1   2   3   4   5   (anova bias, mse, upper cov, width, legend)
      # Row 2:  1   2   6   4   5   (spans continue; anova lower cov)
      # Row 3:  7   8   9  10  11   (coef bias, mse, upper cov, width, legend)
      # Row 4:  7   8  12  10  11   (spans continue; coef lower cov)
      grDevices::pdf(fig_path, width = 15, height = 7)
      graphics::layout(
        matrix(c( 1,  2,  3,  4,  5,
                  1,  2,  6,  4,  5,
                  7,  8,  9, 10, 11,
                  7,  8, 12, 10, 11), nrow = 4L, byrow = TRUE),
        widths  = c(rep(3, 4), 2.2),
        heights = c(1, 1, 1, 1)
      )

      for (ri in row_info) {
        sub_tbl <- sub[sub$table == ri$table, ]
        cex_leg <- min(1.0, 9 / length(ri$terms))

        # Helper: draw lines/points for all terms for a given metric
        draw_lines <- function(metric) {
          for (term in ri$terms) {
            td <- sub_tbl[sub_tbl$term == term, ]
            td <- td[order(td$n), ]
            graphics::lines(td$n,  td[[metric]], col = ri$cols[term], lwd = 2L)
            graphics::points(td$n, td[[metric]], col = ri$cols[term], pch = 16L,
                             cex = 0.8)
          }
        }

        # Shared y-range for both coverage sub-panels
        cov_ylim <- range(
          c(sub_tbl[["upper_coverage"]], sub_tbl[["lower_coverage"]],
            1 - alpha/2 - 0.02, 1.01),
          na.rm = TRUE
        )

        # --- Panel: Bias ---
        bias_rows <- if (mtype == "glm") sub_tbl$n >= 500 else rep(TRUE, nrow(sub_tbl))
        ylim <- range(c(sub_tbl[bias_rows, "bias"], 0), na.rm = TRUE)
        graphics::par(mar = c(2.8, 2.8, 1.8, 0.4), mgp = c(1.7, 0.45, 0))
        graphics::plot(NULL, xlim = range(n_vals), ylim = ylim,
                       xlab = "Sample Size", ylab = "Bias",
                       main = paste(toupper(mtype), ri$main_prefix, "Bias",
                                    paste0("(", vtype, ")")),
                       log = "x", xaxt = "n", bty = "l")
        graphics::axis(1, at = n_vals, labels = n_vals, las = 2L,
                       mgp = c(1.7, 0.35, 0))
        graphics::abline(h = 0, lty = 2L, col = "gray40")
        draw_lines("bias")

        # --- Panel: MSE ---
        mse_rows <- if (mtype == "glm") sub_tbl$n >= 500 else rep(TRUE, nrow(sub_tbl))
        ylim <- if (any(is.finite(sub_tbl[mse_rows, "mse"]))) range(sub_tbl[mse_rows, "mse"], na.rm = TRUE) else c(0, 1)
        graphics::par(mar = c(2.8, 2.8, 1.8, 0.4), mgp = c(1.7, 0.45, 0))
        graphics::plot(NULL, xlim = range(n_vals), ylim = ylim,
                       xlab = "Sample Size", ylab = "MSE",
                       main = paste(toupper(mtype), ri$main_prefix, "MSE",
                                    paste0("(", vtype, ")")),
                       log = "x", xaxt = "n", bty = "l")
        graphics::axis(1, at = n_vals, labels = n_vals, las = 2L,
                       mgp = c(1.7, 0.35, 0))
        draw_lines("mse")

        # --- Panel: Upper Coverage (top half-panel, no x-axis) ---
        graphics::par(mar = c(0.3, 2.8, 1.8, 0.4), mgp = c(1.7, 0.45, 0))
        graphics::plot(NULL, xlim = range(n_vals), ylim = cov_ylim,
                       xlab = "", ylab = "Upper Cov.",
                       main = paste(toupper(mtype), ri$main_prefix, "Coverage",
                                    paste0("(", vtype, ")")),
                       log = "x", xaxt = "n", bty = "l")
        graphics::abline(h = 1 - alpha/2, lty = 2L, col = "gray40")
        draw_lines("upper_coverage")

        # --- Panel: Width ---
        width_rows <- if (mtype == "glm") sub_tbl$n >= 500 else rep(TRUE, nrow(sub_tbl))
        ylim <- if (any(is.finite(sub_tbl[width_rows, "width"]))) range(sub_tbl[width_rows, "width"], na.rm = TRUE) else c(0, 1)
        graphics::par(mar = c(2.8, 2.8, 1.8, 0.4), mgp = c(1.7, 0.45, 0))
        graphics::plot(NULL, xlim = range(n_vals), ylim = ylim,
                       xlab = "Sample Size", ylab = "CI Width",
                       main = paste(toupper(mtype), ri$main_prefix, "CI Width",
                                    paste0("(", vtype, ")")),
                       log = "x", xaxt = "n", bty = "l")
        graphics::axis(1, at = n_vals, labels = n_vals, las = 2L,
                       mgp = c(1.7, 0.35, 0))
        draw_lines("width")

        # --- Legend panel (spans both coverage half-rows) ---
        graphics::par(mar = c(0.5, 0.3, 0.5, 0.3))
        graphics::plot.new()
        graphics::legend(
          "center",
          legend = ri$terms,
          col    = ri$cols[ri$terms],
          lwd    = 2L,
          pch    = 16L,
          bty    = "n",
          cex    = cex_leg,
          title  = paste(ri$main_prefix, "terms"),
          title.font = 2L
        )

        # --- Panel: Lower Coverage (bottom half-panel, with x-axis) ---
        graphics::par(mar = c(2.8, 2.8, 0.3, 0.4), mgp = c(1.7, 0.45, 0))
        graphics::plot(NULL, xlim = range(n_vals), ylim = cov_ylim,
                       xlab = "Sample Size", ylab = "Lower Cov.",
                       main = "",
                       log = "x", xaxt = "n", bty = "l")
        graphics::axis(1, at = n_vals, labels = n_vals, las = 2L,
                       mgp = c(1.7, 0.35, 0))
        graphics::abline(h = 1 - alpha/2, lty = 2L, col = "gray40")
        draw_lines("lower_coverage")
      }

      grDevices::dev.off()
      message("Saved: ", fig_path)
    }
  }

  invisible(summary_table)
}


# ============================================================
#  simCompareMethodsFigures
# ============================================================

#' Per-Term CI Method Comparison Figures
#'
#' Reads simulation summary tables from multiple output directories (one per CI
#' method) and produces one PDF figure per (model type \eqn{\times} variance
#' estimator \eqn{\times} table \eqn{\times} term) combination. Each figure
#' shows the standard performance metrics (Bias, MSE, Upper Coverage, Lower
#' Coverage, CI Width) across sample sizes, with one colored line per CI method.
#'
#' @param output.dirs Named character vector mapping CI method labels to their
#'   simulation output directories. Default:
#'   \code{c(boot = "resiBootSim", normal = "resiAsympNormalSim",
#'   qf = "resiAsympQFSim", cf = "resiAsympCFSim")}.
#'   Directories that do not exist are silently skipped.
#' @param figures.dir Character, directory where comparison figures are saved.
#'   Default: \code{file.path(output.dirs[1], "figures", "method_comparison")}.
#' @param alpha Numeric, nominal CI level used for coverage reference lines.
#'   Default 0.05.
#' @param fixed.knots Logical. Must match the value used in the original
#'   \code{\link{insurancePlasmodeSim}} call; controls how the true RESI values
#'   are re-computed from the full dataset for the coverage-quantile figures.
#'   Default \code{FALSE}.
#'
#' @return Invisibly returns the combined summary \code{data.frame}. Saves to
#'   \code{figures.dir}:
#' \itemize{
#'   \item \code{compare_<model>_<vcov>_<table>.pdf} — per-metric lines across
#'     sample sizes, all terms in one PDF (rows = terms).
#'   \item \code{covquant_<model>_<vcov>_<table>.pdf} — coverage-quantile density
#'     curves, \eqn{(S_{\rm true} - \rm LCI)/(\rm UCI - \rm LCI)}, rows = terms,
#'     columns = sample sizes.
#' }
#' @seealso \code{\link{insurancePlasmodeSim}}, \code{\link{simFigures}}
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot lines points abline legend par axis plot.new mtext
#' @importFrom splines ns
#' @importFrom sandwich vcovHC
#' @importFrom stats lm glm vcov binomial density
#' @export
simCompareMethodsFigures <- function(
    output.dirs = c(boot   = "resiBootSim",
                    normal = "resiAsympNormalSim",
                    qf     = "resiAsympQFSim",
                    cf     = "resiAsympCFSim"),
    figures.dir = NULL,
    alpha       = 0.05,
    fixed.knots = FALSE) {

  method_names <- names(output.dirs)
  if (is.null(method_names) || any(method_names == ""))
    stop("output.dirs must be a named vector (names = CI method labels)")

  if (is.null(figures.dir))
    figures.dir <- file.path(output.dirs[[1]], "figures", "method_comparison")
  dir.create(figures.dir, recursive = TRUE, showWarnings = FALSE)

  # Read and merge tables, skipping missing directories
  tables <- lapply(seq_along(output.dirs), function(i) {
    f <- file.path(output.dirs[[i]], "summary_table.rds")
    if (!file.exists(f)) {
      warning("summary_table.rds not found in: ", output.dirs[[i]], " -- skipping")
      return(NULL)
    }
    tbl <- readRDS(f)
    tbl$ci_method <- method_names[[i]]
    tbl
  })
  tables <- Filter(Negate(is.null), tables)
  if (length(tables) == 0L) stop("No valid summary tables found")
  combined <- do.call(rbind, tables)

  # Keep only methods that were successfully loaded, in supplied order
  avail_methods <- intersect(method_names, unique(combined$ci_method))

  # High-contrast palette and shapes per method
  base_cols   <- c("#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD")
  method_cols <- setNames(base_cols[seq_along(avail_methods)], avail_methods)
  n_vals <- sort(unique(combined$n))
  z_ref  <- qnorm(1 - alpha / 2)

  # Compute true RESI from full insurance dataset for per-term title labels
  true_S_lookup <- tryCatch({
    insurance_full <- RESI::insurance
    if (fixed.knots) {
      .age_knots_ts <- stats::quantile(insurance_full$age, c(1/3, 2/3))
      .age_bk_ts    <- range(insurance_full$age)
      lm_form_ts    <- eval(bquote(
        log10(charges) ~ splines::ns(age, knots = .(.age_knots_ts),
                                     Boundary.knots = .(.age_bk_ts)) *
          sex + bmi + smoker + region))
      glm_form_ts   <- eval(bquote(
        I(charges > 15000) ~ splines::ns(age, knots = .(.age_knots_ts),
                                         Boundary.knots = .(.age_bk_ts)) *
          sex + bmi + smoker + region))
    } else {
      lm_form_ts  <- log10(charges) ~ splines::ns(age, df = 3) *
        sex + bmi + smoker + region
      glm_form_ts <- I(charges > 15000) ~ splines::ns(age, df = 3) *
        sex + bmi + smoker + region
    }
    ts_settings <- list(
      list(mtype="lm",  vtype="parametric", form=lm_form_ts,  fam=NULL,
           robust=FALSE),
      list(mtype="lm",  vtype="robust",     form=lm_form_ts,  fam=NULL,
           robust=TRUE),
      list(mtype="glm", vtype="parametric", form=glm_form_ts,
           fam=stats::binomial(), robust=FALSE),
      list(mtype="glm", vtype="robust",     form=glm_form_ts,
           fam=stats::binomial(), robust=TRUE)
    )
    result_ts <- list()
    for (s_ts in ts_settings) {
      key_ts <- paste(s_ts$mtype, s_ts$vtype, sep = "_")
      m_ts <- tryCatch({
        if (s_ts$mtype == "lm") {
          m <- stats::lm(s_ts$form, data = insurance_full)
          m$call[["formula"]] <- s_ts$form; m
        } else {
          m <- stats::glm(s_ts$form, data = insurance_full, family = s_ts$fam)
          m$call[["formula"]] <- s_ts$form
          m$call[["family"]]  <- s_ts$fam; m
        }
      }, error = function(e) NULL)
      if (is.null(m_ts)) next
      tv_ts <- if (s_ts$robust)
        function(x) sandwich::vcovHC(x, type = "HC0") else stats::vcov
      pe_ts <- tryCatch(
        resi_pe(m_ts, vcovfunc = tv_ts, unbiased = TRUE),
        error = function(e) NULL)
      if (!is.null(pe_ts)) result_ts[[key_ts]] <- pe_ts
    }
    result_ts
  }, error = function(e) {
    warning("Could not compute true RESI for titles: ", conditionMessage(e))
    list()
  })

  get_true_S <- function(mtype, vtype, tbl_name, term) {
    key_ts <- paste(mtype, vtype, sep = "_")
    pe_ts  <- true_S_lookup[[key_ts]]
    if (is.null(pe_ts)) return(NA_real_)
    tbl_ts <- pe_ts[[tbl_name]]
    if (is.null(tbl_ts) || !(term %in% rownames(tbl_ts))) return(NA_real_)
    tbl_ts[term, "RESI"]
  }

  row_h <- 3.2
  leg_h <- 0.55

  for (mtype in c("lm", "glm")) {
    for (vtype in c("parametric", "robust")) {
      for (tbl_name in c("anova", "coefficients")) {

        sub <- combined[
          combined$model == mtype &
          combined$vcov  == vtype &
          combined$table == tbl_name, ]
        if (nrow(sub) == 0L) next

        terms   <- unique(sub$term)
        n_terms <- length(terms)

        # Helper: draw method lines/points for a given metric column
        draw_method_lines <- function(term_data, metric) {
          for (meth in avail_methods) {
            md <- term_data[term_data$ci_method == meth, ]
            md <- md[order(md$n), ]
            if (nrow(md) == 0L) next
            graphics::lines(md$n,  md[[metric]],
                            col = method_cols[[meth]], lwd = 2L)
            graphics::points(md$n, md[[metric]],
                             col = method_cols[[meth]], pch = 16L, cex = 0.8)
          }
        }

        # ======================================================
        # PDF 1: Estimator figures (Bias + MSE only)
        # One row per term, legend at bottom spanning both columns
        # ======================================================
        fig_path_est <- file.path(
          figures.dir,
          paste0("estimator_", mtype, "_", vtype, "_", tbl_name, ".pdf")
        )
        lay_mat_est <- matrix(0L, nrow = n_terms + 1L, ncol = 2L)
        for (ti in seq_len(n_terms)) {
          lay_mat_est[ti, 1L] <- 2L * (ti - 1L) + 1L   # Bias
          lay_mat_est[ti, 2L] <- 2L * (ti - 1L) + 2L   # MSE
        }
        leg_id_est <- 2L * n_terms + 1L
        lay_mat_est[n_terms + 1L, ] <- leg_id_est

        grDevices::pdf(fig_path_est, width = 7,
                       height = row_h * n_terms + leg_h)
        graphics::layout(lay_mat_est,
                         widths  = c(3.2, 3.2),
                         heights = c(rep(row_h, n_terms), leg_h))

        for (ti in seq_along(terms)) {
          term      <- terms[[ti]]
          term_data <- sub[sub$term == term, ]
          true_s    <- get_true_S(mtype, vtype, tbl_name, term)
          term_title <- if (is.finite(true_s))
            sprintf("%s: S=%.3f", term, true_s) else term

          est_rows <- if (mtype == "glm") term_data$n >= 500 else
            rep(TRUE, nrow(term_data))

          # ---- Bias ----
          ylim <- range(c(term_data[est_rows, "bias"], -0.02, 0.01), na.rm = TRUE)
          graphics::par(mar = c(3.2, 3.0, 2.2, 0.5), mgp = c(1.8, 0.5, 0))
          graphics::plot(NULL, xlim = range(n_vals), ylim = ylim,
                         xlab = "Sample Size", ylab = "Bias",
                         main = term_title, log = "x", xaxt = "n", bty = "l")
          graphics::axis(1, at = n_vals, labels = n_vals, las = 2L,
                         mgp = c(1.7, 0.35, 0))
          graphics::abline(h = 0, lty = 2L, col = "gray40")
          draw_method_lines(term_data, "bias")

          # ---- MSE ----
          ylim <- range(term_data[est_rows, "mse"], na.rm = TRUE)
          graphics::par(mar = c(3.2, 3.0, 2.2, 0.5), mgp = c(1.8, 0.5, 0))
          graphics::plot(NULL, xlim = range(n_vals), ylim = ylim,
                         xlab = "Sample Size", ylab = "MSE",
                         main = "", log = "x", xaxt = "n", bty = "l")
          graphics::axis(1, at = n_vals, labels = n_vals, las = 2L,
                         mgp = c(1.7, 0.35, 0))
          draw_method_lines(term_data, "mse")
        }  # end term loop (estimator PDF)

        # ---- Bottom legend (flat/horizontal) ----
        graphics::par(mar = c(0.2, 0.2, 0.2, 0.2))
        graphics::plot.new()
        graphics::legend(
          "center",
          legend = paste0(avail_methods, "    "),
          col    = method_cols[avail_methods],
          pch    = 16L, lwd = 2L, lty = 1L,
          bty    = "n", cex = 0.9,
          ncol   = length(avail_methods),
          title  = "CI method", title.font = 2L
        )
        grDevices::dev.off()
        message("Saved: ", fig_path_est)

        # ======================================================
        # PDF 2: CI comparison (SE cal + stacked coverage + Width)
        # Two sub-rows per term (upper/lower cov stacked), legend at bottom
        # Layout per term:
        #   sub-row 1: [SE(span), UpperCov, Width(span)]
        #   sub-row 2: [SE(span), LowerCov, Width(span)]
        # Drawing order per layout reading: SE → UpperCov → Width → LowerCov
        # ======================================================
        fig_path_ci <- file.path(
          figures.dir,
          paste0("compare_", mtype, "_", vtype, "_", tbl_name, ".pdf")
        )
        n_rows_ci  <- 2L * n_terms + 1L
        lay_mat_ci <- matrix(0L, nrow = n_rows_ci, ncol = 3L)
        panel_id   <- 1L
        for (ti in seq_len(n_terms)) {
          r1 <- 2L * ti - 1L; r2 <- 2L * ti
          se_id   <- panel_id
          ucov_id <- panel_id + 1L
          wid_id  <- panel_id + 2L
          lcov_id <- panel_id + 3L
          lay_mat_ci[r1, 1L] <- se_id;   lay_mat_ci[r2, 1L] <- se_id   # SE spans
          lay_mat_ci[r1, 2L] <- ucov_id                                  # upper cov
          lay_mat_ci[r2, 2L] <- lcov_id                                  # lower cov
          lay_mat_ci[r1, 3L] <- wid_id;  lay_mat_ci[r2, 3L] <- wid_id  # Width spans
          panel_id <- panel_id + 4L
        }
        leg_id_ci <- panel_id
        lay_mat_ci[n_rows_ci, ] <- leg_id_ci

        sub_h <- row_h / 2   # height of each coverage sub-row (half of a full row)
        grDevices::pdf(fig_path_ci, width = 11,
                       height = row_h * n_terms + leg_h)
        graphics::layout(lay_mat_ci,
                         widths  = c(3.5, 3.0, 3.0),
                         heights = c(rep(c(sub_h, sub_h), n_terms), leg_h))

        for (ti in seq_along(terms)) {
          term      <- terms[[ti]]
          term_data <- sub[sub$term == term, ]
          true_s    <- get_true_S(mtype, vtype, tbl_name, term)
          term_title <- if (is.finite(true_s))
            sprintf("%s: S=%.3f", term, true_s) else term

          # Shared coverage y-range (same scale for upper + lower)
          cov_rows <- if (mtype == "glm") term_data$n >= 500 else
            rep(TRUE, nrow(term_data))
          cov_vals <- c(term_data[cov_rows, "upper_coverage"],
                        term_data[cov_rows, "lower_coverage"])
          cov_ylim <- range(c(cov_vals, 1 - alpha/2 - 0.02, 1.01), na.rm = TRUE)

          # Width y-range
          wid_rows <- if (mtype == "glm") term_data$n >= 500 else
            rep(TRUE, nrow(term_data))
          wid_vals <- term_data[wid_rows, "width"]
          wid_ylim <- range(wid_vals[wid_vals <= 10], na.rm = TRUE)

          # Drawing order: SE (spans) → UpperCov → Width (spans) → LowerCov

          # ---- Panel 1: SE calibration (spans both sub-rows, term title here) ----
          se_comp <- do.call(rbind, lapply(avail_methods, function(meth) {
            md <- term_data[term_data$ci_method == meth, ]
            md <- md[order(md$n), ]
            data.frame(
              ci_method  = meth,
              n          = md$n,
              emp_se     = sqrt(md$n) * md$empirical_sd,
              half_width = sqrt(md$n) * md$width / (2 * z_ref),
              stringsAsFactors = FALSE
            )
          }))
          se_vals <- c(se_comp$emp_se, se_comp$half_width)
          lim_se  <- range(se_vals[se_vals <= 10], na.rm = TRUE)
          graphics::par(mar = c(3.2, 3.2, 2.2, 0.5), mgp = c(1.8, 0.5, 0))
          graphics::plot(
            NULL, xlim = lim_se, ylim = lim_se,
            xlab = expression(sqrt(n) %*% " Empirical SE"),
            ylab = expression(sqrt(n) %*% " CI Half-Width / " * z[alpha/2]),
            main = term_title, bty = "l", asp = 1
          )
          graphics::abline(0, 1, lty = 2L, col = "gray40")
          for (meth in avail_methods) {
            md_se <- se_comp[se_comp$ci_method == meth, ]
            md_se <- md_se[order(md_se$n), ]
            if (nrow(md_se) == 0L) next
            graphics::points(md_se$emp_se[1L], md_se$half_width[1L],
                             col = method_cols[[meth]], pch = 16L, cex = 1.1)
            if (nrow(md_se) >= 2L) {
              for (i in seq_len(nrow(md_se) - 1L)) {
                graphics::arrows(
                  md_se$emp_se[i],     md_se$half_width[i],
                  md_se$emp_se[i+1L],  md_se$half_width[i+1L],
                  col = method_cols[[meth]], lwd = 1.5, length = 0.07, angle = 20L
                )
              }
            }
          }

          # ---- Panel 2: Upper coverage (top sub-row, no x-axis) ----
          graphics::par(mar = c(0.3, 3.0, 2.2, 0.5), mgp = c(1.8, 0.45, 0))
          graphics::plot(NULL, xlim = range(n_vals), ylim = cov_ylim,
                         xlab = "", ylab = "Upper Cov.",
                         main = "", log = "x", xaxt = "n", bty = "l")
          graphics::abline(h = 1 - alpha/2, lty = 2L, col = "gray40")
          draw_method_lines(term_data, "upper_coverage")

          # ---- Panel 3: Width (spans both sub-rows) ----
          graphics::par(mar = c(3.2, 3.0, 2.2, 0.5), mgp = c(1.8, 0.5, 0))
          graphics::plot(NULL, xlim = range(n_vals), ylim = wid_ylim,
                         xlab = "Sample Size", ylab = "CI Width",
                         main = "", log = "x", xaxt = "n", bty = "l")
          graphics::axis(1, at = n_vals, labels = n_vals, las = 2L,
                         mgp = c(1.7, 0.35, 0))
          draw_method_lines(term_data, "width")

          # ---- Panel 4: Lower coverage (bottom sub-row, with x-axis) ----
          graphics::par(mar = c(2.8, 3.0, 0.3, 0.5), mgp = c(1.8, 0.45, 0))
          graphics::plot(NULL, xlim = range(n_vals), ylim = cov_ylim,
                         xlab = "Sample Size", ylab = "Lower Cov.",
                         main = "", log = "x", xaxt = "n", bty = "l")
          graphics::axis(1, at = n_vals, labels = n_vals, las = 2L,
                         mgp = c(1.7, 0.35, 0))
          graphics::abline(h = 1 - alpha/2, lty = 2L, col = "gray40")
          draw_method_lines(term_data, "lower_coverage")
        }  # end term loop (CI PDF)

        # ---- Bottom legend (flat/horizontal) ----
        graphics::par(mar = c(0.2, 0.2, 0.2, 0.2))
        graphics::plot.new()
        graphics::legend(
          "center",
          legend = paste0(avail_methods, "    "),
          col    = method_cols[avail_methods],
          pch    = 16L, lwd = 2L, lty = 1L,
          bty    = "n", cex = 0.9,
          ncol   = length(avail_methods),
          title  = "CI method", title.font = 2L
        )
        grDevices::dev.off()
        message("Saved: ", fig_path_ci)
      }
    }
  }

  # ============================================================
  #  Coverage-quantile figures
  #  (S_true - LCI) / (UCI - LCI) density curves
  #  One PDF per (model × vcov × table): rows = terms, columns = n values
  # ============================================================
  message("Generating coverage-quantile figures...")

  insurance    <- RESI::insurance
  ci_lo_col_cq <- paste0(alpha / 2 * 100, "%")
  ci_hi_col_cq <- paste0((1 - alpha / 2) * 100, "%")

  if (fixed.knots) {
    .age_knots    <- stats::quantile(insurance$age, c(1/3, 2/3))
    .age_bk       <- range(insurance$age)
    lm_formula_cq  <- eval(bquote(
      log10(charges) ~ splines::ns(age, knots = .(.age_knots),
                                   Boundary.knots = .(.age_bk)) *
        sex + bmi + smoker + region))
    glm_formula_cq <- eval(bquote(
      I(charges > 15000) ~ splines::ns(age, knots = .(.age_knots),
                                       Boundary.knots = .(.age_bk)) *
        sex + bmi + smoker + region))
  } else {
    lm_formula_cq  <- log10(charges) ~ splines::ns(age, df = 3) *
      sex + bmi + smoker + region
    glm_formula_cq <- I(charges > 15000) ~ splines::ns(age, df = 3) *
      sex + bmi + smoker + region
  }

  cq_settings <- list(
    list(type = "lm",  key = "lm_parametric",
         formula = lm_formula_cq,  family = NULL,
         vcovfunc = stats::vcov,         robust = FALSE),
    list(type = "lm",  key = "lm_robust",
         formula = lm_formula_cq,  family = NULL,
         vcovfunc = sandwich::vcovHC,    robust = TRUE),
    list(type = "glm", key = "glm_parametric",
         formula = glm_formula_cq, family = stats::binomial(),
         vcovfunc = stats::vcov,         robust = FALSE),
    list(type = "glm", key = "glm_robust",
         formula = glm_formula_cq, family = stats::binomial(),
         vcovfunc = sandwich::vcovHC,    robust = TRUE)
  )

  # Compute true RESI from full insurance dataset
  cq_true <- lapply(cq_settings, function(s) {
    full_mod <- tryCatch({
      if (s$type == "lm") {
        m <- stats::lm(s$formula, data = insurance)
        m$call[["formula"]] <- s$formula; m
      } else {
        m <- stats::glm(s$formula, data = insurance, family = s$family)
        m$call[["formula"]] <- s$formula
        m$call[["family"]]  <- s$family; m
      }
    }, error = function(e) NULL)
    if (is.null(full_mod)) return(NULL)
    # Use HC0 for robust (no hat-value correction at full n)
    tv_vcov <- if (s$robust) function(x) sandwich::vcovHC(x, type = "HC0")
               else s$vcovfunc
    tryCatch(resi_pe(full_mod, vcovfunc = tv_vcov, unbiased = TRUE),
             error = function(e) NULL)
  })
  names(cq_true) <- vapply(cq_settings, `[[`, character(1L), "key")

  for (mtype_cq in c("lm", "glm")) {
    for (vtype_cq in c("parametric", "robust")) {
      key_cq  <- paste0(mtype_cq, "_", vtype_cq)
      pe_true <- cq_true[[key_cq]]
      if (is.null(pe_true)) next

      for (tbl_cq in c("anova", "coefficients")) {
        true_t <- pe_true[[tbl_cq]]
        if (is.null(true_t) || nrow(true_t) == 0L) next

        terms_cq <- rownames(true_t)
        terms_cq <- setdiff(terms_cq, c("(Intercept)", "Residuals"))
        if (length(terms_cq) == 0L) next

        n_vals_cq <- sort(unique(combined$n))

        # Load per-replicate (LCI, UCI) for every (method, n, term)
        cq_dat <- lapply(setNames(avail_methods, avail_methods), function(meth) {
          raw_dir <- file.path(output.dirs[[meth]], "sim_raw")
          if (!dir.exists(raw_dir)) return(NULL)
          lapply(setNames(n_vals_cq, as.character(n_vals_cq)), function(n_val) {
            f <- file.path(raw_dir, paste0(key_cq, "_n", n_val, ".rds"))
            if (!file.exists(f)) return(NULL)
            reps <- tryCatch(readRDS(f), error = function(e) NULL)
            if (is.null(reps) || length(reps) == 0L) return(NULL)
            lapply(setNames(terms_cq, terms_cq), function(term) {
              s_true <- true_t[term, "RESI"]
              lci_v <- vapply(reps, function(r) {
                tbl_r <- r[[tbl_cq]]
                if (is.null(tbl_r) || !(term %in% rownames(tbl_r))) return(NA_real_)
                if (!(ci_lo_col_cq %in% colnames(tbl_r)))           return(NA_real_)
                tbl_r[term, ci_lo_col_cq]
              }, numeric(1L))
              uci_v <- vapply(reps, function(r) {
                tbl_r <- r[[tbl_cq]]
                if (is.null(tbl_r) || !(term %in% rownames(tbl_r))) return(NA_real_)
                if (!(ci_hi_col_cq %in% colnames(tbl_r)))           return(NA_real_)
                tbl_r[term, ci_hi_col_cq]
              }, numeric(1L))
              cq_v <- (s_true - lci_v) / (uci_v - lci_v)
              cq_v[!is.finite(cq_v)] <- NA_real_
              cq_v
            })
          })
        })

        # Build PDF
        n_terms_cq  <- length(terms_cq)
        n_sizes_cq  <- length(n_vals_cq)
        n_panels_cq <- n_terms_cq * n_sizes_cq
        lay_cq      <- matrix(seq_len(n_panels_cq),
                              nrow = n_terms_cq, ncol = n_sizes_cq,
                              byrow = TRUE)
        lay_cq      <- cbind(lay_cq, n_panels_cq + 1L)

        col_w_cq <- 2.4
        row_h_cq <- 2.4

        fig_path_cq <- file.path(
          figures.dir,
          paste0("covquant_", mtype_cq, "_", vtype_cq, "_", tbl_cq, ".pdf")
        )
        grDevices::pdf(fig_path_cq,
                       width  = col_w_cq * n_sizes_cq + 2.2,
                       height = row_h_cq * n_terms_cq + 0.5)
        graphics::layout(lay_cq,
                         widths  = c(rep(col_w_cq, n_sizes_cq), 2.2),
                         heights = rep(row_h_cq, n_terms_cq))

        for (ti in seq_along(terms_cq)) {
          term_cq      <- terms_cq[[ti]]
          is_first_row <- (ti == 1L)
          is_last_row  <- (ti == n_terms_cq)

          for (ni in seq_along(n_vals_cq)) {
            n_cq         <- n_vals_cq[[ni]]
            is_first_col <- (ni == 1L)

            # Compute density for each method
            dens_list_cq <- lapply(setNames(avail_methods, avail_methods), function(meth) {
              vals <- cq_dat[[meth]][[as.character(n_cq)]][[term_cq]]
              vals <- vals[is.finite(vals)]
              if (length(vals) < 5L) return(NULL)
              tryCatch(stats::density(vals, n = 256L), error = function(e) NULL)
            })

            xlim_cq <- {
              all_x <- unlist(c(
                lapply(dens_list_cq, function(d) if (!is.null(d)) range(d$x)),
                list(c(0, 1))))
              c(min(all_x) - 0.05, max(all_x) + 0.05)
            }
            max_y_cq <- max(
              vapply(dens_list_cq,
                     function(d) if (!is.null(d)) max(d$y) else 0, numeric(1L)),
              1.0)
            ylim_cq <- c(0, max_y_cq * 1.1)

            t_mar <- if (is_first_row) 1.8 else 0.4
            b_mar <- if (is_last_row)  2.5 else 0.4
            l_mar <- if (is_first_col) 3.2 else 1.2

            graphics::par(mar = c(b_mar, l_mar, t_mar, 0.3), mgp = c(1.5, 0.4, 0))
            graphics::plot(NULL, xlim = xlim_cq, ylim = ylim_cq,
                           xlab = if (is_last_row)  "(S_true - LCI) / Width" else "",
                           ylab = if (is_first_col) "Density" else "",
                           main = if (is_first_row) paste0("n = ", n_cq) else "",
                           bty  = "l",
                           xaxt = if (is_last_row)  "s" else "n",
                           yaxt = if (is_first_col) "s" else "n")

            # Row label (term name) on left margin of first column
            if (is_first_col)
              graphics::mtext(term_cq, side = 2L, line = 2.0,
                              cex = 0.65, las = 0, font = 2L)

            # Reference lines: 0/1 = CI boundaries, 0.5 = midpoint, y=1 = uniform
            graphics::abline(v = c(0, 1), lty = 2L, col = "gray50", lwd = 0.8)
            graphics::abline(v = 0.5,     lty = 3L, col = "gray70", lwd = 0.8)
            graphics::abline(h = 1.0,     lty = 3L, col = "gray80", lwd = 0.8)

            for (meth in avail_methods) {
              d <- dens_list_cq[[meth]]
              if (is.null(d)) next
              graphics::lines(d, col = method_cols[[meth]], lwd = 1.8)
            }
          }
        }

        # Shared legend
        graphics::par(mar = c(0.3, 0.3, 0.3, 0.3))
        graphics::plot.new()
        graphics::legend("center",
                         legend = avail_methods,
                         col    = method_cols[avail_methods],
                         lwd = 2L, lty = 1L,
                         bty = "n", cex = 0.9,
                         title = "CI method", title.font = 2L)

        grDevices::dev.off()
        message("Saved: ", fig_path_cq)
      }
    }
  }

  message("Comparison figures saved to: ", figures.dir)
  invisible(combined)
}
