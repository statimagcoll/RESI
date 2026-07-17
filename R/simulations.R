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

  # upper_coverage: P(hi >= tv) -- CI upper bound is at or above the true value
  # lower_coverage: P(lo <= tv) -- CI lower bound is at or below the true value
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
      message("No raw files found for: ", lbl, " -- skipping")
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

      # Layout: n_terms rows x 2 cols + 1 legend row
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

      # Layout: 4 rows x 5 cols; coverage column (col 3) is split into two
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

# Internal: full delta-method sigma2S (including A/B chain-rule terms) for each
# coefficient in coef_idx (1-based index into coef(model)), using HC-type B.
# Returns named vector sigma2S_Th1_k = quad_k / S_k^2  (= n * resiSE_k^2).
.sigma2S_Th1 <- function(model, coef_idx, type = "HC0") {
  is_lm  <- inherits(model, "lm") && !inherits(model, "glm")
  X      <- model.matrix(model)
  if (any(al <- is.na(coef(model)))) X <- X[, !al, drop = FALSE]
  ef     <- sandwich::estfun(model)
  ef     <- ef[, colnames(X), drop = FALSE]
  p  <- ncol(X); n <- nrow(X)
  e  <- residuals(model, "response")
  sym_fn <- function(M) (M + t(M)) / 2

  if (is_lm) {
    phi   <- summary(model)$sigma^2
    m     <- p + 1L
    theta <- c(phi, coef(model))
    psi_list <- lapply(seq_len(n), function(i) {
      xi <- X[i, , drop = FALSE]; ei <- e[i]
      cbind((ei^2 - phi) / (2*phi^2), matrix(ei * xi / phi, 1))
    })
    psiprime_list <- lapply(seq_len(n), function(i) {
      xi <- X[i, , drop = FALSE]; ei <- e[i]
      rbind(c((phi - 2*ei^2) / (2*phi^3), -ei * as.vector(xi) / phi^2),
            cbind(matrix(-ei * as.vector(xi) / phi^2, p, 1),
                  -crossprod(xi) / phi))
    })
  } else {
    m     <- p
    theta <- coef(model)
    w_vec <- weights(model, type = "working")
    psi_list <- lapply(seq_len(n), function(i)
      matrix(ef[i, , drop = TRUE], 1))
    psiprime_list <- lapply(seq_len(n), function(i)
      -w_vec[i] * crossprod(X[i, , drop = FALSE]))
  }

  A_full <- sym_fn(-Reduce("+", psiprime_list) / n)
  A_inv  <- tryCatch(chol2inv(chol(A_full)), error = function(e2) solve(A_full))

  h    <- hatvalues(model)
  omh  <- pmax(1 - h, .Machine$double.eps)
  sqrtw <- switch(toupper(type),
    "HC0" = rep(1, n),
    "HC1" = rep(sqrt(n / max(n - p, 1L)), n),
    "HC2" = 1 / sqrt(omh),
    "HC3" = 1 / omh,
    rep(1, n)
  )
  B_full <- sym_fn(Reduce("+", lapply(seq_len(n), function(i) {
    psi_i <- psi_list[[i]]
    if (is_lm) psi_i[1L, seq(2L, m)] <- psi_i[1L, seq(2L, m)] * sqrtw[i]
    else       psi_i[1L, ]           <- psi_i[1L, ] * sqrtw[i]
    crossprod(psi_i)
  })) / n)
  cov_th <- sym_fn(A_inv %*% B_full %*% A_inv)

  XtXn    <- crossprod(X) / n
  XtXn_hc <- crossprod(X * sqrtw) / n

  if (is_lm) {
    dA_phi <- rbind(
      c(-1 / phi^3, numeric(p)),
      cbind(numeric(p), -XtXn / phi^2))
    dA_dth <- matrix(0, m * m, m); dA_dth[, 1L] <- as.vector(dA_phi)
    dB_phi <- rbind(
      c((phi^2 - 2 * mean(e^4)) / phi^5,
        -(3 / (2 * phi^4)) * colMeans(X) * mean(e^3)),
      cbind(matrix(-(3 / (2 * phi^4)) * colMeans(X) * mean(e^3), p, 1),
            -XtXn_hc / phi^2))
    dB_dth <- matrix(0, m * m, m); dB_dth[, 1L] <- as.vector(dB_phi)
  } else {
    mu_hat <- fitted(model)
    dA_dth <- matrix(0, m * m, m)
    for (j in seq_len(m)) {
      vj <- w_vec * (1 - 2 * mu_hat) * X[, j]
      dA_dth[, j] <- as.vector(crossprod(X, vj * X) / n)
    }
    dB_dth <- dA_dth
  }

  mcov_n <- n * sandwich::vcovHC(model, type = type)
  nms    <- names(coef(model))[coef_idx]
  result <- setNames(numeric(length(coef_idx)), nms)

  for (ki in seq_along(coef_idx)) {
    k <- coef_idx[ki]
    if (is_lm) {
      L  <- matrix(0, 1, m); L[1L, k + 1L] <- 1  # +1: phi occupies position 1
      Lm <- matrix(0, 1, p); Lm[1L, k]     <- 1
    } else {
      L  <- matrix(0, 1, m); L[1L, k] <- 1
      Lm <- L
    }
    beta_k   <- as.numeric(L %*% theta)
    cov_beta <- as.numeric(L %*% cov_th %*% t(L))
    Ssq      <- beta_k^2 / as.numeric(Lm %*% mcov_n %*% t(Lm))
    if (!is.finite(Ssq) || Ssq <= 0 || !is.finite(cov_beta) || cov_beta <= 0) {
      result[ki] <- NA_real_; next
    }
    La_row <- L %*% A_inv  # 1xm
    Lc_row <- L %*% cov_th # 1xm
    sc <- beta_k / cov_beta
    deriv_th <- sc * L
    deriv_A  <- (sc^2 / 2) *
      (kronecker(La_row, Lc_row) + kronecker(Lc_row, La_row)) %*% dA_dth
    deriv_B  <- -(sc^2 / 2) * kronecker(La_row, La_row) %*% dB_dth
    d_tot    <- deriv_th + deriv_A + deriv_B
    quad     <- max(0, as.numeric(d_tot %*% cov_th %*% t(d_tot)))
    result[ki] <- quad / Ssq
  }
  result
}

# ============================================================
#  simCalibrationFigures
# ============================================================

#' Asymptotic Calibration Check for RESI Variance Estimates
#'
#' Runs a plasmode simulation to verify that the asymptotic normal CI machinery
#' is correctly calibrated for all four model settings (lm/glm x
#' parametric/robust).  For each (setting, sample size) cell the function
#' checks:
#' \enumerate{
#'   \item \strong{Bias(theta)}: \eqn{\hat\theta \to \theta_{\rm true}} --
#'     raw coefficients (and \eqn{\phi = \hat\sigma^2} for \code{lm}) converge
#'     to the full-dataset values.
#'   \item \strong{vcov check A} (estimator consistency):
#'     \eqn{n \cdot \bar V_{\rm analytic} / (n_{\rm full} \cdot V_{\rm true}) \to 1}.
#'   \item \strong{vcov check B} (calibration):
#'     \eqn{n \cdot \widehat{\rm Var}(\hat\beta) / (n_{\rm full} \cdot V_{\rm true}) \to 1}.
#'   \item \strong{Bias(R)}: RESI point estimates converge to the full-dataset
#'     values.
#'   \item \strong{sigma2S check A} (estimator consistency):
#'     \eqn{\bar{\hat\sigma}^2_S / \sigma^2_{S,\rm true} \to 1}.
#'   \item \strong{sigma2S check B} (CI calibration):
#'     \eqn{n \cdot \widehat{\rm Var}(\hat R) / \sigma^2_{S,\rm true} \to 1}.
#' }
#' True values are taken from the full \code{\link{insurance}} dataset using the
#' same definitions as \code{\link{insurancePlasmodeSim}}.  \eqn{\sigma^2_S} is
#' extracted from the asymptotic-normal CI half-width:
#' \eqn{\hat\sigma^2_S = n \cdot ({\rm hw}/z_{\alpha/2})^2}.
#'
#' @param nsim Integer, replicates per (setting, \eqn{n}) cell. Default 500.
#' @param n.vec Integer vector of sample sizes.
#'   Default \code{c(100, 200, 500, 1000, 2000, 5000)}.
#' @param alpha Numeric, nominal CI level. Default 0.05.
#' @param output.dir Character, directory for raw per-cell RDS files.
#'   Default \code{"resiCalibrationSim"}.
#' @param fixed.knots Logical. Fix spline knots at full-dataset quantiles.
#'   Default \code{FALSE}.
#' @param mc.cores.reps Integer, cores for within-cell parallelism. Default 1.
#'
#' @return Invisibly returns the combined metrics \code{data.frame}.
#' @seealso \code{\link{insurancePlasmodeSim}}, \code{\link{simCompareMethodsFigures}}
#' @importFrom parallel mclapply
#' @importFrom splines ns
#' @importFrom sandwich vcovHC
#' @importFrom stats lm glm vcov binomial qnorm coef setNames sd
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot lines points abline legend par axis plot.new layout
#' @export
simCalibrationSim <- function(
    nsim          = 500L,
    n.vec         = c(100L, 200L, 500L, 1000L, 2000L, 5000L),
    alpha         = 0.05,
    output.dir    = "resiCalibrationSim",
    fixed.knots   = FALSE,
    deriv_method  = c("corrected", "original", "population", "population2", "zeroB", "indep_moment"),
    mc.cores.reps = 1L
) {
  deriv_method <- match.arg(deriv_method)
  raw_dir <- file.path(output.dir, "raw", deriv_method)
  dir.create(raw_dir,      recursive = TRUE, showWarnings = FALSE)

  insurance <- RESI::insurance
  n_full    <- nrow(insurance)
  z_ref     <- stats::qnorm(1 - alpha / 2)
  ci_lo_col <- paste0(alpha / 2 * 100,       "%")
  ci_hi_col <- paste0((1 - alpha / 2) * 100, "%")

  .fvars   <- c("sex", "smoker", "region")
  .flevels <- lapply(.fvars, function(v) unique(insurance[[v]]))
  names(.flevels) <- .fvars

  # ---- Formulas ---------------------------------------------------------------
  if (fixed.knots) {
    .age_knots <- stats::quantile(insurance$age, c(1/3, 2/3))
    .age_bk    <- range(insurance$age)
    lm_formula  <- eval(bquote(log10(charges) ~
      splines::ns(age, knots = .(.age_knots), Boundary.knots = .(.age_bk)) *
        sex + bmi + smoker + region))
    glm_formula <- eval(bquote(I(charges > 15000) ~
      splines::ns(age, knots = .(.age_knots), Boundary.knots = .(.age_bk)) *
        sex + bmi + smoker + region))
  } else {
    lm_formula  <- log10(charges) ~
      splines::ns(age, df = 3) * sex + bmi + smoker + region
    glm_formula <- I(charges > 15000) ~
      splines::ns(age, df = 3) * sex + bmi + smoker + region
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

  # ---- True values from full dataset ------------------------------------------
  message("Computing true values from full dataset (n = ", n_full, ")...")
  true_vals <- lapply(model_settings, function(s) {
    .formula <- s$formula; .family <- s$family
    full_mod <- if (s$type == "lm") {
      m <- stats::lm(.formula, data = insurance)
      m$call[["formula"]] <- .formula; m
    } else {
      m <- stats::glm(.formula, data = insurance, family = .family)
      m$call[["formula"]] <- .formula
      m$call[["family"]]  <- .family; m
    }

    # theta_true: phi = sigma^2 = summary()$sigma^2 (same as .resi_precompute)
    phi_true  <- if (s$type == "lm") summary(full_mod)$sigma^2 else NULL
    beta_true <- stats::coef(full_mod)

    # n_full * vcov_true: full covariance matrix (all non-aliased coefficients)
    vcov_full_true_raw <- n_full * s$vcovfunc(full_mod)
    vcov_diag_true     <- diag(vcov_full_true_raw)
    # Submatrix for non-intercept terms (computed later once coef_terms_true is known)

    # R_true using the same population definition as insurancePlasmodeSim:
    #   - robust: HC0 (no hat-value correction)
    #   - parametric lm: divide by sqrt(rdf) = sqrt(n-p), not sqrt(n)
    true_vcovfunc_r <- if (s$vcov_name == "robust") {
      function(x) sandwich::vcovHC(x, type = "HC0")
    } else {
      s$vcovfunc
    }
    rdf_true_r <- if (s$type == "lm" && s$vcov_name == "parametric") {
      full_mod$df.residual
    } else {
      NULL
    }
    pe_true_raw <- resi_pe(full_mod, data = insurance, vcovfunc = true_vcovfunc_r)
    pe_true     <- .simDirectRESI(pe_true_raw, n = n_full, rdf = rdf_true_r)
    coef_tab_true   <- pe_true$coefficients
    coef_terms_true <- rownames(coef_tab_true)[rownames(coef_tab_true) != "(Intercept)"]
    R_true <- coef_tab_true[coef_terms_true, "RESI"]
    names(R_true) <- coef_terms_true   # data.frame[rows,col] drops names

    # sigma2S_true: Sigma_R[k,k] from the delta-method variance of R_hat,
    #   evaluated at the population values (full-model coefficients + HC0 Sigma_theta).
    #   HC0 is used as the "true" population covariance (consistent estimator).
    #   For parametric: vcov_is_model=TRUE triggers the parametric derivative;
    #     type="HC0" supplies HC0-based Sigma_theta via the B matrix.
    #   For robust: vcovfunc=HC0 triggers the robust derivative; type="HC0" same.
    #   sigma2S_true_k = [dR_dtheta * Sigma_theta_HC0 * dR_dtheta^T]_{kk}
    #                  = n_full * (hw / z_ref)^2 from the normal CI
    asym_true <- tryCatch(
      resi_pe_asymptotic(full_mod, vcovfunc = true_vcovfunc_r,
                         ci.method = "normal", type = "HC0",
                         deriv_method = deriv_method),
      error = function(e) NULL
    )
    if (!is.null(asym_true)) {
      tab_a <- asym_true$coefficients
      tab_a <- tab_a[rownames(tab_a) != "(Intercept)", , drop = FALSE]
      hw_true <- (tab_a[, ci_hi_col] - tab_a[, ci_lo_col]) / 2
      names(hw_true) <- rownames(tab_a)   # data.frame[,col] drops names
      sigma2S_true <- n_full * (hw_true / z_ref)^2
    } else {
      sigma2S_true <- setNames(rep(NA_real_, length(R_true)), names(R_true))
    }

    # sigma2S_true_Th1: full delta-method variance (including A/B chain rules)
    #   evaluated at population values with HC0 Sigma_theta.
    #   sigma2S_Th1_k = quad_k / S_k^2 where quad_k uses deriv_theta + deriv_A + deriv_B.
    coef_idx_th1 <- match(coef_terms_true, names(coef(full_mod)))
    sigma2S_true_Th1 <- tryCatch(
      .sigma2S_Th1(full_mod, coef_idx_th1, type = "HC0"),
      error = function(e) setNames(rep(NA_real_, length(coef_terms_true)),
                                    coef_terms_true)
    )

    # vcov_phi_phi_true and vcov_phi_beta_true: sandwich (phi,phi) and (phi,beta_k)
    # elements for the lm settings, using the per-setting HC type as the "true"
    # population reference (so the check targets convergence of the sample
    # estimator to the full-dataset sandwich, consistent with vcov_diag_true).
    # For GLM there is no phi: set to NA.
    precomp_type_tv  <- if (s$vcov_name == "robust") "HC3" else "const"
    precomp_tv       <- NULL   # initialise; assigned below for lm only
    vcov_phi_phi_true  <- NA_real_
    vcov_phi_beta_true <- setNames(rep(NA_real_, length(coef_terms_true)),
                                   coef_terms_true)
    if (s$type == "lm") {
      precomp_tv <- tryCatch(
        .resi_precompute(full_mod, type = precomp_type_tv),
        error = function(e) NULL
      )
      if (!is.null(precomp_tv)) {
        vcov_phi_phi_true <- precomp_tv$cov_theta[1L, 1L]
        beta_names_tv     <- names(coef(full_mod))[!is.na(coef(full_mod))]
        beta_pos_tv       <- match(coef_terms_true, beta_names_tv)
        vcov_phi_beta_true <- precomp_tv$cov_theta[1L, 1L + beta_pos_tv]
        names(vcov_phi_beta_true) <- coef_terms_true
      }
    } else {
      # GLM: precomp_tv needed for vcov_full_true_mat (beta-beta block)
      precomp_tv <- tryCatch(
        .resi_precompute(full_mod, type = precomp_type_tv),
        error = function(e) NULL
      )
    }

    # Chain decomposition of Sigma_R at true-value estimates (HC0 sandwich).
    # Jacobian pieces: J = J_direct + J_Achain + J_Bchain  (each 1 x m for m1=1).
    # Analytic: sig2S_true_X_k = J_X %*% Sig_HC0 %*% t(J_X)  for X in {dir,Ach,Bch}.
    # MC analog: n_s * var(J_X_TV %*% (theta_hat - theta_true))  per replicate.
    # HC0 is used for consistency with sigma2S_true (also HC0).
    precomp_hc0 <- tryCatch(
      .resi_precompute(full_mod, type = "HC0", deriv_method = deriv_method),
      error = function(e) NULL
    )
    na_terms <- setNames(rep(NA_real_, length(coef_terms_true)), coef_terms_true)
    J_dir_list    <- list()
    J_Achain_list <- list()
    J_Bchain_list <- list()
    sig2S_true_dir    <- na_terms
    sig2S_true_Achain <- na_terms
    sig2S_true_Bchain <- na_terms
    theta_true_full   <- NULL
    if (!is.null(precomp_hc0)) {
      Sig_hc0 <- precomp_hc0$cov_theta   # m x m Sigma_theta at HC0
      theta_true_full <- precomp_hc0$theta_hat  # (phi, beta_1,...,beta_p)
      for (tm in coef_terms_true) {
        ct_tm <- tryCatch(.resi_contrast(precomp_hc0,
                                         .get_L_coef(full_mod, tm)),
                          error = function(e) NULL)
        if (!is.null(ct_tm)) {
          J_dir_list[[tm]]    <- ct_tm$dR_direct
          J_Achain_list[[tm]] <- ct_tm$dR_Achain
          J_Bchain_list[[tm]] <- ct_tm$dR_Bchain
          sig2S_true_dir[tm]    <- as.numeric(ct_tm$dR_direct  %*% Sig_hc0 %*% t(ct_tm$dR_direct))
          sig2S_true_Achain[tm] <- as.numeric(ct_tm$dR_Achain %*% Sig_hc0 %*% t(ct_tm$dR_Achain))
          sig2S_true_Bchain[tm] <- as.numeric(ct_tm$dR_Bchain %*% Sig_hc0 %*% t(ct_tm$dR_Bchain))
        }
      }
    }

    list(deriv_method_used   = deriv_method,
         phi_true            = phi_true,
         beta_true           = beta_true,
         vcov_diag_true      = vcov_diag_true,
         # vcov_full_true_mat: full analytic Sigma_theta from precomp_tv (HC3 or const).
         # For lm: (p+1)x(p+1) with phi as first row/col; for glm: pxp.
         # Row/col names: c("phi", coef_terms_true) for lm, coef_terms_true for glm.
         vcov_full_true_mat  = if (!is.null(precomp_tv)) {
           if (s$type == "lm") {
             m_all  <- precomp_tv$m
             beta_names_all <- names(coef(full_mod))[!is.na(coef(full_mod))]
             beta_pos_all   <- match(coef_terms_true, beta_names_all)
             idx_all <- c(1L, 1L + beta_pos_all)
             nm_all  <- c("phi", coef_terms_true)
             mat_all <- precomp_tv$cov_theta[idx_all, idx_all, drop = FALSE]
             rownames(mat_all) <- colnames(mat_all) <- nm_all
             mat_all
           } else {
             beta_names_all <- names(coef(full_mod))[!is.na(coef(full_mod))]
             beta_pos_all   <- match(coef_terms_true, beta_names_all)
             mat_all <- precomp_tv$cov_theta[beta_pos_all, beta_pos_all, drop = FALSE]
             rownames(mat_all) <- colnames(mat_all) <- coef_terms_true
             mat_all
           }
         } else NULL,
         vcov_phi_phi_true   = vcov_phi_phi_true,
         vcov_phi_beta_true  = vcov_phi_beta_true,
         R_true              = R_true,
         sigma2S_true        = sigma2S_true,
         sigma2S_true_Th1    = sigma2S_true_Th1,
         J_dir_list          = J_dir_list,
         J_Achain_list       = J_Achain_list,
         J_Bchain_list       = J_Bchain_list,
         sig2S_true_dir      = sig2S_true_dir,
         sig2S_true_Achain   = sig2S_true_Achain,
         sig2S_true_Bchain   = sig2S_true_Bchain,
         theta_true_full     = theta_true_full,
         coef_terms          = coef_terms_true)
  })
  names(true_vals) <- sapply(model_settings, `[[`, "label")

  # ---- Per-setting, per-n simulations -----------------------------------------
  all_metrics <- list()

  for (s in model_settings) {
    lbl      <- s$label
    tv       <- true_vals[[lbl]]
    is_lm    <- (s$type == "lm")
    .formula <- s$formula; .family <- s$family
    coef_terms <- tv$coef_terms
    message("\nSetting: ", lbl)

    for (n_s in n.vec) {
      cell_label <- paste0(lbl, "_n", n_s)
      rds_path   <- file.path(raw_dir, paste0(cell_label, ".rds"))

      if (file.exists(rds_path)) {
        message("  n = ", n_s, " (loading cached)")
        cell_data <- readRDS(rds_path)
      } else {
        message("  n = ", n_s, " ...")
        reps_raw <- parallel::mclapply(seq_len(nsim), function(i) {
          repeat {
            dat <- insurance[sample(n_full, n_s, replace = TRUE), ]
            ok  <- all(vapply(.fvars, function(v)
              all(.flevels[[v]] %in% dat[[v]]), logical(1L)))
            if (ok && !is_lm) {
              biny <- as.integer(dat$charges > 15000)
              ok   <- length(unique(biny)) > 1L &&
                all(tapply(biny, dat$smoker,
                           function(x) length(unique(x)) > 1L))
            }
            if (ok) break
          }

          mod <- tryCatch({
            m <- if (is_lm) stats::lm(.formula, data = dat)
                 else stats::glm(.formula, data = dat, family = .family)
            m$call[["formula"]] <- .formula
            if (!is_lm) m$call[["family"]] <- .family
            m
          }, error = function(e) NULL)
          if (is.null(mod)) return(NULL)

          tryCatch({
            # theta_hat
            phi_hat  <- if (is_lm) summary(mod)$sigma^2 else NULL
            beta_hat <- stats::coef(mod)[coef_terms]

            # n_s * diag(vcov_hat) for non-intercept terms
            vcov_hat_diag <- n_s * diag(s$vcovfunc(mod))[coef_terms]

            # R_hat and sigma2S_hat via asymptotic normal CI
            asym <- resi_pe_asymptotic(mod, vcovfunc = s$vcovfunc,
                                       ci.method = "normal",
                                       deriv_method = deriv_method)
            tab  <- asym$coefficients
            tab  <- tab[rownames(tab) != "(Intercept)", , drop = FALSE]
            hw   <- (tab[coef_terms, ci_hi_col] - tab[coef_terms, ci_lo_col]) / 2
            names(hw)   <- coef_terms   # data.frame[rows,col] drops names
            R_hat       <- tab[coef_terms, "RESI"]
            names(R_hat) <- coef_terms
            sigma2S_hat <- n_s * (hw / z_ref)^2
            names(sigma2S_hat) <- coef_terms

            # phi-variance and phi-beta covariances from the full sandwich
            # (lm only; requires model-level precomputation)
            vcov_phi_phi  <- NULL
            vcov_phi_beta <- NULL
            if (is_lm) {
              precomp_type_rep <- if (s$vcov_name == "robust") "HC3" else "const"
              precomp_rep <- tryCatch(
                .resi_precompute(mod, type = precomp_type_rep,
                                 deriv_method = deriv_method),
                error = function(e) NULL
              )
              if (!is.null(precomp_rep)) {
                # cov_theta[1,1] = Sigma_theta_{phi,phi} = asymptotic variance of
                # sqrt(n)(phi_hat - phi).  Do NOT multiply by n_s: cov_theta is
                # already the asymptotic (scaled) variance, analogous to
                # n_s * diag(vcovfunc(mod)) for the beta terms.
                vcov_phi_phi <- precomp_rep$cov_theta[1L, 1L]
                beta_names_rep <- names(coef(mod))[!is.na(coef(mod))]
                beta_pos_rep   <- match(coef_terms, beta_names_rep)
                vcov_phi_beta  <- precomp_rep$cov_theta[1L, 1L + beta_pos_rep]
                names(vcov_phi_beta) <- coef_terms
              }
            }

            list(phi_hat        = phi_hat,
                 beta_hat       = beta_hat,
                 theta_hat_full = NULL,     # removed: Jacobian-projection MC not used
                 vcov_hat_diag  = vcov_hat_diag,
                 R_hat          = R_hat,
                 sigma2S_hat    = sigma2S_hat,
                 vcov_phi_phi   = vcov_phi_phi,
                 vcov_phi_beta  = vcov_phi_beta)
          }, error = function(e) NULL)
        }, mc.cores = mc.cores.reps)

        reps <- Filter(Negate(is.null), reps_raw)

        # Store matrices
        cell_data <- list(
          n           = n_s,
          n_success   = length(reps),
          phi_vec     = if (is_lm) sapply(reps, `[[`, "phi_hat") else NULL,
          beta_mat    = do.call(rbind, lapply(reps, `[[`, "beta_hat")),
          vcov_mat    = do.call(rbind, lapply(reps, `[[`, "vcov_hat_diag")),
          R_mat       = do.call(rbind, lapply(reps, `[[`, "R_hat")),
          sig2S_mat   = do.call(rbind, lapply(reps, `[[`, "sigma2S_hat")),
          phi_var_vec = if (is_lm && !is.null(reps[[1L]]$vcov_phi_phi))
                          sapply(reps, `[[`, "vcov_phi_phi") else NULL,
          phi_cov_mat = if (is_lm && !is.null(reps[[1L]]$vcov_phi_beta))
                          do.call(rbind, lapply(reps, `[[`, "vcov_phi_beta")) else NULL
        )
        saveRDS(cell_data, rds_path)
      }

      if (cell_data$n_success == 0L) next

      # ---- Compute metrics ---------------------------------------------------
      n_s_eff <- cell_data$n  # actual n used

      # (1) Bias(phi) for lm
      phi_bias <- if (!is.null(cell_data$phi_vec))
        mean(cell_data$phi_vec, na.rm = TRUE) - tv$phi_true else NA_real_

      # (2) Bias(beta_k)
      beta_bias <- colMeans(cell_data$beta_mat, na.rm = TRUE) -
        tv$beta_true[coef_terms]

      # (3) vcov check A: mean(n_s * vcov_hat_kk) / (n_full * vcov_true_kk)
      vcov_A <- colMeans(cell_data$vcov_mat, na.rm = TRUE) /
        tv$vcov_diag_true[coef_terms]

      # (4) vcov check B: n_s * Var_empirical(beta_k) / (n_full * vcov_true_kk)
      vcov_B <- (n_s_eff * apply(cell_data$beta_mat, 2L, stats::var, na.rm = TRUE)) /
        tv$vcov_diag_true[coef_terms]

      # (5) Bias(R_k)
      R_bias <- colMeans(cell_data$R_mat, na.rm = TRUE) - tv$R_true[coef_terms]

      # (6) sigma2S check A: mean(sigma2S_hat_k) / sigma2S_true_k
      sig2S_A <- colMeans(cell_data$sig2S_mat, na.rm = TRUE) /
        tv$sigma2S_true[coef_terms]

      # (7) sigma2S check B: n_s * Var_empirical(R_k) / sigma2S_true_k
      sig2S_B <- (n_s_eff * apply(cell_data$R_mat, 2L, stats::var, na.rm = TRUE)) /
        tv$sigma2S_true[coef_terms]

      # (8) Coefficient of variation of the variance estimator V_hat_kk across
      #     replicates -- quantifies the "variance of the variance estimator".
      cv_Vhat <- apply(cell_data$vcov_mat, 2L, function(x)
        stats::sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))

      # (9) phi variance check A: mean(Sigma_theta_phi_phi_hat) / Sigma_theta_phi_phi_true
      #     cov_theta[1,1] IS Sigma_theta (asymptotic variance of sqrt(n)*phi_hat),
      #     so no n factor is needed -- analogous to n_s*diag(vcovfunc) for beta.
      #     (scalar, same for all terms in this cell; NA for GLM)
      phi_var_A <- if (!is.null(cell_data$phi_var_vec) &&
                        is.finite(tv$vcov_phi_phi_true) &&
                        tv$vcov_phi_phi_true > 0)
        mean(cell_data$phi_var_vec, na.rm = TRUE) / tv$vcov_phi_phi_true
      else NA_real_

      # (10) phi variance check B: n_s * Var_empirical(phi_hat) / Sigma_theta_phi_phi_true
      #      Empirical counterpart: n_s * var(phi_vec) estimates Sigma_theta_{phi,phi}.
      phi_var_B <- if (!is.null(cell_data$phi_vec) &&
                        is.finite(tv$vcov_phi_phi_true) &&
                        tv$vcov_phi_phi_true > 0)
        n_s_eff * stats::var(cell_data$phi_vec, na.rm = TRUE) / tv$vcov_phi_phi_true
      else NA_real_

      # (11) phi-beta covariance check A: per-term ratio (NA for GLM or near-zero truth)
      phi_cov_A <- if (!is.null(cell_data$phi_cov_mat) &&
                        all(is.finite(tv$vcov_phi_beta_true)))
        colMeans(cell_data$phi_cov_mat, na.rm = TRUE) / tv$vcov_phi_beta_true
      else setNames(rep(NA_real_, length(coef_terms)), coef_terms)

      # (12) Raw MC variances of the estimator components (no Jacobian involved).
      #   var_R_MC_k   = n_s * var(R_hat_k)     -- actual MC variance of RESI
      #   var_b_MC_k   = n_s * var(beta_hat_k)  -- MC variance of numerator
      #   var_se_MC_k  = n_s * var(vcov_mat_k)  -- MC variance of scaled sandwich vcov
      # Compare var_R_MC to sigma2S_true (analytic total) to get the total ratio.
      # Compare var_b_MC to vcov_diag_true to get the direct-channel ratio (= vcov_B).
      # If sig2S_B >> vcov_B, the A/B chain or nonlinear terms are responsible.
      var_R_MC <- n_s_eff * apply(cell_data$R_mat,    2L, stats::var, na.rm = TRUE)
      var_b_MC <- n_s_eff * apply(cell_data$beta_mat, 2L, stats::var, na.rm = TRUE)
      var_se_MC <- n_s_eff * apply(cell_data$vcov_mat, 2L, stats::var, na.rm = TRUE)

      row <- data.frame(
        model        = s$type,
        vcov         = s$vcov_name,
        deriv_method = deriv_method,
        n            = n_s_eff,
        n_success    = cell_data$n_success,
        term         = coef_terms,
        phi_bias     = phi_bias,
        beta_bias    = beta_bias,
        vcov_A       = vcov_A,
        vcov_B       = vcov_B,
        R_bias       = R_bias,
        sig2S_A      = sig2S_A,
        sig2S_B      = sig2S_B,
        cv_Vhat      = cv_Vhat,
        phi_var_A    = phi_var_A,
        phi_var_B    = phi_var_B,
        phi_cov_A    = phi_cov_A,
        var_R_MC     = var_R_MC,
        var_b_MC     = var_b_MC,
        var_se_MC    = var_se_MC,
        row.names    = NULL,
        stringsAsFactors = FALSE
      )
      all_metrics[[length(all_metrics) + 1L]] <- row
    }
  }

  summary_df <- do.call(rbind, all_metrics)
  saveRDS(summary_df, file.path(output.dir, paste0("calibration_summary_",  deriv_method, ".rds")))
  saveRDS(true_vals,  file.path(output.dir, paste0("calibration_truevals_", deriv_method, ".rds")))
  message("Simulation complete (deriv_method = '", deriv_method, "'). Results saved to: ", output.dir)
  invisible(summary_df)
}

#' Calibration Figures for RESI Variance Estimates
#'
#' Reads output from \code{\link{simCalibrationSim}} and produces one PDF per
#' model setting showing convergence of point estimates, covariance estimates,
#' RESI estimates, and RESI variance estimates to their population targets.
#'
#' @param output.dir Character, directory containing output from
#'   \code{\link{simCalibrationSim}}. Default \code{"resiCalibrationSim"}.
#' @param figures.dir Character, output directory for PDFs.
#'   Default \code{file.path(output.dir, "figures")}.
#' @param alpha Numeric, nominal CI level. Default 0.05.
#'
#' @return Invisibly returns the summary data frame. Saves PDF figures to
#'   \code{figures.dir}.
#' @seealso \code{\link{simCalibrationSim}}, \code{\link{insurancePlasmodeSim}}
#' @importFrom splines ns
#' @importFrom sandwich vcovHC
#' @importFrom stats lm glm vcov binomial setNames
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot lines points abline legend par axis plot.new layout
#' @export
simCalibrationFigures <- function(
    output.dir   = "resiCalibrationSim",
    figures.dir  = NULL,
    alpha        = 0.05,
    deriv_method = c("corrected", "original", "population", "population2", "zeroB", "indep_moment")
) {
  deriv_method <- match.arg(deriv_method)
  if (is.null(figures.dir))
    figures.dir <- file.path(output.dir, "figures", deriv_method)
  dir.create(figures.dir, recursive = TRUE, showWarnings = FALSE)

  summary_path  <- file.path(output.dir, paste0("calibration_summary_",  deriv_method, ".rds"))
  truevals_path <- file.path(output.dir, paste0("calibration_truevals_", deriv_method, ".rds"))
  if (!file.exists(summary_path))
    stop("'", summary_path, "' not found. Run simCalibrationSim(deriv_method='", deriv_method, "') first.")
  if (!file.exists(truevals_path))
    stop("'", truevals_path, "' not found. Run simCalibrationSim(deriv_method='", deriv_method, "') first.")

  summary_df <- readRDS(summary_path)
  true_vals  <- readRDS(truevals_path)

  # Reconstruct minimal model_settings from true_vals names (label -> type/vcov_name)
  model_settings <- lapply(names(true_vals), function(lbl) {
    parts <- strsplit(lbl, "_", fixed = TRUE)[[1L]]
    list(type      = parts[1L],
         label     = lbl,
         vcov_name = parts[2L])
  })

  # ---- Figures ---------------------------------------------------------------
  # Abbreviate long term names for display in legends
  .abbrev <- function(x) {
    x <- sub("splines::ns\\(age, df = [0-9]+\\)", "ns(age)", x)
    x <- sub("regionnorthwest", "reg:NW", x)
    x <- sub("regionsoutheast", "reg:SE", x)
    x <- sub("regionsouthwest", "reg:SW", x)
    x <- sub("regionnortheast", "reg:NE", x)
    x <- sub(":sexmale$",       ":sex",   x)
    x
  }

  .sim_pal <- c("#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD",
                "#8C564B", "#E377C2", "#17BECF", "#BCBD22", "#7F7F7F",
                "#AEC7E8", "#98DF8A")
  n_vals_plot <- sort(unique(summary_df$n))

  for (s in model_settings) {
    lbl      <- s$label
    is_lm    <- (s$type == "lm")
    tv       <- true_vals[[lbl]]
    coef_terms <- tv$coef_terms
    n_terms    <- length(coef_terms)
    term_cols  <- setNames(.sim_pal[seq_len(n_terms)], coef_terms)

    sub <- summary_df[summary_df$model == s$type &
                        summary_df$vcov  == s$vcov_name, ]
    if (nrow(sub) == 0L) next

    fig_path <- file.path(figures.dir, paste0("calibration_", lbl, ".pdf"))

    # Layout: [bias_theta | vcov_ratios] / [bias_R | sig2S_ratios] / [legend]
    # Panel 1: Bias of theta (phi + beta_k lines; phi as a special dotted line)
    # Panel 2: vcov ratios  (solid=A, dashed=B, one color per term)
    # Panel 3: Bias of R_k
    # Panel 4: sigma2S ratios (solid=A, dashed=B, dotted=Th1)
    # Panel 5: legend (spans full width)
    is_glm <- !is_lm
    grDevices::pdf(fig_path, width = 9, height = 8.5)
    graphics::layout(
      matrix(c(1, 2, 3, 4, 5, 5), nrow = 3L, byrow = TRUE),
      widths  = c(4.0, 4.0),
      heights = c(3.0, 3.0, 2.0)
    )

    # --- Panel 1: Bias(theta) ---
    # Collect bias range for beta + phi
    all_bias_theta <- unlist(lapply(coef_terms, function(tm) {
      sub[sub$term == tm, "beta_bias"]
    }))
    if (is_lm) {
      phi_biases <- sub[sub$term == coef_terms[1L], "phi_bias"]
      all_bias_theta <- c(all_bias_theta, phi_biases)
    }
    ylim1 <- range(c(all_bias_theta, 0), na.rm = TRUE)
    ylim1 <- ylim1 + diff(ylim1) * c(-0.06, 0.06)
    if (is_glm) ylim1 <- c(-0.2, 0.2)   # clip glm bias panel

    graphics::par(mar = c(3.2, 3.5, 2.2, 0.5), mgp = c(2.0, 0.5, 0))
    graphics::plot(NULL, xlim = range(n_vals_plot), ylim = ylim1,
                   xlab = "n", ylab = "Bias",
                   main = paste0(lbl, ": Bias of theta"),
                   log = "x", xaxt = "n", bty = "l")
    graphics::axis(1L, at = n_vals_plot, labels = n_vals_plot, las = 2L,
                   mgp = c(1.7, 0.35, 0))
    graphics::abline(h = 0, lty = 2L, col = "gray40")

    for (tm in coef_terms) {
      d <- sub[sub$term == tm, ]
      d <- d[order(d$n), ]
      graphics::lines(d$n,  d$beta_bias, col = term_cols[[tm]], lwd = 1.5)
      graphics::points(d$n, d$beta_bias, col = term_cols[[tm]], pch = 16L, cex = 0.7)
    }
    if (is_lm) {
      d <- sub[sub$term == coef_terms[1L], ]
      d <- d[order(d$n), ]
      graphics::lines(d$n,  d$phi_bias, col = "#000000", lwd = 2L, lty = 3L)
      graphics::points(d$n, d$phi_bias, col = "#000000", pch = 17L, cex = 0.8)
    }

    # --- Panel 2: normalized-bias covariance calibration ---
    # For each n, load raw cell, build full theta matrix (phi + beta for lm),
    # and compute the normalized bias of every covariance element:
    #   diff[j,k] = (Sigma_analytic[j,k] - n_s*cov_MC[j,k]) /
    #               sqrt(n_s*var_MC[j,j] * n_s*var_MC[k,k])
    # = (Sigma_A[j,k] - Sigma_MC[j,k]) / sqrt(Sigma_MC[j,j] * Sigma_MC[k,k])
    # Diagonal: diff[k,k] = Sigma_A[k,k]/Sigma_MC[k,k] - 1 = 1/vcov_B_k - 1
    # Off-diagonal: normalized difference in covariance (target 0)
    # Includes phi row/col for lm; colored lines = diagonal, gray = off-diagonal.
    raw_dir_p   <- file.path(output.dir, "raw", deriv_method)
    Sig_A_mat   <- tv$vcov_full_true_mat    # analytic Sigma_theta (includes phi if lm)
    lm_has_phi  <- is_lm && !is.null(Sig_A_mat) && "phi" %in% rownames(Sig_A_mat)
    all_labels  <- if (!is.null(Sig_A_mat)) rownames(Sig_A_mat) else coef_terms
    n_all       <- length(all_labels)

    diag_diff_list <- vector("list", length(n_vals_plot))
    off_diff_list  <- vector("list", length(n_vals_plot))
    for (i_n in seq_along(n_vals_plot)) {
      nv <- n_vals_plot[i_n]
      rp <- file.path(raw_dir_p, paste0(lbl, "_n", nv, ".rds"))
      if (!file.exists(rp)) next
      cd <- readRDS(rp)
      if (is.null(cd$beta_mat) || nrow(cd$beta_mat) < 3L) next
      n_s_p <- cd$n
      bm    <- cd$beta_mat
      if (!all(coef_terms %in% colnames(bm))) next
      bm <- bm[, coef_terms, drop = FALSE]
      # Full theta matrix: prepend phi for lm
      theta_mc <- if (lm_has_phi && !is.null(cd$phi_vec))
        cbind(phi = cd$phi_vec, bm) else bm
      if (is.null(Sig_A_mat)) next
      Sig_MC <- n_s_p * stats::cov(theta_mc, use = "pairwise.complete.obs")
      # Normalized bias: (Sig_A - Sig_MC) / sqrt(Sig_MC_jj * Sig_MC_kk)
      sd_MC   <- sqrt(pmax(diag(Sig_MC), .Machine$double.eps))
      Sig_MC_nm <- Sig_A_mat  # same dimension; use A for name alignment
      if (nrow(Sig_MC) == nrow(Sig_A_mat)) {
        diff_mat <- (Sig_A_mat - Sig_MC) /
          outer(sd_MC, sd_MC)  # element-wise divide by product of MC SDs
      } else {
        diff_mat <- matrix(NA_real_, nrow(Sig_A_mat), ncol(Sig_A_mat))
      }
      rownames(diff_mat) <- colnames(diff_mat) <- all_labels
      diag_diff_list[[i_n]] <- diag(diff_mat)
      off_diff_list[[i_n]]  <- diff_mat[upper.tri(diff_mat)]
    }

    diag_diff_mat <- do.call(rbind, lapply(diag_diff_list, function(x)
      if (is.null(x)) rep(NA_real_, n_all) else x[all_labels]))
    off_diff_vals <- unlist(off_diff_list)
    all_diff_vals <- c(as.vector(diag_diff_mat), off_diff_vals)
    ylim2 <- range(c(all_diff_vals, 0), na.rm = TRUE)
    ylim2 <- ylim2 + diff(ylim2) * c(-0.06, 0.06)
    if (!is.finite(ylim2[1L])) ylim2 <- c(-1, 2)
    if (is_glm) ylim2[2L] <- min(ylim2[2L], 3)

    graphics::par(mar = c(3.2, 3.5, 2.2, 0.5), mgp = c(2.0, 0.5, 0))
    graphics::plot(NULL, xlim = range(n_vals_plot), ylim = ylim2,
                   xlab = "n",
                   ylab = expression((Sigma[A] - Sigma[MC]) / sqrt(Sigma[MC,jj] * Sigma[MC,kk])),
                   main = paste0(lbl, ": vcov norm. bias  (target = 0)"),
                   log = "x", xaxt = "n", bty = "l")
    graphics::axis(1L, at = n_vals_plot, labels = n_vals_plot, las = 2L,
                   mgp = c(1.7, 0.35, 0))
    graphics::abline(h = 0, lty = 2L, col = "gray40")

    # Off-diagonal: thin gray
    n_pairs <- if (n_all > 1L) ncol(utils::combn(n_all, 2L)) else 0L
    for (pair in seq_len(n_pairs)) {
      vals <- vapply(off_diff_list, function(x)
        if (is.null(x) || length(x) < pair) NA_real_ else x[pair], numeric(1L))
      if (any(is.finite(vals)))
        graphics::lines(n_vals_plot, vals, col = "gray70", lwd = 0.8)
    }

    # Diagonal: colored per term (beta), black for phi
    for (k_idx in seq_len(n_all)) {
      lbl_k <- all_labels[k_idx]
      vals  <- diag_diff_mat[, k_idx]
      if (!any(is.finite(vals))) next
      if (lbl_k == "phi") {
        graphics::lines(n_vals_plot,  vals, col = "#000000", lwd = 2L)
        graphics::points(n_vals_plot, vals, col = "#000000", pch = 15L, cex = 0.8)
      } else if (lbl_k %in% names(term_cols)) {
        graphics::lines(n_vals_plot,  vals, col = term_cols[[lbl_k]], lwd = 1.5)
        graphics::points(n_vals_plot, vals, col = term_cols[[lbl_k]], pch = 16L, cex = 0.7)
      }
    }

    # --- Panel 3: Bias(R) ---
    all_R_bias <- unlist(lapply(coef_terms, function(tm) sub[sub$term == tm, "R_bias"]))
    ylim3 <- range(c(all_R_bias, 0), na.rm = TRUE)
    ylim3 <- ylim3 + diff(ylim3) * c(-0.06, 0.06)
    if (is_glm) ylim3 <- c(-0.1, 0.1)   # clip glm bias panel

    graphics::par(mar = c(3.2, 3.5, 2.2, 0.5), mgp = c(2.0, 0.5, 0))
    graphics::plot(NULL, xlim = range(n_vals_plot), ylim = ylim3,
                   xlab = "n", ylab = "Bias",
                   main = paste0(lbl, ": Bias of R"),
                   log = "x", xaxt = "n", bty = "l")
    graphics::axis(1L, at = n_vals_plot, labels = n_vals_plot, las = 2L,
                   mgp = c(1.7, 0.35, 0))
    graphics::abline(h = 0, lty = 2L, col = "gray40")

    for (tm in coef_terms) {
      d <- sub[sub$term == tm, ]
      d <- d[order(d$n), ]
      graphics::lines(d$n,  d$R_bias, col = term_cols[[tm]], lwd = 1.5)
      graphics::points(d$n, d$R_bias, col = term_cols[[tm]], pch = 16L, cex = 0.7)
    }

    # --- Panel 4: sigma2S ratios (A solid, B dashed) ---
    all_sig2S <- unlist(lapply(coef_terms, function(tm) {
      d <- sub[sub$term == tm, ]
      c(d$sig2S_A, d$sig2S_B)
    }))
    ylim4 <- range(c(all_sig2S, 1), na.rm = TRUE)
    ylim4 <- ylim4 + diff(ylim4) * c(-0.06, 0.06)
    if (is_glm) ylim4 <- c(0.8, 3)   # clip glm ratio panel

    graphics::par(mar = c(3.2, 3.5, 2.2, 0.5), mgp = c(2.0, 0.5, 0))
    graphics::plot(NULL, xlim = range(n_vals_plot), ylim = ylim4,
                   xlab = "n",
                   ylab = "Ratio  (target = 1)",
                   main = paste(lbl, expression(sigma[S]^2), "calibration"),
                   log = "x", xaxt = "n", bty = "l")
    graphics::axis(1L, at = n_vals_plot, labels = n_vals_plot, las = 2L,
                   mgp = c(1.7, 0.35, 0))
    graphics::abline(h = 1, lty = 2L, col = "gray40")

    for (tm in coef_terms) {
      d    <- sub[sub$term == tm, ]
      d    <- d[order(d$n), ]
      # sig2S_Th1 = sig2S_B * sigma2S_true / sigma2S_true_Th1
      r_Th1 <- tv$sigma2S_true[tm] / tv$sigma2S_true_Th1[tm]
      sig2S_Th1_vals <- d$sig2S_B * r_Th1
      graphics::lines(d$n, d$sig2S_A,       col = term_cols[[tm]], lwd = 1.5, lty = 1L)
      graphics::lines(d$n, d$sig2S_B,       col = term_cols[[tm]], lwd = 1.5, lty = 2L)
      graphics::lines(d$n, sig2S_Th1_vals,  col = term_cols[[tm]], lwd = 1.5, lty = 3L)
      graphics::points(d$n, d$sig2S_A,      col = term_cols[[tm]], pch = 16L, cex = 0.7)
      graphics::points(d$n, d$sig2S_B,      col = term_cols[[tm]], pch = 1L,  cex = 0.7)
      graphics::points(d$n, sig2S_Th1_vals, col = term_cols[[tm]], pch = 2L,  cex = 0.7)
    }

    # --- Panel 5: Legend ---
    short_terms <- .abbrev(coef_terms)
    phi_leg <- if (is_lm) list(
      legend = c(short_terms, "phi (sigma^2)"),
      col    = c(unname(term_cols[coef_terms]), "#000000"),
      pch    = c(rep(16L, n_terms), 15L),
      lty    = c(rep(1L, n_terms), 1L)
    ) else list(
      legend = short_terms,
      col    = unname(term_cols[coef_terms]),
      pch    = rep(16L, n_terms),
      lty    = rep(1L, n_terms)
    )
    leg_ncol <- ceiling(length(phi_leg$legend) / 3L)

    graphics::par(mar = c(0.2, 0.5, 0.2, 0.5))
    graphics::plot.new()

    graphics::legend("left",
                     legend = phi_leg$legend,
                     col    = phi_leg$col,
                     pch    = phi_leg$pch,
                     lty    = phi_leg$lty,
                     lwd    = 1.5,
                     bty    = "n", cex = 0.80,
                     ncol   = leg_ncol,
                     title  = "Term", title.font = 2L)
    graphics::legend("right",
                     legend = c("Diag: vcov_B_k - 1  (variance excess)",
                                "Off-diag: cor_MC - cor_A  (gray)",
                                "sigma2S: Empirical / Th1 (C)"),
                     col    = c("black", "gray70", "gray30"),
                     lty    = c(1L, 1L, 3L),
                     pch    = c(16L, NA_integer_, 2L),
                     lwd    = c(1.5, 0.8, 1.5),
                     bty    = "n", cex = 0.85,
                     title  = "Line key", title.font = 2L)

    grDevices::dev.off()
    message("Saved: ", fig_path)
  }

  message("Calibration figures saved to: ", figures.dir)

  # ---- Structured return value -----------------------------------------------
  # For each setting: matrices of size (n_vals x terms) for each metric,
  # plus a phi_bias vector (length n_vals) for lm settings.
  structured <- lapply(model_settings, function(s) {
    lbl <- s$label
    sub <- summary_df[summary_df$model == s$type &
                        summary_df$vcov  == s$vcov_name, ]
    if (nrow(sub) == 0L) return(NULL)
    tv         <- true_vals[[lbl]]
    terms_here <- tv$coef_terms
    n_here     <- sort(unique(sub$n))

    # Helper: reshape one metric column into a (n_vals x terms) matrix
    .mat <- function(col) {
      m <- matrix(NA_real_, nrow = length(n_here), ncol = length(terms_here),
                  dimnames = list(n = as.character(n_here), term = terms_here))
      for (ni in seq_along(n_here)) {
        for (tm in terms_here) {
          v <- sub[sub$n == n_here[ni] & sub$term == tm, col]
          if (length(v) == 1L) m[ni, tm] <- v
        }
      }
      m
    }

    out <- list(
      n              = n_here,
      n_success      = setNames(
        vapply(n_here, function(ni) {
          v <- sub[sub$n == ni & sub$term == terms_here[1L], "n_success"]
          if (length(v)) v[1L] else NA_integer_
        }, integer(1L)), as.character(n_here)),
      # (1) Bias of raw coefficients
      beta_bias      = .mat("beta_bias"),
      # (2) vcov: analytical estimator ratio
      vcov_A         = .mat("vcov_A"),
      # (3) vcov: empirical-variance calibration ratio
      vcov_B         = .mat("vcov_B"),
      # (4) Bias of R
      R_bias         = .mat("R_bias"),
      # (5) sigma2S: analytical estimator ratio
      sig2S_A          = .mat("sig2S_A"),
      # (6) sigma2S: empirical-variance calibration ratio (key CI check)
      sig2S_B          = .mat("sig2S_B"),
      # (7) sigma2S: empirical variance / Th1 full-delta-method reference
      sig2S_Th1        = local({
        r <- tv$sigma2S_true[terms_here] / tv$sigma2S_true_Th1[terms_here]
        m <- .mat("sig2S_B")
        sweep(m, 2L, r, "*")
      }),
      # Reference values used for normalisation
      R_true              = tv$R_true,
      sigma2S_true        = tv$sigma2S_true,
      sigma2S_true_Th1    = tv$sigma2S_true_Th1,
      vcov_diag_true      = tv$vcov_diag_true[terms_here]
    )
    if (s$type == "lm") {
      phi_row <- sub[sub$term == terms_here[1L], ]
      phi_row <- phi_row[order(phi_row$n), ]
      out$phi_bias <- setNames(phi_row$phi_bias, as.character(phi_row$n))
      out$phi_true <- tv$phi_true
    }
    out
  })
  names(structured) <- sapply(model_settings, `[[`, "label")
  structured        <- Filter(Negate(is.null), structured)

  invisible(list(summary_table = summary_df, by_setting = structured))
}


#' Per-Term CI Method Comparison Figures
#'
#' Reads simulation summary tables from multiple output directories (one per CI
#' method) and produces PDF figures comparing CI methods across sample sizes.
#' For each (model type x variance estimator x table) combination two PDFs are
#' saved: an estimator comparison (Bias, MSE) and a CI comparison (SE
#' calibration, coverage, width).
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
#'   \code{\link{insurancePlasmodeSim}} call. Default \code{FALSE}.
#'
#' @return Invisibly returns the combined summary \code{data.frame}. Saves
#'   \code{estimator_<model>_<vcov>_<table>.pdf},
#'   \code{compare_<model>_<vcov>_<table>.pdf}, and
#'   \code{covquant_<model>_<vcov>_<table>.pdf} to \code{figures.dir}.
#' @seealso \code{\link{insurancePlasmodeSim}}, \code{\link{simFigures}},
#'   \code{\link{simEstimatorFigures}}
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot lines points abline legend par axis plot.new mtext
#' @importFrom splines ns
#' @importFrom sandwich vcovHC
#' @importFrom stats lm glm vcov binomial density setNames sd
#' @export
simCompareMethodsFigures <- function(
    output.dirs = c(Bootstrap   = "resiBootSim",
                    Normal = "resiAsympNormalSim",
                    QF    = "resiAsympQFSim",
                    CF     = "resiAsympCFSim"),
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
        # Drawing order per layout reading: SE -> UpperCov -> Width -> LowerCov
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
          # No n-based clipping: use all sample sizes for both lm and glm.
          cov_vals <- c(term_data[["upper_coverage"]], term_data[["lower_coverage"]])
          cov_ylim <- range(c(cov_vals, 1 - alpha/2 - 0.02, 1.01), na.rm = TRUE)

          # Width y-range: clip extreme values (1.5 for glm, 10 for lm)
          max_wid  <- if (mtype == "glm") 1.5 else 10
          wid_vals <- term_data[["width"]]
          wid_ylim <- range(wid_vals[wid_vals <= max_wid], na.rm = TRUE)

          # Drawing order: SE (spans) -> UpperCov -> Width (spans) -> LowerCov

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
          max_se  <- if (mtype == "glm") 3 else 10
          lim_se  <- range(se_vals[se_vals <= max_se], na.rm = TRUE)
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
  #  One PDF per (model x vcov x table): rows = terms, columns = n values
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
