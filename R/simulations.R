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
#  simBiasWidthFigures
# ============================================================

#' Bias / CI Width Ratio Figures Across CI Methods
#'
#' Reads simulation summary tables from multiple output directories (one per CI
#' method) and produces one PDF per (model type \eqn{\times} variance estimator)
#' combination. Each figure has one row per CI method showing \code{bias / CI width}
#' vs sample size for anova and coefficient tables. Values near \eqn{\pm 0.5}
#' indicate that bias alone can shift the CI off the true value.
#'
#' @param output.dirs Named character vector mapping CI method labels to their
#'   simulation output directories. Default:
#'   \code{c(boot = "resiBootSim", normal = "resiAsympNormalSim", qf = "resiAsympQFSim")}.
#'   Directories that do not exist are silently skipped.
#' @param figures.dir Character, directory where figures are saved. Default:
#'   \code{file.path(output.dirs[[1]], "figures")}.
#' @param alpha Numeric, unused; kept for API consistency. Default 0.05.
#'
#' @return Invisibly returns the combined summary \code{data.frame}. Saves PDF
#'   figures named \code{sim_biaswidth_<model>_<vcov>.pdf} to \code{figures.dir}.
#' @seealso \code{\link{simFigures}}, \code{\link{simCompareMethodsFigures}}
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot lines points abline legend par axis plot.new layout
#' @export
simBiasWidthFigures <- function(
    output.dirs = c(boot   = "resiBootSim",
                    normal = "resiAsympNormalSim",
                    qf     = "resiAsympQFSim"),
    figures.dir = NULL,
    alpha       = 0.05) {

  method_names <- names(output.dirs)
  if (is.null(method_names) || any(method_names == ""))
    stop("output.dirs must be a named vector (names = CI method labels)")

  if (is.null(figures.dir))
    figures.dir <- file.path(output.dirs[[1]], "figures")
  dir.create(figures.dir, recursive = TRUE, showWarnings = FALSE)

  .sim_pal <- c(
    "#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD",
    "#8C564B", "#E377C2", "#17BECF", "#BCBD22", "#7F7F7F",
    "#AEC7E8", "#98DF8A"
  )

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
  combined$bias_width_ratio <- combined$bias / combined$width

  avail_methods <- intersect(method_names, unique(combined$ci_method))
  n_methods     <- length(avail_methods)
  n_vals        <- sort(unique(combined$n))

  for (mtype in c("lm", "glm")) {
    for (vtype in c("parametric", "robust")) {
      sub <- combined[combined$model == mtype & combined$vcov == vtype, ]
      if (nrow(sub) == 0L) next

      anova_terms <- unique(sub$term[sub$table == "anova"])
      coef_terms  <- unique(sub$term[sub$table == "coefficients"])
      anova_cols  <- setNames(.sim_pal[seq_along(anova_terms)], anova_terms)
      coef_cols   <- setNames(.sim_pal[seq_along(coef_terms)],  coef_terms)
      cex_anova   <- min(1.0, 9 / length(anova_terms))
      cex_coef    <- min(1.0, 9 / length(coef_terms))

      fig_path <- file.path(figures.dir,
        paste0("sim_biaswidth_", mtype, "_", vtype, ".pdf"))

      # Layout: n_methods rows x 4 cols [anova_plot, coef_plot, leg_anova, leg_coef]
      # Cols 3-4 (legends) span all rows via repeated panel index.
      leg_anova_idx <- n_methods * 2L + 1L
      leg_coef_idx  <- n_methods * 2L + 2L
      layout_mat <- matrix(0L, nrow = n_methods, ncol = 4L)
      for (i in seq_len(n_methods)) {
        layout_mat[i, 1L] <- (i - 1L) * 2L + 1L
        layout_mat[i, 2L] <- (i - 1L) * 2L + 2L
        layout_mat[i, 3L] <- leg_anova_idx
        layout_mat[i, 4L] <- leg_coef_idx
      }

      grDevices::pdf(fig_path, width = 12, height = 3.5 * n_methods)
      graphics::layout(layout_mat,
                       widths  = c(4, 4, 1.8, 1.8),
                       heights = rep(1, n_methods))

      for (mi in seq_along(avail_methods)) {
        meth  <- avail_methods[[mi]]
        sub_m <- sub[sub$ci_method == meth, ]

        for (tbl_info in list(
          list(table = "anova",        terms = anova_terms, cols = anova_cols,
               prefix = "Anova"),
          list(table = "coefficients", terms = coef_terms,  cols = coef_cols,
               prefix = "Coef")
        )) {
          sub_tbl <- sub_m[sub_m$table == tbl_info$table, ]

          yr   <- range(sub_tbl$bias_width_ratio, na.rm = TRUE, finite = TRUE)
          ylim <- range(c(yr, 0), na.rm = TRUE)

          graphics::par(mar = c(3.2, 3.2, 2.0, 0.5), mgp = c(1.9, 0.5, 0))
          graphics::plot(NULL, xlim = range(n_vals), ylim = ylim,
                         xlab = "Sample Size", ylab = "Bias / CI Width",
                         main = paste(toupper(mtype), tbl_info$prefix, "Bias/Width",
                                      paste0("(", vtype, ", ", meth, ")")),
                         log = "x", xaxt = "n", bty = "l")
          graphics::axis(1, at = n_vals, labels = n_vals, las = 2L,
                         mgp = c(1.7, 0.35, 0))
          graphics::abline(h = 0,             lty = 2L, col = "gray40")
          graphics::abline(h = c(-0.5, 0.5),  lty = 3L, col = "gray70")

          for (term in tbl_info$terms) {
            td <- sub_tbl[sub_tbl$term == term, ]
            td <- td[order(td$n), ]
            if (nrow(td) == 0L) next
            graphics::lines(td$n,  td$bias_width_ratio,
                            col = tbl_info$cols[term], lwd = 2L)
            graphics::points(td$n, td$bias_width_ratio,
                             col = tbl_info$cols[term], pch = 16L, cex = 0.8)
          }
        }
      }

      # Shared legends (drawn last; span all rows via layout)
      graphics::par(mar = c(0.5, 0.3, 0.5, 0.3))
      graphics::plot.new()
      graphics::legend("center", legend = anova_terms,
                       col = anova_cols[anova_terms],
                       lwd = 2L, pch = 16L, bty = "n", cex = cex_anova,
                       title = "Anova terms", title.font = 2L)

      graphics::par(mar = c(0.5, 0.3, 0.5, 0.3))
      graphics::plot.new()
      graphics::legend("center", legend = coef_terms,
                       col = coef_cols[coef_terms],
                       lwd = 2L, pch = 16L, bty = "n", cex = cex_coef,
                       title = "Coef terms", title.font = 2L)

      grDevices::dev.off()
      message("Saved: ", fig_path)
    }
  }

  invisible(combined)
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
#'   \code{c(boot = "resiBootSim", normal = "resiAsympNormalSim", qf = "resiAsympQFSim")}.
#'   Directories that do not exist are silently skipped.
#' @param figures.dir Character, directory where comparison figures are saved.
#'   Default: \code{file.path(output.dirs[1], "figures", "method_comparison")}.
#' @param alpha Numeric, nominal CI level used for coverage reference lines.
#'   Default 0.05.
#'
#' @return Invisibly returns the combined summary \code{data.frame}. Saves one PDF
#'   per (model \eqn{\times} vcov \eqn{\times} table \eqn{\times} term) to
#'   \code{figures.dir}, named
#'   \code{compare_<model>_<vcov>_<table>_<term>.pdf}.
#' @seealso \code{\link{insurancePlasmodeSim}}, \code{\link{simFigures}}
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot lines points abline legend par axis plot.new
#' @export
simCompareMethodsFigures <- function(
    output.dirs = c(boot   = "resiBootSim",
                    normal = "resiAsympNormalSim",
                    qf     = "resiAsympQFSim"),
    figures.dir = NULL,
    alpha       = 0.05) {

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

  # Panel order per row: Bias, MSE, SE comparison, Upper Cov, Lower Cov, Width
  metrics       <- c("bias", "mse", "upper_coverage", "lower_coverage", "width")
  metric_labels <- c("Bias", "MSE", "Upper Cov.", "Lower Cov.", "CI Width")
  metric_refs   <- list(
    bias           = list(h = 0,           col = "gray40"),
    mse            = NULL,
    upper_coverage = list(h = 1 - alpha/2, col = "gray40"),
    lower_coverage = list(h = 1 - alpha/2, col = "gray40"),
    width          = NULL
  )
  # For GLM, truncate non-coverage y-axes to n >= 500
  is_coverage <- c(bias = FALSE, mse = FALSE, upper_coverage = TRUE,
                   lower_coverage = TRUE, width = FALSE)

  n_data_panels <- 6L   # Bias, MSE, SE comparison, Upper Cov, Lower Cov, Width

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

        # One PDF per (model x vcov x table), rows = terms
        fig_path <- file.path(
          figures.dir,
          paste0("compare_", mtype, "_", vtype, "_", tbl_name, ".pdf")
        )

        # Layout: n_terms rows x (n_data_panels + 1) cols;
        # last column is a single legend panel spanning all rows
        legend_id <- n_data_panels * n_terms + 1L
        lay_mat   <- matrix(0L, nrow = n_terms, ncol = n_data_panels + 1L)
        for (ti in seq_len(n_terms))
          lay_mat[ti, seq_len(n_data_panels)] <-
            (ti - 1L) * n_data_panels + seq_len(n_data_panels)
        lay_mat[, n_data_panels + 1L] <- legend_id   # shared legend column

        row_h <- 3.2
        grDevices::pdf(fig_path, width = 22, height = row_h * n_terms + 0.3)
        graphics::layout(
          lay_mat,
          widths  = c(rep(2.8, n_data_panels), 2.0),
          heights = rep(row_h, n_terms)
        )

        for (ti in seq_along(terms)) {
          term      <- terms[[ti]]
          term_data <- sub[sub$term == term, ]

          # y-range helper: for GLM non-coverage metrics, use only n >= 500;
          # for SE and width, ignore values > 10 when computing limits
          ylim_for <- function(m) {
            rows <- if (mtype == "glm" && !is_coverage[[m]])
                      term_data$n >= 500 else rep(TRUE, nrow(term_data))
            vals <- term_data[rows, m]
            if (m == "width") vals <- vals[vals <= 10]
            range(vals, na.rm = TRUE)
          }

          # ---- Panel 1: Bias (carries row title = term name) ----
          m   <- "bias"
          ref <- metric_refs[[m]]
          ylim <- ylim_for(m)
          if (!is.null(ref)) ylim <- range(c(ylim, ref$h - 0.02, ref$h + 0.01))
          graphics::par(mar = c(3.2, 3.0, 2.2, 0.5), mgp = c(1.8, 0.5, 0))
          graphics::plot(NULL, xlim = range(n_vals), ylim = ylim,
                         xlab = "Sample Size", ylab = "Bias",
                         main = term, log = "x", xaxt = "n", bty = "l")
          graphics::axis(1, at = n_vals, labels = n_vals, las = 2L,
                         mgp = c(1.7, 0.35, 0))
          graphics::abline(h = 0, lty = 2L, col = "gray40")
          for (meth in avail_methods) {
            md <- term_data[term_data$ci_method == meth, ]
            md <- md[order(md$n), ]
            if (nrow(md) == 0L) next
            graphics::lines(md$n,  md$bias, col = method_cols[[meth]], lwd = 2L)
            graphics::points(md$n, md$bias, col = method_cols[[meth]],
                             pch = 16L, cex = 0.8)
          }

          # ---- Panel 2: MSE ----
          m   <- "mse"
          ylim <- ylim_for(m)
          graphics::par(mar = c(3.2, 3.0, 2.2, 0.5), mgp = c(1.8, 0.5, 0))
          graphics::plot(NULL, xlim = range(n_vals), ylim = ylim,
                         xlab = "Sample Size", ylab = "MSE",
                         main = "", log = "x", xaxt = "n", bty = "l")
          graphics::axis(1, at = n_vals, labels = n_vals, las = 2L,
                         mgp = c(1.7, 0.35, 0))
          for (meth in avail_methods) {
            md <- term_data[term_data$ci_method == meth, ]
            md <- md[order(md$n), ]
            if (nrow(md) == 0L) next
            graphics::lines(md$n,  md$mse, col = method_cols[[meth]], lwd = 2L)
            graphics::points(md$n, md$mse, col = method_cols[[meth]],
                             pch = 16L, cex = 0.8)
          }

          # ---- Panel 3: SE comparison (method-colored arrows, sequential n) ----
          se_comp <- do.call(rbind, lapply(avail_methods, function(meth) {
            md <- term_data[term_data$ci_method == meth, ]
            md <- md[order(md$n), ]
            data.frame(
              ci_method    = meth,
              n            = md$n,
              empirical_se = sqrt(md$n) * md$empirical_sd,
              estimated_se = sqrt(md$n) * md$width / (2 * z_ref),
              stringsAsFactors = FALSE
            )
          }))
          se_vals <- c(se_comp$empirical_se, se_comp$estimated_se)
          lim_se  <- range(se_vals[se_vals <= 10], na.rm = TRUE)
          se_title <- if (ti == 1L) "SE Calibration" else ""
          graphics::par(mar = c(3.2, 3.0, 2.2, 0.5), mgp = c(1.8, 0.5, 0))
          graphics::plot(
            NULL, xlim = lim_se, ylim = lim_se,
            xlab = expression(sqrt(n) %*% " Empirical SE"),
            ylab = expression(sqrt(n) %*% " Estimated SE"),
            main = se_title, bty = "l", asp = 1
          )
          graphics::abline(0, 1, lty = 2L, col = "gray40")
          for (meth in avail_methods) {
            md_se <- se_comp[se_comp$ci_method == meth, ]
            md_se <- md_se[order(md_se$n), ]
            if (nrow(md_se) == 0L) next
            # Starting point marker
            graphics::points(md_se$empirical_se[1L], md_se$estimated_se[1L],
                             col = method_cols[[meth]], pch = 16L,
                             cex = 1.1)
            # Arrows connecting sequential n values
            if (nrow(md_se) >= 2L) {
              for (i in seq_len(nrow(md_se) - 1L)) {
                graphics::arrows(
                  md_se$empirical_se[i],   md_se$estimated_se[i],
                  md_se$empirical_se[i+1L], md_se$estimated_se[i+1L],
                  col = method_cols[[meth]], lwd = 1.5,
                  length = 0.07, angle = 20L
                )
              }
            }
          }

          # ---- Panels 4-6: Upper Cov, Lower Cov, Width ----
          for (mi in 3:5) {
            m    <- metrics[[mi]]
            ref  <- metric_refs[[m]]
            ylim <- ylim_for(m)
            if (!is.null(ref)) ylim <- range(c(ylim, ref$h - 0.02, ref$h + 0.01))
            graphics::par(mar = c(3.2, 3.0, 2.2, 0.5), mgp = c(1.8, 0.5, 0))
            graphics::plot(NULL, xlim = range(n_vals), ylim = ylim,
                           xlab = "Sample Size", ylab = metric_labels[[mi]],
                           main = "", log = "x", xaxt = "n", bty = "l")
            graphics::axis(1, at = n_vals, labels = n_vals, las = 2L,
                           mgp = c(1.7, 0.35, 0))
            if (!is.null(ref))
              graphics::abline(h = ref$h, lty = 2L, col = ref$col)
            for (meth in avail_methods) {
              md <- term_data[term_data$ci_method == meth, ]
              md <- md[order(md$n), ]
              if (nrow(md) == 0L) next
              graphics::lines(md$n,  md[[m]], col = method_cols[[meth]], lwd = 2L)
              graphics::points(md$n, md[[m]], col = method_cols[[meth]],
                               pch = 16L, cex = 0.8)
            }
          }
        }  # end term loop

        # Shared legend panel (drawn once after all rows, spans full column height)
        graphics::par(mar = c(0.3, 0.3, 0.3, 0.3))
        graphics::plot.new()
        graphics::legend(
          "center",
          legend = avail_methods,
          col    = method_cols[avail_methods],
          pch    = 16L,
          lwd    = 2L, lty = 1L,
          bty = "n", cex = 0.95,
          title = "CI method", title.font = 2L
        )

        grDevices::dev.off()
        message("Saved: ", fig_path)
      }
    }
  }

  message("Comparison figures saved to: ", figures.dir)
  invisible(combined)
}
