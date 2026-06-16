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

  data.frame(
    bias     = colMeans(diff_mat,           na.rm = TRUE),
    mse      = colMeans(diff_mat ^ 2L,      na.rm = TRUE),
    coverage = colMeans(in_ci,              na.rm = TRUE),
    width    = colMeans(hi_mat - lo_mat,    na.rm = TRUE),
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
#'   Default 500. Use 10 for initial testing.
#' @param alpha Numeric, CI significance level. Default 0.05.
#' @param output.dir Character, path to the directory where all results are saved.
#'   Created if it does not exist. Default \code{"resiBootSim"}.
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
#' @seealso \code{\link{simFigures}}, \code{\link{resi}}, \code{\link{resi_pe}}
#' @importFrom parallel mclapply
#' @importFrom splines ns
#' @importFrom sandwich vcovHC
#' @importFrom stats lm glm vcov binomial sample
#' @export
insurancePlasmodeSim <- function(nsim              = 1000L,
                                  n.vec             = c(50, 100, 200, 500, 1000, 2000, 5000),
                                  nboot             = 500L,
                                  alpha             = 0.05,
                                  output.dir        = "resiBootSim",
                                  mc.cores.settings = 1L,
                                  mc.cores.reps     = 1L) {

  insurance <- RESI::insurance

  # Pre-compute expected factor levels for plasmode resampling check
  .fvars   <- c("sex", "smoker", "region")
  .flevels <- lapply(.fvars, function(v) unique(insurance[[v]]))
  names(.flevels) <- .fvars

  ci_lo <- paste0(alpha / 2 * 100, "%")
  ci_hi <- paste0((1 - alpha / 2) * 100, "%")

  lm_formula  <- log10(charges) ~ splines::ns(age, df = 3) * sex + bmi + smoker + region
  glm_formula <- I(charges > 15000) ~ splines::ns(age, df = 3) * sex + bmi + smoker + region

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
    resi_pe(full_mod, data = insurance, vcovfunc = s$vcovfunc)
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
             alpha = alpha, store.boot = FALSE),
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
#'   names. Default \code{"bootstrap"}.
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
                        ci.label   = "bootstrap") {

  # High-contrast matte palette (matplotlib tab10 + extensions)
  .sim_pal <- c(
    "#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD",
    "#8C564B", "#E377C2", "#17BECF", "#BCBD22", "#7F7F7F",
    "#AEC7E8", "#98DF8A"
  )

  summary_table <- readRDS(file.path(output.dir, "summary_table.rds"))
  figures_dir   <- file.path(output.dir, "figures")
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  n_vals        <- sort(unique(summary_table$n))
  metrics       <- c("bias", "mse", "coverage", "width")
  metric_labels <- c("Bias", "MSE", "CI Coverage", "CI Width")

  for (mtype in c("lm", "glm")) {
    for (vtype in c("parametric", "robust")) {

      sub <- summary_table[summary_table$model == mtype &
                             summary_table$vcov  == vtype, ]
      if (nrow(sub) == 0L) next

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

      # Layout: 2 rows × 5 cols; col 5 of each row is the legend panel.
      # Panels fill sequentially: 1-4 (anova metrics), 5 (anova legend),
      # 6-9 (coef metrics), 10 (coef legend).
      grDevices::pdf(fig_path, width = 15, height = 7)
      graphics::layout(
        matrix(seq_len(10), nrow = 2L, byrow = TRUE),
        widths = c(rep(3, 4), 2.2)
      )

      for (ri in row_info) {
        sub_tbl <- sub[sub$table == ri$table, ]

        # Legend cex: as large as possible while fitting all terms
        cex_leg <- min(1.0, 9 / length(ri$terms))

        for (m_idx in seq_along(metrics)) {
          metric <- metrics[m_idx]
          mlab   <- metric_labels[m_idx]

          y_all  <- sub_tbl[[metric]]
          ylim   <- range(y_all, na.rm = TRUE)
          if (metric == "coverage") {
            ylim <- range(c(ylim, 1 - alpha - 0.05, 1 + 0.02))
          }
          if (metric == "bias") {
            ylim <- range(c(ylim, 0))  # ensure zero is visible
          }

          graphics::par(mar = c(2.8, 2.8, 1.8, 0.4), mgp = c(1.7, 0.45, 0))
          graphics::plot(
            NULL,
            xlim = range(n_vals), ylim = ylim,
            xlab = "Sample Size",
            ylab = mlab,
            main = paste(toupper(mtype), ri$main_prefix, mlab,
                         paste0("(", vtype, ")")),
            log  = "x",
            xaxt = "n",
            bty  = "l"
          )
          graphics::axis(1, at = n_vals, labels = n_vals, las = 2L,
                         mgp = c(1.7, 0.35, 0))

          if (metric == "coverage") {
            graphics::abline(h = 1 - alpha, lty = 2L, col = "gray40")
          }
          if (metric == "bias") {
            graphics::abline(h = 0, lty = 2L, col = "gray40")
          }

          for (term in ri$terms) {
            td <- sub_tbl[sub_tbl$term == term, ]
            td <- td[order(td$n), ]
            graphics::lines(td$n,  td[[metric]], col = ri$cols[term], lwd = 2)
            graphics::points(td$n, td[[metric]], col = ri$cols[term], pch = 16,
                             cex = 0.8)
          }
        }

        # Legend panel — minimal margins, maximise font size
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
      }

      grDevices::dev.off()
      message("Saved: ", fig_path)
    }
  }

  invisible(summary_table)
}
