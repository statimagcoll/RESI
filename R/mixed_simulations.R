# [AI-generated: GitHub Copilot, 2026-11-15]
# [Prompt summary: Lyu-Welsh limiting RESI targets and hybrid plasmode
#  calibration of mixed-model RESI on mlmRev::egsingle.]

# Mixed-model RESI calibration: Lyu & Welsh (2022) limiting targets and the
# hybrid plasmode simulation for a three-level lmer model with a school random
# intercept and child-level random effects, e.g. `y ~ X + (1 | school/child)`
# or `y ~ X + (1 | school) + (year | child:school)` (developed on
# mlmRev::egsingle).
#
# Design convention (matches the plasmode): schools are drawn uniformly from
# the observed school profiles and students within a school in proportion to
# the observed school size, so a design with average m students per school
# gives school s a weight m * n_s / mean(n_s).
#
# RESI targets are beta_k / sqrt(n_k Var_kk) (coefficients) and
# sqrt(beta_T' (n_T Var_TT)^{-1} beta_T) (ANOVA terms), with n_k = g schools
# for between-school terms (constant within school, incl. the intercept) and
# n_k = N students for within-school terms, matching resi_pe(level = 2).

.mixedVarComps <- function(fit, school, child) {
  vc <- lme4::VarCorr(fit)
  if (length(vc[[school]]) != 1L) {
    stop("Only a random intercept is supported at the school level (", school, ").")
  }
  G_child <- as.matrix(vc[[child]])
  attributes(G_child) <- list(dim = dim(G_child), dimnames = dimnames(G_child))
  list(tau2_school = as.numeric(vc[[school]]),
       G_child = G_child, tau2_child = G_child[1L, 1L],
       sigma2 = stats::sigma(fit)^2)
}

# Child-level random-effects design Z (observations x q) of `fit`.
.mixedChildZ <- function(fit, child) {
  mm <- lme4::getME(fit, "mmList")
  hit <- which(sub(".*\\|\\s*", "", names(mm)) == child)
  if (length(hit) != 1L) stop("Expected one random-effects term for ", child, ".")
  as.matrix(mm[[hit]])
}

.mixedBetweenColumns <- function(X, school_id) {
  vapply(seq_len(ncol(X)), function(j) {
    all(tapply(X[, j], school_id, function(x) diff(range(x)) < 1e-10))
  }, logical(1L))
}

.mixedFixedTermLabels <- function(fit) {
  attr(stats::terms(reformulas::nobars(stats::formula(fit))), "term.labels")
}

# Per-school sufficient statistics of the GLS information for a three-level
# model with a school random intercept integrated out:
#   A_s = sum_j X_j' W_j^{-1} X_j,  b_s = sum_j X_j' W_j^{-1} 1,
#   c_s = sum_j 1' W_j^{-1} 1,      W_j = sigma2 I + Z_j G Z_j',
# so school s contributes A_s - b_s b_s' / (1 / tau2_school + c_s). Z holds the
# child-level random-effects columns (default: a random intercept, G = tau2)
# and W_j^{-1} = (I - Z_j M_j Z_j') / sigma2 with M_j = (sigma2 I + G Z_j'Z_j)^{-1} G.
lyuSchoolStats <- function(X, school_id, child_id, G_child, sigma2,
                           Z = matrix(1, nrow(X), 1L)) {
  G_child <- as.matrix(G_child)
  child_id <- droplevels(as.factor(child_id))
  school_of_child <- tapply(as.character(school_id), child_id, `[`, 1L)
  rows <- split(seq_len(nrow(X)), child_id)
  q <- ncol(Z)
  per_child <- lapply(rows, function(r) {
    Xj <- X[r, , drop = FALSE]
    Zj <- Z[r, , drop = FALSE]
    M <- solve(sigma2 * diag(q) + G_child %*% crossprod(Zj), G_child)
    XZ <- crossprod(Xj, Zj)
    z1 <- colSums(Zj)
    list(A = (crossprod(Xj) - XZ %*% M %*% t(XZ)) / sigma2,
         b = (colSums(Xj) - drop(XZ %*% M %*% z1)) / sigma2,
         c = (length(r) - drop(z1 %*% M %*% z1)) / sigma2)
  })
  schools <- unique(as.character(school_id))
  stats <- lapply(schools, function(s) {
    js <- which(school_of_child == s)
    list(A = Reduce(`+`, lapply(per_child[js], `[[`, "A")),
         b = Reduce(`+`, lapply(per_child[js], `[[`, "b")),
         c = sum(vapply(per_child[js], `[[`, numeric(1L), "c")),
         n_child = length(js))
  })
  names(stats) <- schools
  stats
}

.lyuResiTables <- function(beta, V, n, assign, term_labels) {
  coef_resi <- beta / sqrt(n * diag(V))
  anova_resi <- vapply(seq_along(term_labels), function(t) {
    cols <- which(assign == t)
    if (length(unique(n[cols])) != 1L) stop("Term spans between and within blocks: ", term_labels[t])
    Vt <- n[cols[1]] * V[cols, cols, drop = FALSE]
    sqrt(drop(crossprod(beta[cols], solve(Vt, beta[cols]))))
  }, numeric(1L))
  list(coefficients = data.frame(RESI = coef_resi, row.names = names(beta)),
       anova = data.frame(RESI = anova_resi, row.names = term_labels))
}

#' Lyu-Welsh limiting RESI for a three-level lmer fit
#'
#' Evaluates the g, m -> Inf limit (Lyu & Welsh, 2022) of the per-unit
#' covariance of the fixed effects and the corresponding RESI at the fitted
#' parameters. Between-school columns: n_g Var(beta_b) -> tau2_school C^{-1},
#' C = g^{-1} sum_s x_s x_s'. Within-school columns: the school intercept acts
#' as a fixed effect, so I_w / N = N^{-1} sum_s (A_s - b_s b_s' / c_s). The
#' between and within blocks are asymptotically independent.
#' @param fit A `lmerMod` fit with a school random intercept and child-level
#'   random effects (intercept and optionally slopes), e.g.
#'   `y ~ X + (1 | school/child)`.
#' @param school,child Names of the school and child grouping factors in
#'   `lme4::getME(fit, "flist")`.
#' @return List with `coefficients` and `anova` RESI tables, the logical
#'   `between` indicator per column, and the per-unit limiting covariance `V`.
#' @keywords internal
lyuLimitResi <- function(fit, school = "schoolid", child = "childid:schoolid") {
  X <- lme4::getME(fit, "X")
  fl <- lme4::getME(fit, "flist")
  vc <- .mixedVarComps(fit, school, child)
  beta <- lme4::fixef(fit)
  school_id <- fl[[school]]
  between <- .mixedBetweenColumns(X, school_id)
  Xb <- X[!duplicated(school_id), between, drop = FALSE]
  Vb <- vc$tau2_school * solve(crossprod(Xb) / nrow(Xb))

  st <- lyuSchoolStats(X[, !between, drop = FALSE], school_id, fl[[child]],
                       vc$G_child, vc$sigma2, .mixedChildZ(fit, child))
  N <- sum(vapply(st, `[[`, numeric(1L), "n_child"))
  Iw <- Reduce(`+`, lapply(st, function(s) s$A - tcrossprod(s$b) / s$c))
  Vw <- solve(Iw / N)

  V <- matrix(0, ncol(X), ncol(X), dimnames = list(colnames(X), colnames(X)))
  V[between, between] <- Vb
  V[!between, !between] <- Vw
  out <- .lyuResiTables(beta, V, rep(1, ncol(X)), attr(X, "assign"),
                        .mixedFixedTermLabels(fit))
  out$between <- setNames(between, colnames(X))
  out$V <- V
  out
}

# Information contributions I_s(m) of each observed school profile at
# average m students per school.
.mixedSchoolInfo <- function(fit, m, school, child) {
  X <- lme4::getME(fit, "X")
  fl <- lme4::getME(fit, "flist")
  vc <- .mixedVarComps(fit, school, child)
  school_id <- fl[[school]]
  st <- lyuSchoolStats(X, school_id, fl[[child]], vc$G_child, vc$sigma2, .mixedChildZ(fit, child))
  n_s <- vapply(st, `[[`, numeric(1L), "n_child")
  r <- m / mean(n_s)
  list(I_s = lapply(st, function(s) {
    r * s$A - (r * s$b) %*% t(r * s$b) / (1 / vc$tau2_school + r * s$c)
  }), between = .mixedBetweenColumns(X, school_id), X = X)
}

#' Fixed-m RESI under the expected hybrid plasmode design
#'
#' Exact RESI and per-unit SE when each school profile is represented in
#' proportion to its sampling probability with average `m` students per
#' school (the number of schools g cancels). Converges to [lyuLimitResi()] as
#' m -> Inf; equals `resi_pe(level = 2)` on the observed design at
#' m = mean school size.
#' @inheritParams lyuLimitResi
#' @param m Average number of students per school.
#' @return As [lyuLimitResi()] plus `se_unit`, sqrt(n_k Var_kk).
#' @keywords internal
expectedDesignResi <- function(fit, m, school = "schoolid", child = "childid:schoolid") {
  info <- .mixedSchoolInfo(fit, m, school, child)
  g <- length(info$I_s)
  V <- solve(Reduce(`+`, info$I_s))
  n <- ifelse(info$between, g, m * g)
  beta <- lme4::fixef(fit)
  out <- .lyuResiTables(beta, V, n, attr(info$X, "assign"), .mixedFixedTermLabels(fit))
  out$se_unit <- setNames(sqrt(n * diag(V)), names(beta))
  out
}

#' Design-averaged per-unit SE for the hybrid plasmode
#'
#' sqrt(n_k E_X[V(X)]_kk) at g schools and average m students per school,
#' averaging the conditional covariance V(X) = I(X)^{-1} over `n_rep`
#' resampled sets of g school profiles. Because beta-hat is conditionally
#' unbiased this is the finite-(g, m) reference for the Monte Carlo SD of
#' beta-hat; it exceeds `expectedDesignResi()$se_unit` since
#' E[I^{-1}] >= E[I]^{-1}.
#' @inheritParams expectedDesignResi
#' @param g Number of schools.
#' @param n_rep Number of resampled school sets.
#' @keywords internal
designAveragedSe <- function(fit, g, m, n_rep = 2000L,
                             school = "schoolid", child = "childid:schoolid") {
  info <- .mixedSchoolInfo(fit, m, school, child)
  v <- replicate(n_rep, tryCatch(
    diag(solve(Reduce(`+`, info$I_s[sample.int(length(info$I_s), g, TRUE)]))),
    error = function(e) rep(NA_real_, ncol(info$X))))
  n <- ifelse(info$between, g, m * g)
  setNames(sqrt(n * rowMeans(v, na.rm = TRUE)), colnames(info$X))
}

.mixedSimParams <- function(fit) {
  fe <- lme4::fixef(fit)
  names(fe) <- paste0("beta_", names(fe))
  vc <- as.data.frame(lme4::VarCorr(fit))
  vc_names <- ifelse(vc$grp == "Residual", "var_Residual",
                     ifelse(is.na(vc$var2), paste0("var_", vc$grp, "_", vc$var1),
                            paste0("cov_", vc$grp, "_", vc$var1, "_", vc$var2)))
  c(fe, setNames(vc$vcov, vc_names))
}

.hybridResiSummary <- function(reps, truth_resi, settings) {
  sources <- list(coefficients = "resi_coef", anova = "resi_anova")
  do.call(rbind, unlist(lapply(names(sources), function(table_name) {
    rep_name <- sources[[table_name]]
    truth <- truth_resi[[table_name]]
    lapply(rownames(truth), function(term) {
      values <- vapply(reps, function(rep) rep[[rep_name]][term, "RESI"], numeric(1L))
      replicate_n <- vapply(reps, function(rep) rep[[rep_name]][term, "n"], numeric(1L))
      first <- reps[[1L]][[rep_name]][term, , drop = FALSE]
      data.frame(
        as.list(settings), table = table_name, term = term,
        n = mean(replicate_n), n_sd = stats::sd(replicate_n),
        variation_level = first$`Variation level`,
        n_level = first$`n level`,
        truth_resi = truth[term, "RESI"],
        mc_mean_resi = mean(values, na.rm = TRUE),
        bias_resi = truth[term, "RESI"] - mean(values, na.rm = TRUE),
        mc_sd_resi = stats::sd(values, na.rm = TRUE),
        mc_se_resi = stats::sd(values, na.rm = TRUE) / sqrt(sum(is.finite(values))),
        n_success = sum(is.finite(values)),
        check.names = FALSE, row.names = NULL)
    })
  }), recursive = FALSE))
}

#' Hybrid plasmode calibration simulation for three-level mixed-model RESI
#'
#' Fits `formula` (ML) to `data` and treats the fitted fixed effects and
#' variance components as the data-generating truth. Each replicate
#' (i) resamples `n_schools` school covariate profiles with replacement and,
#' within each, students with replacement in proportion to the observed school
#' size so that the average is `n_per_site`; (ii) draws school and student
#' random effects with replacement from the BLUP pools, each centred and
#' rescaled to the fitted variance component, plus Gaussian residuals;
#' (iii) refits by REML and computes `resi_pe(level = 2)` with the model-based
#' covariance. RESI and per-unit SE targets are the Lyu-Welsh g, m -> Inf
#' limit at the full-data ML estimates ([lyuLimitResi()]).
#'
#' @param formula Model formula with a school random intercept and child-level
#'   random effects, e.g. `y ~ X + (1 | school/child)` or
#'   `y ~ X + (1 | school) + (year | child:school)`.
#' @param data Data frame containing the variables in `formula`.
#' @param school,child Names of the school and child id columns of `data`.
#' @param n_schools,n_per_site Vectors crossed to form the simulation grid.
#' @param nsim Number of replicates per grid cell.
#' @param mc.cores Number of cores for [parallel::mclapply()] over replicates.
#' @param seed Seed (with the "L'Ecuyer-CMRG" generator) for reproducibility.
#' @param control `lme4::lmerControl()` used for all fits.
#' @param rep_fun Optional `function(fit, design)` evaluated on each replicate
#'   fit; its return value is stored as `extra` in the replicate output.
#' @param keep_reps Logical; return the per-replicate output.
#' @return A list with `summary` (parameter-level calibration per cell),
#'   `resi_summary` (RESI calibration per cell and term), `truth_fit`,
#'   `truth_resi`, and optionally `reps` (per cell).
#' @export
hybridPlasmodeSim <- function(formula, data, school = "schoolid", child = "childid",
                              n_schools = c(15, 30, 60, 100),
                              n_per_site = c(25, 100, 200),
                              nsim = 1000L, mc.cores = 1L, seed = 2026,
                              control = lme4::lmerControl(optimizer = "bobyqa",
                                                          optCtrl = list(maxfun = 2e5)),
                              rep_fun = NULL, keep_reps = !is.null(rep_fun)) {
  data <- as.data.frame(data)
  data[[school]] <- factor(data[[school]])
  data[[child]] <- factor(data[[child]])
  school_grp <- school
  child_grp <- paste0(child, ":", school)
  response <- all.vars(formula[[2L]])
  fixed_rhs <- stats::delete.response(stats::terms(reformulas::nobars(formula)))

  truth_fit <- lme4::lmer(formula, data = data, REML = FALSE, control = control)
  true_params <- .mixedSimParams(truth_fit)
  truth_resi <- lyuLimitResi(truth_fit, school_grp, child_grp)
  beta <- lme4::fixef(truth_fit)
  vc <- .mixedVarComps(truth_fit, school_grp, child_grp)
  resid_sd <- stats::sigma(truth_fit)

  # BLUP pools centred to mean zero and rescaled so their (1/G) covariance
  # equals the fitted component, undoing BLUP shrinkage; draws are with
  # replacement (whole rows for vector-valued child effects).
  msqrt <- function(S, inverse = FALSE) {
    e <- eigen(S, symmetric = TRUE)
    d <- pmax(e$values, 0)
    d <- if (inverse) ifelse(d > 1e-12 * max(d), 1 / sqrt(d), 0) else sqrt(d)
    e$vectors %*% (d * t(e$vectors))
  }
  rescale_blup <- function(blup, target_cov) {
    centered <- scale(as.matrix(blup), scale = FALSE)
    centered %*% msqrt(crossprod(centered) / nrow(centered), inverse = TRUE) %*%
      msqrt(as.matrix(target_cov))
  }
  re <- lme4::ranef(truth_fit)
  school_pool <- drop(rescale_blup(re[[school_grp]][[1L]], vc$tau2_school))
  child_pool <- rescale_blup(re[[child_grp]], vc$G_child)
  child_bar <- Filter(function(b) identical(deparse(b[[3L]]), child_grp),
                      reformulas::findbars(formula))
  if (length(child_bar) != 1L) stop("Expected one random-effects term for ", child_grp, ".")
  child_z_formula <- stats::as.formula(call("~", child_bar[[1L]][[2L]]))

  school_profiles <- levels(data[[school]])
  rows_by_school <- split(seq_len(nrow(data)), data[[school]])
  children_by_school <- lapply(rows_by_school, function(idx) unique(as.character(data[[child]][idx])))
  rows_by_child <- split(seq_len(nrow(data)), as.character(data[[child]]))
  mean_site_size <- mean(lengths(children_by_school))

  sample_design <- function(g, m) {
    source_schools <- sample(school_profiles, g, replace = TRUE)
    pieces <- lapply(seq_len(g), function(s) {
      child_ids <- children_by_school[[source_schools[s]]]
      k <- max(1L, as.integer(round(m * length(child_ids) / mean_site_size)))
      selected <- sample(child_ids, k, replace = TRUE)
      idx <- rows_by_child[selected]
      list(rows = unlist(idx, use.names = FALSE),
           child = rep(paste0(s, "_", seq_len(k)), lengths(idx)),
           school = rep(as.character(s), sum(lengths(idx))))
    })
    design <- data[unlist(lapply(pieces, `[[`, "rows")), , drop = FALSE]
    # Levels in order of appearance: random effects are indexed by level code.
    school_ids <- unlist(lapply(pieces, `[[`, "school"))
    child_ids <- unlist(lapply(pieces, `[[`, "child"))
    design[[school]] <- factor(school_ids, levels = unique(school_ids))
    design[[child]] <- factor(child_ids, levels = unique(child_ids))
    rownames(design) <- NULL
    design
  }

  generate_outcome <- function(design) {
    school_effect <- sample(school_pool, nlevels(design[[school]]), replace = TRUE)
    child_effect <- child_pool[sample.int(nrow(child_pool), nlevels(design[[child]]),
                                          replace = TRUE), , drop = FALSE]
    child_part <- rowSums(stats::model.matrix(child_z_formula, data = design) *
                            child_effect[design[[child]], , drop = FALSE])
    fixed_part <- stats::model.matrix(fixed_rhs, data = design) %*% beta
    drop(fixed_part) + school_effect[design[[school]]] + child_part +
      stats::rnorm(nrow(design), sd = resid_sd)
  }

  fe_names <- paste0("beta_", names(beta))
  between <- truth_resi$between
  # Lyu-Welsh per-unit SE reference sqrt(n_k V_kk).
  truth_se_unit <- setNames(sqrt(diag(truth_resi$V)), fe_names)
  param_names <- names(true_params)
  is_fe <- param_names %in% fe_names

  grid <- expand.grid(n_schools = n_schools, n_per_site = n_per_site)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  rows <- resi_rows <- all_reps <- vector("list", nrow(grid))
  for (cell in seq_len(nrow(grid))) {
    g <- grid$n_schools[cell]
    m <- grid$n_per_site[cell]
    reps <- parallel::mclapply(seq_len(nsim), function(i) {
      design <- sample_design(g, m)
      design[[response]] <- generate_outcome(design)
      fit <- tryCatch(suppressWarnings(suppressMessages(
        lme4::lmer(formula, data = design, REML = TRUE, control = control))),
        error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      # A resampled design can miss a rare factor level; such replicates are
      # counted as failed fits rather than estimating a different model.
      if (length(lme4::fixef(fit)) != length(beta)) return(NULL)
      point <- resi_pe(fit, anova = TRUE, vcovfunc = stats::vcov, level = 2)
      list(params = .mixedSimParams(fit),
           fe_se = setNames(sqrt(diag(as.matrix(stats::vcov(fit)))), fe_names),
           resi_coef = point$coefficients, resi_anova = point$anova,
           extra = if (!is.null(rep_fun)) rep_fun(fit, design))
    }, mc.cores = mc.cores, mc.set.seed = TRUE)
    reps <- Filter(Negate(is.null), reps)

    param_mat <- do.call(rbind, lapply(reps, function(r) r$params[param_names]))
    mc_mean <- colMeans(param_mat, na.rm = TRUE)
    mc_sd <- apply(param_mat, 2L, stats::sd, na.rm = TRUE)
    mean_model_se <- colMeans(do.call(rbind, lapply(reps, function(r) r$fe_se[fe_names])), na.rm = TRUE)
    n_scaling <- setNames(ifelse(between, g, g * m), fe_names)
    mc_sd_unit <- mc_sd[fe_names] * sqrt(n_scaling)

    rows[[cell]] <- data.frame(
      n_schools = g, n_per_site = m, n_success = length(reps),
      param = param_names,
      truth = true_params[param_names],
      mc_mean = mc_mean[param_names],
      bias = true_params[param_names] - mc_mean[param_names],
      mc_sd = mc_sd[param_names],
      mean_model_se = ifelse(is_fe, mean_model_se[param_names], NA_real_),
      se_ratio = ifelse(is_fe, mean_model_se[param_names] / mc_sd[param_names], NA_real_),
      se_n_level = ifelse(is_fe, ifelse(between[match(param_names, fe_names)], "schools", "students"), NA),
      n_scaling = unname(n_scaling[param_names]),
      truth_se_scaled = unname(truth_se_unit[param_names]),
      mc_error_sd_scaled = unname(mc_sd_unit[param_names]),
      scaled_mc_error_sd_ratio = unname(mc_sd_unit[param_names] / truth_se_unit[param_names]),
      row.names = NULL, stringsAsFactors = FALSE)
    resi_rows[[cell]] <- .hybridResiSummary(reps, truth_resi,
                                            list(n_schools = g, n_per_site = m))
    if (keep_reps) all_reps[[cell]] <- reps
    message("Finished cell ", cell, "/", nrow(grid), ": n_schools = ", g,
            ", n_per_site = ", m, " (", length(reps), " successful fits)")
  }
  out <- list(summary = do.call(rbind, rows), resi_summary = do.call(rbind, resi_rows),
              truth_fit = truth_fit, truth_resi = truth_resi, grid = grid)
  if (keep_reps) out$reps <- all_reps
  out
}
