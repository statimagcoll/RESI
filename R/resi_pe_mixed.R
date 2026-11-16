# [AI-generated: GitHub Copilot, 2026-09-29]
# [Prompt summary: mixed-model RESI with an explicit, term-specific normalizing n
#  independent of the cluster used by the robust variance estimator.]

.resiMixedModelDesign <- function(model.full, engine) {
  if (engine == "lmerMod") {
    X <- lme4::getME(model.full, "X")
    groups <- as.list(lme4::getME(model.full, "flist"))
    data <- tryCatch({
      data_arg <- model.full@call$data
      if (is.null(data_arg)) stop("no stored data call")
      d <- as.data.frame(eval(data_arg, envir = environment(formula(model.full))))
      mf <- model.frame(model.full)
      used_rows <- rownames(mf)
      if (!is.null(used_rows) && all(used_rows %in% rownames(d))) {
        d[used_rows, , drop = FALSE]
      } else {
        as.data.frame(mf)
      }
    }, error = function(e) as.data.frame(model.frame(model.full)))
  } else {
    data <- as.data.frame(nlme::getData(model.full))
    fixed_formula <- formula(model.full)
    X <- stats::model.matrix(stats::terms(fixed_formula), data = data)
    dims <- summary(model.full)$dims
    inner_to_outer <- rev(seq_len(dims$Q))
    groups <- lapply(inner_to_outer, function(k) nlme::getGroups(model.full, level = k))
    names(groups) <- names(dims$ngrps)[seq_len(dims$Q)]
  }
  if (nrow(X) != length(groups[[1]])) {
    stop("The fixed-effect design and random-effect grouping rows do not align.")
  }
  list(X = X, assign = attr(X, "assign"), terms = attr(stats::terms(
    if (engine == "lmerMod") reformulas::nobars(formula(model.full)) else formula(model.full)
  ), "term.labels"), groups = groups, data = data)
}

.resiMixedHasVariation <- function(x) {
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  if (nrow(x) < 2L) return(FALSE)
  any(vapply(seq_len(ncol(x)), function(j) {
    vals <- x[, j]
    diff(range(vals, na.rm = TRUE)) > 1e-10 * max(1, max(abs(vals), na.rm = TRUE))
  }, logical(1L)))
}

.resiMixedTermLevel <- function(term_matrix, groups, term_name) {
  n_levels <- length(groups)
  if (n_levels == 0L) stop("Mixed-model RESI requires at least one grouping factor.")

  for (k in seq_len(n_levels)) {
    group <- droplevels(as.factor(groups[[k]]))
    row_groups <- split(seq_along(group), group, drop = TRUE)
    if (any(vapply(row_groups, function(idx) .resiMixedHasVariation(term_matrix[idx, , drop = FALSE]), logical(1L)))) {
      return(list(variation_level = paste0("observation within ", names(groups)[[k]]),
                  n_level = names(groups)[[k]]))
    }

    group_means <- do.call(rbind, lapply(row_groups, function(idx) {
      colMeans(term_matrix[idx, , drop = FALSE])
    }))
    if (k == n_levels) {
      if (.resiMixedHasVariation(group_means)) {
        return(list(variation_level = names(groups)[[k]], n_level = names(groups)[[k]]))
      }
    } else {
      parent <- droplevels(as.factor(groups[[k + 1L]]))
      parent_by_group <- vapply(row_groups, function(idx) {
        values <- unique(as.character(parent[idx]))
        if (length(values) != 1L) NA_character_ else values[[1L]]
      }, character(1L))
      if (anyNA(parent_by_group)) {
        stop("Cannot infer a term's variation level for crossed/non-nested grouping factors; supply `n` explicitly.")
      }
      sibling_groups <- split(seq_len(nrow(group_means)), parent_by_group, drop = TRUE)
      if (any(vapply(sibling_groups, function(idx) {
        .resiMixedHasVariation(group_means[idx, , drop = FALSE])
      }, logical(1L)))) {
        return(list(variation_level = names(groups)[[k]], n_level = names(groups)[[k]]))
      }
    }
  }
  stop("Cannot infer a variation level for term '", term_name,
       "' because it is constant across the observed design; supply `n` explicitly.")
}

# Descriptive decomposition of a term's design-matrix variance across the
# nested grouping levels: grand mean -> outermost group means -> ... ->
# innermost group means -> observations. Sums of squares are taken over
# observations, so the components are orthogonal and add to the total; the
# observation-level remainder (1 - sum of returned proportions) is omitted.
# Multi-column terms pool sums of squares across columns. Returns NA when the
# term is constant or the grouping factors are not nested.
.resiMixedVarianceProportions <- function(term_matrix, groups) {
  term_matrix <- as.matrix(term_matrix)
  out <- setNames(rep(NA_real_, length(groups)), names(groups))
  centered <- sweep(term_matrix, 2L, colMeans(term_matrix))
  total <- sum(centered^2)
  if (!is.finite(total) || total <= 1e-12 * max(1, nrow(term_matrix))) return(rev(out))
  group_means <- lapply(groups, function(g) {
    g <- droplevels(as.factor(g))
    apply(term_matrix, 2L, function(x) ave(x, g))
  })
  for (k in seq_along(groups)) {
    if (k < length(groups)) {
      inner <- droplevels(as.factor(groups[[k]]))
      outer <- droplevels(as.factor(groups[[k + 1L]]))
      if (any(tapply(as.character(outer), inner, function(x) length(unique(x))) > 1L)) {
        return(rev(out))
      }
      parent <- group_means[[k + 1L]]
    } else {
      parent <- matrix(colMeans(term_matrix), nrow(term_matrix), ncol(term_matrix), byrow = TRUE)
    }
    out[[k]] <- sum((group_means[[k]] - parent)^2) / total
  }
  rev(out)
}

.resiMixedScaleMap <- function(design, n = NULL, level = NULL) {
  X <- as.matrix(design$X)
  assign <- design$assign
  term_labels <- design$terms
  groups <- design$groups
  group_counts <- vapply(groups, function(g) nlevels(droplevels(as.factor(g))), integer(1L))
  hierarchy_names <- rev(names(groups))
  hierarchy_counts <- rev(group_counts)
  term_levels <- setNames(vector("list", length(term_labels)), term_labels)

  for (term_id in seq_along(term_labels)) {
    cols <- which(assign == term_id)
    term_levels[[term_labels[[term_id]]]] <- .resiMixedTermLevel(
      X[, cols, drop = FALSE], groups, term_labels[[term_id]])
  }
  term_levels <- c(list("(Intercept)" = list(
    variation_level = "overall/intercept",
    n_level = names(groups)[[1L]]
  )), term_levels)
  variation_levels <- vapply(term_levels, `[[`, character(1L), "variation_level")
  n_levels <- vapply(term_levels, `[[`, character(1L), "n_level")

  if (is.null(level)) {
    level <- 1L
    message("No `level` supplied; using level 1 (", hierarchy_names[[1L]],
            ") to normalize within-", hierarchy_names[[1L]], " terms.")
  }
  if (is.character(level) && length(level) == 1L) {
    if (!level %in% hierarchy_names) {
      stop("`level` must name a random-effect grouping factor or be a numeric hierarchy level.")
    }
    level_index <- match(level, hierarchy_names)
  } else if (is.numeric(level) && length(level) == 1L &&
             is.finite(level) && level == as.integer(level) &&
             level >= 1L && level <= length(hierarchy_names)) {
    level_index <- as.integer(level)
  } else {
    stop("`level` must be one numeric hierarchy level or a random-effect grouping name.")
  }
  selected_level <- hierarchy_names[[level_index]]
  selected_n <- hierarchy_counts[[level_index]]

  if (!is.null(n)) {
    if (length(n) != 1L || !is.numeric(n) || is.na(n) ||
        !is.finite(n) || n <= 0) {
      stop("`n` must be one positive finite numeric normalizing count.")
    }
    selected_n <- n[[1L]]
    selected_level <- "user-specified n"
    # When n is user-supplied, use it for all terms (backward compatibility)
    term_ns <- setNames(rep(selected_n, length(term_levels)), names(term_levels))
    n_levels[] <- selected_level
  } else {
    # Lyu & Welsh (2022): terms constant within the outermost cluster
    # (intercept, cluster-level covariates) are normalized by the number of
    # clusters g; all other terms by the number of units at `level`.
    cluster_count <- hierarchy_counts[[1L]]
    within_cluster_count <- hierarchy_counts[[level_index]]
    top_cluster_name <- hierarchy_names[[1L]]
    is_between_cluster <- variation_levels %in% c("overall/intercept", top_cluster_name)
    
    term_ns <- setNames(
      ifelse(
        is_between_cluster,
        cluster_count,
        within_cluster_count
      ),
      names(term_levels)
    )
    n_levels <- ifelse(is_between_cluster, top_cluster_name, selected_level)
    names(n_levels) <- names(term_levels)
  }

  coefficient_terms <- character(ncol(X))
  for (j in seq_len(ncol(X))) {
    term_id <- assign[[j]]
    coefficient_terms[[j]] <- if (term_id == 0L) "(Intercept)" else term_labels[[term_id]]
  }
  names(coefficient_terms) <- colnames(X)

  prop_names <- paste("Var prop", rev(names(groups)))
  coefficient_var_prop <- t(vapply(seq_len(ncol(X)), function(j) {
    .resiMixedVarianceProportions(X[, j, drop = FALSE], groups)
  }, numeric(length(groups))))
  if (length(groups) == 1L) coefficient_var_prop <- t(coefficient_var_prop)
  dimnames(coefficient_var_prop) <- list(colnames(X), prop_names)
  term_var_prop <- t(vapply(names(term_levels), function(term) {
    cols <- if (term == "(Intercept)") which(assign == 0L) else which(assign == match(term, term_labels))
    .resiMixedVarianceProportions(X[, cols, drop = FALSE], groups)
  }, numeric(length(groups))))
  if (length(groups) == 1L) term_var_prop <- t(term_var_prop)
  dimnames(term_var_prop) <- list(names(term_levels), prop_names)

    list(variation_level = variation_levels, n_level = n_levels,
      n = term_ns, coefficient_term = coefficient_terms,
      coefficient_var_prop = coefficient_var_prop, term_var_prop = term_var_prop,
      group_counts = group_counts, level = level_index,
      level_name = selected_level)
}

.resiMixedCluster <- function(cluster, design, model.full, engine) {
  if (is.null(cluster)) return(NULL)
  if (is.character(cluster) && length(cluster) == 1L) {
    if (cluster %in% names(design$groups)) return(design$groups[[cluster]])
    if (cluster %in% names(design$data)) return(design$data[[cluster]])
    stop("`cluster` must name a model grouping factor or a column in the fitted data.")
  }
  if (length(cluster) != nrow(design$X)) {
    stop("A supplied `cluster` vector must have one value per model observation.")
  }
  cluster
}

.resiMixedVcov <- function(model.full, vcovfunc, vcov.args, cluster) {
  if (!is.null(cluster) && "cluster" %in% names(vcov.args)) {
    stop("Specify the robust-variance cluster either with `cluster` or `vcov.args$cluster`, not both.")
  }
  if (identical(vcovfunc, clubSandwich::vcovCR) && !("type" %in% names(vcov.args))) {
    vcov.args$type <- "CR3"
  }
  if (!is.null(cluster)) vcov.args$cluster <- cluster
  fast_cr2 <- identical(vcovfunc, clubSandwich::vcovCR) &&
    identical(vcov.args$type, "CR2") && all(names(vcov.args) %in% c("type", "cluster")) &&
    !is.null(vcov.args$cluster) &&
    ((inherits(model.full, "lmerMod") && all(stats::weights(model.full) == 1)) ||
       (inherits(model.full, "lme") && is.null(.mixedCiLmeUnsupported(model.full))))
  function(object) {
    if (fast_cr2) {
      out <- tryCatch(.mixedCiCR2(object, vcov.args$cluster), error = function(e) NULL)
      if (!is.null(out)) return(out)
    }
    do.call(vcovfunc, c(list(object), vcov.args))
  }
}

.resi_pe_mixed <- function(model.full, engine, anova, vcovfunc, Anova.args,
                           vcov.args, unbiased, n, level, cluster, include_level) {
  design <- .resiMixedModelDesign(model.full, engine)
  scale_map <- .resiMixedScaleMap(design, n, level)
  cluster_vec <- .resiMixedCluster(cluster, design, model.full, engine)
  vcovfun <- .resiMixedVcov(model.full, vcovfunc, vcov.args, cluster_vec)
  coefficient_summary <- if (engine == "lmerMod") {
    as.matrix(summary(model.full)$coefficients)
  } else {
    as.matrix(summary(model.full)$tTable)
  }
  estimate_name <- if (engine == "lmerMod") "Estimate" else "Value"
  estimate <- coefficient_summary[, estimate_name]
  coefficient_terms <- scale_map$coefficient_term[rownames(coefficient_summary)]
  coefficient_n <- unname(scale_map$n[coefficient_terms])

  naive <- identical(vcovfunc, stats::vcov)
  if (naive) {
    covariance <- stats::vcov(model.full)
    se <- sqrt(diag(as.matrix(covariance)))[rownames(coefficient_summary)]
  } else {
    covariance <- vcovfun(model.full)
    se <- sqrt(diag(covariance))[rownames(coefficient_summary)]
  }
  statistic <- estimate / se
  coefficient_table <- as.data.frame(coefficient_summary)
  if (!naive) coefficient_table$Robust.SE <- unname(se)
  coefficient_table$Wald <- statistic^2
  coefficient_table$RESI <- z2S(statistic, coefficient_n, unbiased)
  coefficient_table$Term <- unname(coefficient_terms)
  coefficient_table$n <- coefficient_n
  if (include_level) {
     coefficient_table$`Variation level` <- unname(scale_map$variation_level[coefficient_terms])
     coefficient_table$`n level` <- unname(scale_map$n_level[coefficient_terms])
     coefficient_table <- cbind(coefficient_table, as.data.frame(
       scale_map$coefficient_var_prop[rownames(coefficient_summary), , drop = FALSE],
       check.names = FALSE))
  }

  output <- list(estimates = setNames(coefficient_table$RESI, rownames(coefficient_table)),
                 coefficients = coefficient_table, naive.var = naive)

  if (anova) {
    if (naive) {
      anova_table <- suppressMessages(do.call(car::Anova,
        c(list(mod = model.full), Anova.args)))
    } else {
      anova_table <- suppressMessages(do.call(car::Anova,
        c(list(mod = model.full, vcov. = vcovfun), Anova.args)))
    }
    anova_terms <- rownames(anova_table)
    missing_terms <- setdiff(anova_terms, names(scale_map$n))
    if (length(missing_terms)) {
      stop("Could not map ANOVA terms to fixed-effect terms: ",
           paste(missing_terms, collapse = ", "))
    }
    anova_n <- unname(scale_map$n[anova_terms])
    anova_table$RESI <- chisq2S(anova_table$Chisq, anova_table$Df, anova_n)
    anova_table$n <- anova_n
    if (include_level) {
      anova_table$`Variation level` <- unname(scale_map$variation_level[anova_terms])
      anova_table$`n level` <- unname(scale_map$n_level[anova_terms])
      for (prop_name in colnames(scale_map$term_var_prop)) {
        anova_table[[prop_name]] <- unname(scale_map$term_var_prop[anova_terms, prop_name])
      }
    }
    output$estimates <- c(output$estimates, setNames(anova_table$RESI, anova_terms))
    output$anova <- anova_table
    class(output$anova) <- c("anova_resi", class(output$anova))
  }
  output$model.full <- list(call = if (engine == "lmerMod") model.full@call else model.full$call,
                            formula = formula(model.full))
  class(output) <- c("resi", "list")
  output
}

#' @describeIn resi_pe Mixed-model RESI point estimation
#' @param n Optional positive scalar normalizing count applied to every
#'   coefficient and ANOVA term. If `NULL` (default), the count is
#'   term-specific (Lyu & Welsh, 2022): terms constant within the outermost
#'   grouping factor (including the intercept) use the number of outermost
#'   groups, and all other terms use the number of groups at `level`.
#' @param level Random-effect hierarchy level or grouping-factor name whose
#'   group count normalizes terms that vary within the outermost grouping
#'   factor. Numeric levels are ordered from outermost to innermost grouping
#'   (for `egsingle`, school is level 1 and student is level 2). If omitted,
#'   level 1 is used and a message names the grouping factor.
#' @param cluster Optional grouping-factor name or observation-length vector
#'   specifying the cluster used only by `vcovfunc`. `NULL` leaves the
#'   variance function's default unchanged. For `clubSandwich::vcovCR`, this
#'   is commonly the outermost grouping factor; it need not equal the level
#'   used to infer `n`.
#' @param include_level Logical; add `Variation level` and `n level` columns
#'   to coefficient and ANOVA tables (default `TRUE`). `Variation level`
#'   reports the lowest observed design level at which a term changes;
#'   `n level` reports the grouping factor whose count is the term's `n`. These can differ (e.g., `year` varies by occasion within a
#'   student, while its default `n` counts students). Explicit `n` overrides
#'   retain the design-derived variation label and mark `n level` as
#'   user-specified. Also adds one `Var prop <group>` column per grouping
#'   factor (outermost first): the share of the fixed-effect design column's
#'   (coefficients) or term's (ANOVA; pooled across its columns) observation-
#'   level sum of squares lying between units of that grouping factor within
#'   the next-outer one. The remainder, `1 - sum`, is the share varying within
#'   the innermost group. Values are `NA` for the intercept and for
#'   non-nested grouping factors.
#' @importFrom clubSandwich vcovCR
#' @export
resi_pe.lmerMod <- function(model.full, anova = TRUE,
                            vcovfunc = clubSandwich::vcovCR,
                            Anova.args = list(), vcov.args = list(),
                            unbiased = TRUE, n = NULL, level = NULL, cluster = NULL,
                            include_level = TRUE, ...) {
  .resi_pe_mixed(model.full, "lmerMod", anova, vcovfunc, Anova.args,
                 vcov.args, unbiased, n, level, cluster, include_level)
}

#' @describeIn resi_pe Mixed-model RESI point estimation
#' @export
resi_pe.lme <- function(model.full, anova = TRUE,
                        vcovfunc = clubSandwich::vcovCR,
                        Anova.args = list(), vcov.args = list(),
         unbiased = TRUE, n = NULL, level = NULL, cluster = NULL,
                        include_level = TRUE, ...) {
  .resi_pe_mixed(model.full, "lme", anova, vcovfunc, Anova.args,
       vcov.args, unbiased, n, level, cluster, include_level)
}
