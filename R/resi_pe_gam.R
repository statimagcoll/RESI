#' @describeIn resi_pe RESI point estimation for generalized additive models
#'   (\code{mgcv::gam} and \code{mgcv::bam} fits)
#' @details For \code{gam} and \code{bam} models, the RESI is computed from
#' the Wald statistics reported by \code{mgcv::summary.gam}. These use the
#' Bayesian posterior covariance \code{Vp}, so \code{vcovfunc} is not used.
#' Parametric coefficients use \code{\link{t2S}} (or \code{\link{z2S}} when
#' the scale parameter is known). Smooth terms use \code{\link{f2S}} with the
#' reference degrees of freedom \code{Ref.df} and the residual degrees of
#' freedom \code{n - sum(edf)} (or \code{\link{chisq2S}} when the scale is
#' known). The smooth-term RESI estimates the RESI of the population
#' smooth function (net of the other terms); when that function is well
#' represented by the basis, it does not depend on the basis size or the
#' penalty. Penalized terms without an unpenalized null space
#' (\code{bs = "re"}, shrinkage bases, or \code{select = TRUE}) are tested by
#' mgcv with a different statistic, so their RESI is \code{NA}. The overall
#' test is not computed for these models.
#' @export
resi_pe.gam <- function(model.full, data, anova = TRUE, coefficients = TRUE,
                        unbiased = TRUE, ...) {
  if (!requireNamespace("mgcv", quietly = TRUE))
    stop("\nPackage 'mgcv' is required for gam and bam models")
  s <- mgcv::summary.gam(model.full)
  n <- length(model.full$y)
  rdf <- model.full$df.residual
  output <- list(model.full = list(call = model.full$call,
                                   formula = formula(model.full)),
                 model.reduced = NULL)
  est <- c()

  if (coefficients) {
    pt <- as.data.frame(s$p.table)
    pt$RESI <- if ("t value" %in% colnames(pt)) {
      suppressWarnings(t2S(pt[, "t value"], rdf, n, unbiased))
    } else {
      suppressWarnings(z2S(pt[, "z value"], n, unbiased))
    }
    output$coefficients <- pt
    est <- c(est, stats::setNames(pt$RESI, rownames(pt)))
  }

  if (anova) {
    stat_name <- if (isTRUE(model.full$scale.estimated)) "F" else "Chisq"
    pv_name <- if (stat_name == "F") "Pr(>F)" else "Pr(>Chisq)"
    stat_col <- if (stat_name == "F") "F" else "Chi.sq"
    rows <- list()
    if (!is.null(s$pTerms.table) && nrow(s$pTerms.table) > 0) {
      p <- s$pTerms.table
      rows$p <- data.frame(edf = NA_real_, Df = p[, "df"], stat = p[, stat_col],
                           pval = p[, "p-value"], row.names = rownames(p))
    }
    if (!is.null(s$s.table) && nrow(s$s.table) > 0) {
      sm <- s$s.table
      rows$s <- data.frame(edf = sm[, "edf"], Df = sm[, "Ref.df"], stat = sm[, stat_col],
                           pval = sm[, "p-value"], row.names = rownames(sm))
    }
    tab <- do.call(rbind, unname(rows))
    if (!is.null(tab)) {
      tab$RESI <- if (stat_name == "F") {
        f2S(tab$stat, tab$Df, rdf, n)
      } else {
        chisq2S(tab$stat, tab$Df, n)
      }
      # mgcv tests smooths without an unpenalized null space with a
      # different statistic (Wood, 2013b), which has no RESI conversion
      # (bs = "re", shrinkage bases, select = TRUE), so these get NA.
      re_test <- vapply(model.full$smooth, function(sm) {
        !isTRUE(sm$fixed) && identical(as.integer(sm$null.space.dim), 0L)
      }, logical(1L))
      no_null <- vapply(model.full$smooth, `[[`, character(1L), "label")[re_test]
      if (length(no_null)) {
        tab[rownames(tab) %in% no_null, "RESI"] <- NA_real_
        warning("\nRESI is not computed for smooth terms without an unpenalized ",
                "null space: ", paste(no_null, collapse = ", "))
      }
      colnames(tab)[3:4] <- c(stat_name, pv_name)
      attr(tab, "heading") <- "Analysis of Deviance Table (Wald tests from mgcv::summary.gam)\n"
      class(tab) <- c("anova_resi", "anova", "data.frame")
      output$anova <- tab
      est <- c(est, stats::setNames(tab$RESI, rownames(tab)))
    }
  }

  output$estimates <- est
  output$naive.var <- TRUE
  class(output) <- c("resi", "list")
  output
}

#' @describeIn resi RESI point estimates for generalized additive models.
#'   Confidence intervals are not yet available for \code{gam} and \code{bam}
#'   fits; use \code{\link{resi_pe}}.
#' @export
resi.gam <- function(model.full, ...) {
  stop("\nConfidence intervals are not yet implemented for gam and bam models. ",
       "Use resi_pe() for point estimates.")
}
