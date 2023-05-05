#' @export
print.resi <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nAnalysis of effect sizes based on RESI:")
  cat("\nConfidence level = ", x$alpha)
  if (is.null(x$model.reduced$formula)){
    cat("\nCall:  ", paste(deparse(x$model.full$call), sep = "\n", collapse = "\n"),  "\n",sep = "")
  }
  else{
    cat(ifelse(x$model.reduced$formula == as.formula(paste(format(formula(x$model.full)[[2]]), "~ 1")), "\nCall:  ", "\nFull Model:"),
        paste(deparse(x$model.full$call), sep = "\n", collapse = "\n"), "\n",sep = "")
    if (!(x$model.reduced$formula == as.formula(paste(format(formula(x$model.full)[[2]]), "~ 1")))) cat("Reduced Model:", paste(deparse(x$model.reduced$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  }

  # coefficients table
  if (!is.null(x$coefficients)){
    cat("\nCoefficient Table \n")
    print(round(x$coefficients, digits = digits))
  }

  # anova table
  if (!is.null(x$anova)){
    cat("\n\n")
    print(round(x$anova, digits = digits))
  }

  # overall
  if (!(is.null(x$model.reduced))){
    if (is.null(x$model.reduced$formula)){
      cat("\nOverall RESI comparing model to intercept-only model:\n\n")
      print(round(x$overall, digits = digits))
    } else{
      if ((x$model.reduced$formula == as.formula(paste(format(formula(x$model.full)[[2]]), "~ 1")))){
        cat("\nOverall RESI comparing model to intercept-only model:\n\n")
      }
      else{
        cat("\nOverall RESI comparing full model to reduced model:\n\n")}
      overall = as.data.frame(x$overall)
      overall = overall[nrow(overall),]
      rownames(overall) = NULL
      print(round(overall, digits = digits))
      }
    }

  cat("\nNotes:")
  if (x$naive.var) cat("\n1. The RESI was calculated using the naive covariance estimator.")
  else  cat("\n1. The RESI was calculated using a robust covariance estimator.")
  if (x$boot.method == "nonparam") cat("\n2. Confidence intervals (CIs) constructed using", x$nboot,"non-parametric bootstraps. \n")
  if (x$boot.method == "bayes") cat("\n2. Credible intervals constructed using", x$nboot,"Bayesian bootstraps. \n")
  # if(nzchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep = "")
  # report number of failed bootstraps for nls model
  if (!(is.null(x$nfail))) cat("3. The bootstrap was successful in", x$nboot - x$nfail, "out of", x$nboot, "attempts. \n" )

  invisible(x)
}

#' @export
print.summary_resi <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nAnalysis of effect sizes based on RESI:")
  cat("\nConfidence level = ", x$alpha)
  cat("\nCall:  ", paste(deparse(x$model.full$call), sep = "\n", collapse = "\n"),  "\n",sep = "")
  cat("\nCoefficient Table \n")
  print(round(x$coefficients, digits = digits))
  invisible(x)
}
