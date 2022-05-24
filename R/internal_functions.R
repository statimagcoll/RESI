# Internal functions

#' Transferation from Wald test statistics to squared RESI
#' @export
#' @param chisq The chi-square statistic for the parameter(s) of interest.
#' @param df Number of degrees of freedom of the chi-square statistic.
#' @param rdf Model residual degrees of freedom or number of independent samples.
#' @return Returns a scalar or vector argument of the squared robust effect size index estimate.
chisq2Ssq = function(chisq, df, rdf){
  S = (chisq - df)/rdf
}


#' Non-parametric bootstrap sampling
#' @export
#' @param data data frame; The data frame that need bootstrapping
#' @param id.var character; for clustered/longitudinal data, the name of id variable used as sampling unit.
#' @return Returns a data frame containing bootstrapped data.
boot.samp <- function(data, id.var = NULL) {
  params = as.list(match.call()[-1])
  if (is.matrix(data)) data = as.data.frame(data)
  if (is.null(id.var)) {
    boot.ind = sample(1:nrow(data), replace = TRUE)
    boot.data = data[boot.ind, ]
  } else {
    boot.ind = sample(unique(data[, id.var]), replace = TRUE)
    boot.data = data[unlist(lapply(boot.ind, function(x) which(x == data[, id.var]))), ]
  }
  return(boot.data)
}

#' Bayesian bootstrap sampling (Rubin, 1981)
#' @export
#' @param data data frame; The data frame that need bootstrapping
#' @return Returns a data frame containing weights generated via Bayesian bootstraps.
bayes.samp <- function(data) {
  if (is.matrix(data)) data = as.data.frame(data)
  n = nrow(data)
  repeat{
    # Generate the random numbers from unif(0, 1)
    u = runif(n-1)
    u.sort = sort(u)
    g = c(u.sort, 1) - c(0, u.sort)
    if (sum(g == 0) == 0) break
  } # end `repeat`
  boot.data = cbind(data, g)
  return(boot.data)
}



print.resi <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (x$robust.var) cat("\n Analysis of Effect sizes (ANOES) based on RESI using robust sandwich covariance estimator: ")
  if (! x$robust.var) cat("\n Analysis of Effect sizes (ANOES) based on RESI using naive covariance estimator: ")
  cat("\n significance level = ", x$alpha)
  cat("\n Call:  ",
      paste(deparse(x$input.object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if(length(coef(x$input.object))) {
    print.default(format(x$resi, digits = digits),
                  print.gap = 2, quote = FALSE)
  } else cat("No coefficients\n\n")

  if (x$boot.method == "nonparam") cat("\n Confidence intervals (CIs) constructed using", x$nboot,"non-parametric bootstraps")
  if (x$boot.method == "bayes") cat("\n Credible intervals constructed using", x$nboot,"Bayesian bootstraps")
  if(nzchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep = "")

  invisible(x)
}

