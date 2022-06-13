# Internal functions

#' Transformation from Wald test statistics to squared RESI
#' @param chisq The chi-square statistic for the parameter(s) of interest.
#' @param df Number of degrees of freedom of the chi-square statistic.
#' @param rdf Model residual degrees of freedom or number of independent samples.
#' @return Returns a scalar or vector argument of the squared robust effect size index estimate.
#' @export
chisq2Ssq = function(chisq, df, rdf){
  S = (chisq - df)/rdf
}

#' Compute the robust effect size index estimator from chi-squared statistic.
#'
#' This function computes the robust effect size index from Vandekar, Rao, & Blume (2020).
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' For mixed effects models, RESI is conditional on the average correlation
#' structure within subjects.
#' @param chisq The chi-square statistic for the parameter of interest.
#' @param df Number of degrees of freedom of the chi-square statistic.
#' @param n Number of independent samples.
#' @keywords power
#' @return Returns a scalar or vector argument of the the robust effect size index estimate.
#' @export
chisq2S <- function(chisq, df, n){
  S = (chisq - df)/n
  sqrt(ifelse(S<0, 0, S))
}

#' Compute the robust effect size index estimator from t statistic
#'
#' This function computes the robust effect size index from Vandekar, Rao, & Blume (2020).
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' @param t The t statistic for the parameter of interest.
#' @param rdf Model residual degrees of freedom or number of independent samples.
#' @keywords power
#' @return Returns a scalar or vector argument of the the robust effect size index estimate.
#' @export
t2S <- function(t, rdf){
  2*t/rdf*exp(lgamma(rdf/2)- lgamma((rdf-1)/2))
}

#' Compute the robust effect size index estimator from Z statistic
#'
#' This function computes the robust effect size index from Vandekar, Rao, & Blume (2020).
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' @param z The Z statistic for the parameter of interest.
#' @param n Number of independent samples.
#' @keywords power
#' @return Returns a scalar or vector argument of the the robust effect size index estimate.
#' @export
z2S <- function(z, n){
  z/sqrt(n)
}


#' Compute the robust effect size index estimator from F-statistic
#'
#' This function computes the robust effect size index from Vandekar, Rao, & Blume (2020).
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' @param f The chi-square statistic for the parameter of interest.
#' @param df Number of degrees of freedom of the chi-square statistic.
#' @param rdf Model residual degrees of freedom or number of independent samples.
#' @keywords power
#' @return Returns a scalar or vector argument of the the robust effect size index estimate.
#' @export
f2S <- function(f, df, rdf){
  S = (f*df*(rdf-2)/rdf - df)/rdf
  sqrt(ifelse(S<0, 0, S))
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

#' @export
print.resi <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\n Analysis of Effect sizes (ANOES) based on RESI: ")
  cat("\n Significance level = ", x$alpha)
  cat( ifelse(is.null(x$input.object.reduced), "\n Call:  ", "\n Full Model:"),
      paste(deparse(x$input.object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (!is.null(x$input.object.reduced)) cat("\n Reduced Model:", paste(deparse(x$input.object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if(length(coef(x$input.object))) {
    print.default(format(round(x$resi, digits = digits)),
                  print.gap = 2, quote = FALSE)
  } else cat("No coefficients\n\n")

  cat("\n Notes:")
  if (x$robust.var) cat("\n 1. The RESI is calculated using the robust (sandwich) covariance estimator.")
  else cat("\n 1. The RESI is calculated using the naive covariance estimator.")
  if (x$boot.method == "nonparam") cat("\n 2. Confidence intervals (CIs) constructed using", x$nboot,"non-parametric bootstraps \n")
  if (x$boot.method == "bayes") cat("\n 2. Credible intervals constructed using", x$nboot,"Bayesian bootstraps \n")
  if(nzchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep = "")

  invisible(x)
}


#' Convert S to Cohen's d
#'
#' Converts the robust effect size index to Cohen's d using formula from Vandekar, Rao, & Blume (2020).
#' @param S The robust effect size index.
#' @param pi The sampling proportions.
#' @keywords power
#' @return Returns an estimate the robust effect size index
#' @details The pi parameter comes from the fact that Cohen's d doesn't account for unequal sample proportions in the population, but S does.
#' The default is set to a natural value 1/2, which corresponds to a case control design, for example, where sampling proportions always are controlled by the experimenter.
#' @export
S2d = function(S, pi=0.5){
  S * sqrt(1/pi + 1/(1-pi) )
}

