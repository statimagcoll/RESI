#' Compute the robust effect size index estimator from chi-squared statistic.
#'
#' This function computes the robust effect size index from Vandekar, Rao, & Blume (2020).
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' For mixed effects models, RESI is conditional on the average correlation structure within subjects.
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

### A few conversion functions for other indices

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


