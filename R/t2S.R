#' Compute the robust effect size index estimate from t statistic (default)
#'
#' This function computes the robust effect size index from Vandekar, Rao, & Blume (2020).
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' @param t The t statistic for the parameter of interest.
#' @param n Number of independent samples.
#' @param rdf Model residual degrees of freedom/degrees of freedom of the T statistic.
#' @details This function computes S, the RESI, from a T statistic using the default
#' formula. There is another function, \code{\link{t2S_alt}}, that uses an alternative
#' formula. This function's formula is derived by solving the expected value of the
#' T statistic for S. It is unbiased and consistent.
#' @keywords power
#' @return Returns a scalar or vector argument of the the robust effect size index estimate.
#' @examples
#' # to obtain t values, first fit a lm
#' mod = lm(charges ~ region * age + bmi + sex, data = RESI::insurance)
#' # run lmtest::coeftest to get t values, using a robust variance-covariance formula
#' ts = lmtest::coeftest(mod, vcov. = sandwich::vcovHC)[,'t value']
#'
#' # get RESI estimates
#' t2S(ts, n = nrow(RESI::insurance), rdf = mod$df.residual)
#'
#' @export
t2S <- function(t, n, rdf){
  (t*sqrt(2))/(sqrt(n*rdf))*exp(lgamma(rdf/2) - lgamma((rdf-1)/2))
}
