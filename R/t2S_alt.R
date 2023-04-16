#' Compute the robust effect size index estimate from t statistic (alternative)
#'
#' This function computes the robust effect size index from Vandekar, Tao, & Blume (2020).
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' @param t The t statistic for the parameter of interest.
#' @param rdf Model residual degrees of freedom or number of independent samples.
#' @param n Number of independent samples.
#' @details This function computes S, the RESI, from a T statistic using the alternative
#' formula. There is another function, \code{\link{t2S}}, that is the default for
#' \code{\link{resi}}. This function's formula is derived by squaring the T statistic
#' and using the \code{\link{f2S}} formula. This function may be appealing for its
#' intuitive relationship to the F statistic; the absolute value of RESI estimates
#' using this formula will be equal to a RESI estimate using an F statistic for
#' the same model. However, this estimator does have finite sample bias, which is an
#' important consideration for the coverage of the bootstrapping that \code{resi} uses.
#'
#' The formula for this conversion is:
#'
#' \eqn{ \sqrt(max(0, (t^2 * (rdf - 2)/rdf - 1)/rdf))}
#' @examples
#' # to obtain t values, first fit a lm
#' mod = lm(charges ~ region * age + bmi + sex, data = RESI::insurance)
#' # run lmtest::coeftest to get t values, using a robust variance-covariance formula
#' ts = lmtest::coeftest(mod, vcov. = sandwich::vcovHC)[,'t value']
#'
#' # get RESI estimates
#' t2S_alt(ts, rdf = mod$df.residual, n = nrow(RESI::insurance))
#' @return Returns a scalar or vector argument of the the robust effect size index estimate.
#' @export
t2S_alt <- function(t, rdf, n){
  Ssq = (t^2*(rdf-2)/rdf - 1)/n
  sqrt(ifelse(Ssq<0, 0, Ssq))*(sqrt(t^2)/t)
}
