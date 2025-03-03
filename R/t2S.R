#' Compute the robust effect size index estimate from t statistic (default)
#'
#' This function computes the robust effect size index from Vandekar, Tao, & Blume (2020).
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' @param t The t statistic for the parameter of interest.
#' @param rdf Model residual degrees of freedom/degrees of freedom of the t statistic.
#' @param n Number of independent samples.
#' @param unbiased Logical, whether to use unbiased or alternative estimator. See details.
#' @details This function computes S, the RESI, from a t statistic. The formula for the
#' unbiased estimator (default) is derived by solving the expected value of the
#' t statistic for S. It is unbiased and consistent.
#'
#' The formula for the unbiased conversion is:
#'
#' \eqn{S = (t * \sqrt(2) * \Gamma(rdf/2)) / (\sqrt(n * rdf) * \Gamma((rdf - 1)/2))}
#'
#' The formula for the alternative estimator is derived by squaring the t statistic
#' and using the \code{\link{f2S}} formula. This estimator may be appealing for its
#' intuitive relationship to the F statistic; the absolute value of RESI estimates
#' using this formula will be equal to a RESI estimate using an F statistic for
#' the same model. However, this estimator does have finite sample bias, which is an
#' important consideration for the coverage of the bootstrapping that \code{resi} uses.
#'
#' The formula for the alternative conversion is:
#'
#' \eqn{ \sqrt(max(0, (t^2 * (rdf - 2)/rdf - 1)/rdf))}
#'
#' @return Returns a scalar or vector argument of the robust effect size index estimate.
#' @examples
#' # to obtain t values, first fit a lm
#' mod = lm(charges ~ region * age + bmi + sex, data = RESI::insurance)
#' # run lmtest::coeftest to get t values, using a robust variance-covariance formula
#' ts = lmtest::coeftest(mod, vcov. = sandwich::vcovHC)[,'t value']
#'
#' # get RESI estimates using unbiased estimator
#' t2S(ts, n = nrow(RESI::insurance), rdf = mod$df.residual)
#'
#' # get RESI estimates using alternative estimator
#' t2S(ts, n = nrow(RESI::insurance), rdf = mod$df.residual, unbiased = FALSE)
#' @export
t2S = function(t, rdf, n, unbiased = TRUE){
  if (unbiased){
    S = (t*sqrt(2))/(sqrt(n*rdf))*exp(lgamma(rdf/2) - lgamma((rdf-1)/2))
  } else{
    Ssq = (t^2*(rdf-2)/rdf - 1)/n
    S = sqrt(ifelse(Ssq<0, 0, Ssq))*(sqrt(t^2)/t)
  }
  return(S)
}
