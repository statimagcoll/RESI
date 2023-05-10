#' Compute the robust effect size index estimate from Z statistic
#'
#' This function computes the robust effect size index from Vandekar, Tao, & Blume (2020).
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' @param z The Z statistic for the parameter of interest.
#' @param n Number of independent samples.
#' @param unbiased Logical, whether to use unbiased or alternative estimator. See details.
#' @details This function computes S, the RESI, from a Z statistic. The formula for the
#' unbiased estimator (default) is derived by solving the expected value of the
#' Z statistic for S. It is unbiased and consistent.
#'
#' The formula for the unbiased conversion is:
#'
#' \eqn{S = Z/\sqrt(n)}
#'
#' The formula for the alternative estimator is derived by squaring the Z statistic
#' and using the \code{\link{chisq2S}} formula. This estimator may be appealing for its
#' intuitive relationship to the Chi-square statistic; the absolute value of RESI estimates
#' using this formula will be equal to a RESI estimate using a Chi-square statistic for
#' the same model. However, this estimator does have finite sample bias, which is an
#' important consideration for the coverage of the bootstrapping that \code{resi} uses.
#'
#' The formula for the alternative conversion is:
#'
#' \eqn{ \sqrt(max(0, (Z^2 - 1)/n)) * sign(Z)}
#' @return Returns a scalar or vector argument of the the robust effect size index estimate.
#' @examples
#' # to obtain example z values, first fit a glm
#' mod = glm(charges ~ region * age + bmi + sex, data = RESI::insurance)
#' # run coeftest to get z values using a robust variance-covariance function
#' zs = lmtest::coeftest(mod, vcov. = sandwich::vcovHC)[,'z value']
#'
#' # get RESI estimates using unbiased estimator
#' z2S(zs, n = nrow(RESI::insurance))
#'
#' # get RESI estimates usng alternative estimator
#' z2S(zs, n = nrow(RESI::insurance), unbiased = FALSE)
#' @export
z2S = function(z, n, unbiased = TRUE){
  if (unbiased){
    S = z/sqrt(n)
  } else{
    S = ifelse((z^2-1)<0, 0, sqrt((z^2-1)/n)*(sqrt(z^2)/z))
  }
  return(S)
}
