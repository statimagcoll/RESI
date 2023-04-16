#' Compute the robust effect size index estimate from Z statistic
#'
#' This function computes the robust effect size index from Vandekar, Tao, & Blume (2020).
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' @param z The Z statistic for the parameter of interest.
#' @param n Number of independent samples.
#' @details This function computes S, the RESI, from a Z statistic using the default
#' formula. There is another function, \code{\link{z2S_alt}}, that uses an alternative
#' formula. This function's formula is derived by solving the expected value of the
#' Z statistic for S. It is unbiased and consistent.
#'
#' The formula for this conversion is:
#'
#' \eqn{ S = Z/\sqrt(n)}
#' @return Returns a scalar or vector argument of the the robust effect size index estimate.
#' @examples
#' # to obtain example z values, first fit a glm
#' mod = glm(charges ~ region * age + bmi + sex, data = RESI::insurance)
#' # run coeftest to get z values using a robust variance-covariance function
#' zs = lmtest::coeftest(mod, vcov. = sandwich::vcovHC)[,'z value']
#'
#' # get RESI estimates
#' z2S(zs, n = nrow(RESI::insurance))
#' @export
z2S <- function(z, n){
  z/sqrt(n)
}
