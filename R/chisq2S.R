#' Compute the robust effect size index estimate from chi-squared statistic.
#'
#' This function computes the robust effect size index from Vandekar, Tao, & Blume (2020).
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' For mixed effects models, RESI is conditional on the average correlation
#' structure within subjects.
#' @param chisq The chi-square statistic for the parameter of interest.
#' @param df Number of degrees of freedom of the chi-square statistic.
#' @param n Number of independent samples.
#' @return Returns a scalar or vector argument of the the robust effect size index estimate.
#' @details The formula for converting a Chi-square statistic to RESI is:
#'
#' \eqn{ S = \sqrt(max( 0, (chisq - df)/n))}
#' @examples
#' # obtain Chi-sq value by fitting an lm and running a Wald test
#' mod = lm(charges ~ region * age + bmi + sex, data = RESI::insurance)
#'
#' # run a Wald test with robust variance
#' wt = lmtest::waldtest(mod, vcov = sandwich::vcovHC, test = 'Chisq')
#'
#' # get Chi-sq value and degrees of freedom
#' chisq = wt$Chisq[2]
#' df = abs(wt$Df[2])
#'
#' # run chisq2S to convert to RESI
#' chisq2S(chisq, df = df, n = nrow(mod$model))
#'
#' @export
chisq2S <- function(chisq, df, n){
  S = (chisq - df)/n
  sqrt(ifelse(S<0, 0, S))
}
