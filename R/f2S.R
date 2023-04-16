#' Compute the robust effect size index estimate from F-statistic
#'
#' This function computes the robust effect size index from Vandekar, Tao, & Blume (2020).
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' @param f The F statistic for the parameter of interest.
#' @param df Number of degrees of freedom of the F statistic.
#' @param rdf Model residual degrees of freedom.
#' @param n Number of independent samples.
#' @details The formula for converting an F statistic to S is:
#'
#' \eqn{ S = \sqrt(max(0, (f * df * (rdf - 2)/rdf - df)/n))}
#'
#' The estimator is derived by setting the statistic equal to the expected value of
#' the test statistic and solving for S.
#'
#' @return Returns a scalar or vector argument of the the robust effect size index estimate.
#' # to obtain example F values, first fit a lm
#' mod = lm(charges ~ region * age + bmi + sex, data = RESI::insurance)
#'
#' # run Anova, using a robust variance-covariance function
#' # get the F values and Df values
#' fs = car::Anova(mod, vcov. = sandwich::vcovHC)$F
#' dfs = car::Anova(mod, vcov. = sandwich::vcovHC)$Df
#'
#' # get RESI estimates
#' f2S(fs, df = dfs, rdf = mod$df.residual, n = nrow(RESI::insurance))
#' @export
f2S <- function(f, df, rdf, n){
  S = (f*df*(rdf-2)/rdf - df)/n
  sqrt(ifelse(S<0, 0, S))
}
