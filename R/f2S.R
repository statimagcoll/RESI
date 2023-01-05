#' Compute the robust effect size index estimate from F-statistic
#'
#' This function computes the robust effect size index from Vandekar, Rao, & Blume (2020).
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' @param f The F statistic for the parameter of interest.
#' @param df Number of degrees of freedom of the F statistic.
#' @param rdf Model residual degrees of freedom or number of independent samples.
#' @details The formula for converting an F statistic to S is:
#'
#' \eqn{ S = \sqrt(max(0, (f * df * (rdf - 2)/rdf - 1)/rdf))}
#'
#' The estimator is derived by setting the statistic equal to the expected value of
#' the test statistic and solving for S. A modification of dividing by the residual degrees
#' of freedom rather than the total sample size is made based on the results of
#' simulations that showed this adjustment results in smaller finite-sample bias.
#'
#' @return Returns a scalar or vector argument of the the robust effect size index estimate.
#' # to obtain example F values, first fit a glm
#' mod = glm(charges ~ region * age + bmi + sex, data = RESI::insurance)
#'
#' # run Anova, using a robust variance-covariance function
#' # get the F values and Df values
#' fs = car::Anova(mod, vcov. = sandwich::vcovHC)$F
#' dfs = car::Anova(mod, vcov. = sandwich::vcovHC)$Df
#'
#' # get RESI estimates
#' f2S(fs, df = dfs, rdf = mod$df.residual)
#' @export
f2S <- function(f, df, rdf){
  S = (f*df*(rdf-2)/rdf - df)/rdf
  sqrt(ifelse(S<0, 0, S))
}
