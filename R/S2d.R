#' Convert S to Cohen's d
#'
#' Converts the robust effect size index (S) to Cohen's d using the formula from
#' Vandekar, Tao, & Blume (2020).
#' @param S Numeric, the robust effect size index.
#' @param pi Numeric, the sampling proportions.
#' @return Returns an estimate of Cohen's \emph{d} based on the RESI.
#' @details The pi parameter comes from the fact that Cohen's d doesn't account
#' for unequal sample proportions in the population, but S does.
#'
#' The default is set to a natural value 1/2, which corresponds to a case
#' control design, for example, where sampling proportions always are
#' controlled by the experimenter.
#'
#' The formula for the conversion is:
#'
#' \eqn{ d = | S * \sqrt(1/\pi + 1/(1 - \pi)) |}
#' @examples
#' # fit a simple linear regression with a binary predictor
#' mod = lm(charges ~ sex, data = RESI::insurance)
#'
#' # calculate t-value
#' t = summary(mod)$coefficients[2, "t value"]
#'
#' # calculate RESI (S)
#' S = t2S(t, n = 1338, rdf = 1336)
#'
#' # determine sample proportions
#' pi = length(which(RESI::insurance[,"sex"]=="male"))/1338
#'
#' # convert S to Cohen's d
#' S2d(S = S, pi = pi)
#'
#' @export
S2d = function(S, pi=0.5){
  abs(S * sqrt(1/pi + 1/(1-pi)))
}
