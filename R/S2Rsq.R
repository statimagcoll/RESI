#' Covert S to R^2
#'
#' Converts robust effect size index (S) to R^2, the partial
#' coefficient of determination, using the formula from Vandekar, Tao, & Blume (2020).
#' @param S Numeric, the robust effect size index.
#' @details The formula for the conversion is:
#'
#' \eqn{ R^2 = S^2 / (1 + S^2)}
#' @return Returns an estimate of R^2 based on the RESI
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
#' # convert S to R^2
#' S2Rsq(S)
#'
#' @export
S2Rsq <- function(S){
  S^2/(1+S^2)
}
