#' Covert S to Cohen's \emph{f}^2
#'
#' Converts robust effect size index (S) to Cohen's \emph{f}^2
#' (effect size for multiple regression) using the formula from Vandekar, Rao, & Blume (2020).
#' @param S Numeric,the robust effect size index.
#' @details The formula for the conversion is:
#'
#' \eqn{ f^2 = S^2}
#' @return Returns an estimate of Cohen's \emph{f}^2 based on the RESI
#' @examples
#' # fit a linear regression model with continuous outcome and predictor
#' mod = lm(charges ~ age, data = RESI::insurance)
#'
#' # obtain t value for calculating RESI
#' t = summary(mod)$coefficients[2, "t value"]
#'
#' # calculate RESI
#' S = t2S(t, n = 1338, rdf = 1336)
#'
#' # convert to f^2
#' S2fsq(S)
#'
#' @export
S2fsq <- function(S){
  S^2
}
