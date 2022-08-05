#' Covert Cohen's \emph{d} to |S|
#'
#' Converts Cohen's \emph{d} robust effect size index (S) using the formula from
#' Vandekar, Rao, & Blume (2020).
#' @param d Numeric, value of Cohen's \emph{d}.
#' @param pi Numeric, the sampling proportions.
#' @return Returns an estimate the robust effect size index
#' @details The pi parameter comes from the fact that Cohen's d doesn't account
#' for unequal sample proportions in the population, but S does.
#'
#' The default is set to a natural value 1/2, which corresponds to a case
#' control design, for example, where sampling proportions always are
#' controlled by the experimenter.
#'
#' The formula to convert Cohen's \emph{d} to S is:
#'
#' \eqn{S = d/\sqrt( 1/\pi + 1/ (1 - \pi))}
#' @examples
#' # Consider an experiment with equal sampling proportions and a medium effect size
#' # corresponding to a Cohen's d of 0.5.
#' # convert to RESI (S)
#' d2S(d = 0.5)
#'
#' # This corresponds to a RESI of 0.25.
#' @export
d2S <- function(d, pi = 0.5){
  abs(d)/sqrt(1/pi + 1/(1-pi))
}
