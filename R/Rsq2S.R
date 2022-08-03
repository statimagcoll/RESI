#' Covert R^2 to S
#'
#' Converts R^2, the partial coefficient of determination, to
#' robust effect size index (S) using the formula from Vandekar, Rao, & Blume (2020).
#' @param Rsq Numeric, R^2
#' @keywords power
#' @return Returns an estimate of R^2 based on the RESI
#' @examples
#' # consider a moderate effect size of R^2 = 0.1
#' Rsq2S(0.1)
#' # this corresponds to a RESI of 0.333
#' @export
Rsq2S <- function(Rsq){
  sqrt((-Rsq)/(Rsq-1))
}
