#' Covert Cohen's \emph{f}^2 to S
#'
#' Converts Cohen's \emph{f}^2 to robust effect size index (S)
#' using the formula from Vandekar, Tao, & Blume (2020).
#' @param fsq Numeric, value of Cohen's \emph{f}^2.
#' @return Returns an estimate the robust effect size index
#' @details The formula for the conversion is:
#'
#' \eqn{S =  \sqrt(f^2)}
#' @examples
#' # consider a moderate effect size of f^2 = 0.3
#' fsq2S(0.3)
#' # This corresponds to a RESI of 0.5477226
#' @export
fsq2S <- function(fsq){
  if (any(fsq < 0)){
    stop("\nfsq must be non-negative")
  }
  sqrt(fsq)
}
