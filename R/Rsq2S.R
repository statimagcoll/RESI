#' Covert R^2 to S
#'
#' Converts R^2, the partial coefficient of determination, to
#' robust effect size index (S) using the formula from Vandekar, Tao, & Blume (2020).
#' @param Rsq Numeric, R^2
#' @return Returns an estimate of R^2 based on the RESI
#' @details The formula for the conversion is:
#'
#' \eqn{S = \sqrt((-R^2)/(R^2 - 1))}
#' @examples
#' # consider a moderate effect size of R^2 = 0.1
#' Rsq2S(0.1)
#' # this corresponds to a RESI of 0.333
#' @export
Rsq2S <- function(Rsq){
  if (any(Rsq < 0 | Rsq > 1)){
    stop("\nRsq must be between 0 and 1")
  }
  S = c()
  for (i in 1:length(Rsq)){
    if (Rsq[i] == 1){
      S[i] = Inf
    } else{
      S[i] = sqrt((-Rsq[i])/(Rsq[i]-1))
    }
  }
  return(S)
}
