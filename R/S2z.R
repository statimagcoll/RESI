#' Convert RESI (S) estimate to Z statistic
#'
#' Converts the robust effect size index (S) to Z statistic.
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' @param S The value of the RESI estimate.
#' @param n Number of independent samples.
#' @param unbiased Logical, whether the unbiased or alternative estimator was used to compute RESI estimate. Default is TRUE.
#' @return Returns a scalar or vector argument of the Chi-square statistic.
#' @details The formula for converting a RESI estimate to a corresponding Z statistic depends on which estimator
#' was used to compute the RESI estimate (unbiased vs. alternative, see \code{\link{z2S}}). For the unbiased estimator,
#' the RESI can be positive or negative and there is a 1-1 transformation from S to Z. The formula for converting S (unbiased) to the
#' Z statistic is:
#'
#' \eqn{\sqrt(n)*S}
#'
#' For the alternative formula, if the RESI estimate is 0, the Z statistic is only known within an interval, [-1, 1].
#' For a non-zero S, the formula is:
#'
#' \eqn{\sqrt{S^2}/S\sqrt(n*abs(S) + 1)}
#'
#' @examples
#' # convert S estimates with corresponding degrees of freedom to
#' Z statistics estimates (using unbiased formula)
#' S_ests = c(-0.2, 0, 0.1)
#' S2z(S = S_ests, n = 300, unbiased = TRUE)
#'
#' # convert S estimates with corresponding degrees of freedom to
#' Z statistics estimates (using alernative formula)
#' S_ests = c(-0.2, 0, 0.1)
#' S2z(S = S_ests, n = 300, unbiased = FALSE)
#' @export
S2z = function(S, n, unbiased = TRUE){
  if (unbiased){
    z = sqrt(n)*S
  } else{
    z = sqrt(S^2)/S * sqrt(n*abs(S) + 1)
    if (any(S == 0)){
      if(length(S) == 1){
        stop(paste0("Function is not 1-1 for S = 0, Z statistic is between -1 and 1"))
      } else{
        z[S == 0]= NA
        warning("Function is not 1-1 for S = 0, Z statistic is between -1 and 1")
      }
    }
  }
  return(z)
}
