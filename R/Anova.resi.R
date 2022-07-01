#' Anova method for resi objects
#'
#' After running the \code{\link{resi}} function on a fitted model, this function can be used to print the Anova-style table component. If the resi function was run with the `store.boot = TRUE` option to store the full matrix of bootstrapped estimates, the user can specify a different alpha level for this function's confidence intervals.
#' @param object an object resulting from resi function
#' @param alpha an optional new specification for the confidence level. Can be vector-valued
#' @param ... ignored
#' @return Returns an `anova` object containing the computed Anova-style table
#' @export
Anova.resi <- function(object, alpha = NULL, ...){
  if(is.null(object$anova)){
    stop('\nresi function was not run with anova = TRUE option')
  }

  if (is.null(alpha)){
    output = object$anova
  }
  else{
    if (is.null(object$boot.results)){
      stop('\nresi function was not run with store.boot = TRUE option')
    }
    output = object$anova[,1:(ncol(object$anova)-2)]
    CIs = apply(object$boot.results[,(ncol(object$boot.results)-nrow(object$anova)+1):ncol(object$boot.results)], 2,  quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
    CIs = t(CIs)
    output[1:nrow(CIs), c(paste(alpha/2*100, '%', sep=''), paste((1-alpha/2)*100, '%', sep=''))] = CIs
  }
  class(output) = c('anova.resi',class(output))
  output
}