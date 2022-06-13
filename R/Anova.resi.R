#' Anova method for resi objects
#'
#' After running the \code{\link{resi}} function on a fitted model, this function can be used to print the Anova-style table component. If the resi function was run with the `store.boot = TRUE` option to store the full matrix of bootstrapped estimates, the user can specify a different alpha level for this function's confidence intervals.
#' @param resi.obj an object resulting from resi function
#' @param alpha an optional new specification for the confidence level. Can be vector-valued
#' @return Returns an `anova` object containing the computed Anova-style table
#' @export
Anova.resi <- function(resi.obj, alpha = NULL){
  if(is.null(resi.obj$anova)){
    stop('\n resi function was not run with anova = TRUE option')
  }

  if (is.null(alpha)){
    output = resi.obj$anova
  }
  else{
    output = resi.obj$anova[,1:(ncol(resi.obj$anova)-2)]
    CIs = apply(resi.obj$boot.results[,(ncol(resi.obj$boot.results)-nrow(resi.obj$anova)+1):ncol(resi.obj$boot.results)], 2,  quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
    CIs = t(CIs)
    output[1:nrow(CIs), c(paste(alpha/2*100, '%', sep=''), paste((1-alpha/2)*100, '%', sep=''))] = CIs
  }
  output
}
