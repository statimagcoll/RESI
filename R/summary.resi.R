#' Summary method for resi objects
#'
#' After running the \code{\link{resi}} function on a fitted model, this function can be used to print the summary component. If the resi function was run with the `store.boot = TRUE` option to store the full matrix of bootstrapped estimates, the user can specify a different alpha level for this function's confidence intervals.
#' @param object an object resulting from resi function
#' @param alpha an optional new specification for the confidence level. Can be vector-valued
#' @return Returns a `summary.resi` object containing the computed summary table
#' @export
summary.resi <- function(object, alpha = NULL){
  if(is.null(object$coefficients)){
    stop('\nresi function was not run with summary = TRUE option')
  }

  output = list(alpha = alpha, model.full = object$model.full)
  if (is.null(alpha)){
    output$alpha = object$alpha
    output$coefficients = object$coefficients
  }
  else{
    if (is.null(object$boot.results)){
      stop('\nresi was not run with store.boot = TRUE option')
    }
    output$coefficients = object$coefficients[,1:(ncol(object$coefficients)-2)]
    CIs = apply(object$boot.results[,2:(1+nrow(object$coefficients))], 2,  quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
    CIs = t(CIs)
    output$coefficients[1:nrow(CIs), c(paste(alpha/2*100, '%', sep=''), paste((1-alpha/2)*100, '%', sep=''))] = CIs
  }
  class(output) = c('summary.resi')
  output
}

