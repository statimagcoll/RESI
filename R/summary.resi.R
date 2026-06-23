#' Summary method for resi objects
#'
#' After running the \code{\link{resi}} function on a fitted model, this function can
#' be used to print the coefficients table component. If the resi function was run with
#' the `store.boot = TRUE` option to store the full matrix of bootstrapped estimates,
#' the user can specify a different alpha level for this function's confidence intervals.
#' @param object an object resulting from resi function
#' @param alpha an optional new specification for the confidence level. Can be vector-valued
#' @param ... ignored
#' @return Returns a `summary_resi` object containing the computed coefficients table
#' @examples
#' # fit a model
#' mod = lm(charges ~ bmi + sex, data = RESI::insurance)
#'
#' # run resi with the store.boot = TRUE option
#' resi_obj = resi(mod, nboot = 100, store.boot = TRUE, alpha = 0.01)
#'
#' # run summary, specifying a different alpha level if desired
#' summary(resi_obj, alpha = 0.05)
#' @export
summary.resi <- function(object, alpha = NULL, ...){
  if(is.null(object$coefficients)){
    stop("\nresi function was not run with coefficients = TRUE option")
  }

  output = list(alpha = alpha, model.full = object$model.full)
  if (is.null(alpha)){
    output$alpha = object$alpha
    output$coefficients = object$coefficients
  }
  else{
    if (!(all(alpha %in% object$alpha))){
      # Scenario 1: model type does not support bootstrapping (e.g. lme),
      # OR store.boot was not used — in either case boot.results is NULL
      # and we cannot recompute quantiles at a different alpha.
      if (is.null(object$boot.results)){
        stop(paste0(
          "\nCannot compute CIs at a different alpha level. Either:\n",
          "  (a) bootstrapping is not yet supported for this model type, or\n",
          "  (b) resi() was not run with store.boot = TRUE.\n",
          "Re-run resi() with the desired alpha value directly."
        ))
      }
    }
    if(is.null(object$boot.results)){
      output$coefficients = object$coefficients[c(1:(which(colnames(object$coefficients)
                                                           == "RESI")),
                                                  which(colnames(object$coefficients)%in%
                                                          c(paste(alpha/2*100, "%", sep=""),
                                                            paste((1-alpha/2)*100, "%", sep=""))))]
    }
    else{
      output$coefficients = object$coefficients[,1:(which(colnames(object$coefficients) == "RESI"))]
      boot.results = object$boot.results$t
      CIs = apply(boot.results[,2:(1+nrow(object$coefficients))], 2,
                  quantile, probs = sort(c(alpha/2, 1-alpha/2)), na.rm = TRUE)
      CIs = t(CIs)
      output$coefficients[1:nrow(CIs), c(paste(alpha/2*100, "%", sep=""),
                                         paste((1-rev(alpha)/2)*100, "%", sep=""))] = CIs
    }
  }
  class(output) = c("summary_resi")
  output
}

