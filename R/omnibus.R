#' Omnibus (Overall) Wald Test for resi objects
#'
#' After running the \code{\link{resi}} function on a fitted model, this function
#' can be used to print the overall Wald test component. If the resi function was run
#' with the `store.boot = TRUE` option to store the full matrix of bootstrapped estimates,
#' the user can specify a different alpha level for this function's confidence intervals.
#' @param object an object resulting from resi function
#' @param alpha an optional new specification for the confidence level. Can be vector-valued
#' @param ... ignored
#' @return Returns a `omnibus_resi` object containing the computed omnibus Wald test
#' @examples
#' # fit a model
#' mod = lm(charges ~ bmi + sex, data = RESI::insurance)
#'
#' # run resi with the store.boot = TRUE option
#' resi_obj = resi(mod, nboot = 100, store.boot = TRUE, alpha = 0.01)
#'
#' # run summary, specifying a different alpha level if desired
#' omnibus(resi_obj, alpha = 0.05)
#' @export
omnibus <- function(object, alpha = NULL, ...){
  if(!"resi" %in% class(object)){
    stop("\nThis function only works on resi objects")
  }
  if(is.null(object$overall)){
    stop("\nThe resi object did not compute an omnibus test")
  }

  output = list(alpha = alpha, model.full = object$model.full)
  if (is.null(alpha)){
    output$alpha = object$alpha
    output$overall = object$overall
  }
  else{
    if (!(all(alpha %in% object$alpha))){
      if (is.null(object$boot.results)){
        stop("\nresi function was not run with store.boot = TRUE option")}
    }
    if(is.null(object$boot.results)){
      output$overall = object$overall[,c(1:(which(colnames(object$overall)
                                                           == "RESI")),
                                                  which(colnames(object$overall)%in%
                                                          c(paste(alpha/2*100, "%", sep=""),
                                                            paste((1-alpha/2)*100, "%", sep=""))))]
    }
    else{
      output$overall = object$overall[,1:(which(colnames(object$overall) == "RESI"))]
      boot.results = object$boot.results$t
      CIs = quantile(boot.results[,1], probs = sort(c(alpha/2, 1-alpha/2)),
                     na.rm = TRUE)
      output$overall[nrow(object$overall), c(paste(alpha/2*100, "%", sep=""),
                                         paste((1-rev(alpha)/2)*100, "%", sep=""))] = CIs
    }
  }
  class(output) = c("omnibus_resi")
  output
}

