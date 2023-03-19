#' Anova method for resi objects
#'
#' After running the \code{\link{resi}} function on a fitted model, this function can be used to print the Anova-style table component. If the resi function was run with the `store.boot = TRUE` option to store the full matrix of bootstrapped estimates, the user can specify a different alpha level for this function's confidence intervals.
#' @param mod an object resulting from resi function
#' @param alpha an optional new specification for the confidence level. Can be vector-valued
#' @param ... ignored
#' @return Returns an `anova` object containing the computed Anova-style table
#' @examples
#' # fit a model
#' mod = lm(charges ~ bmi + sex, data = RESI::insurance)
#'
#' # run resi with the store.boot = TRUE option
#' resi.obj = resi(mod, nboot = 100, store.boot = TRUE, alpha = 0.01)
#'
#' # run Anova, specifying a different alpha level if desired
#' car::Anova(resi.obj, alpha = 0.05)
#' @export
Anova.resi <- function(mod, alpha = NULL, ...){
  if(is.null(mod$anova)){
    stop('\nresi function was not run with anova = TRUE option')
  }

  if (is.null(alpha)){
    output = mod$anova
  }
  else{
    if (!(all(alpha %in% mod$alpha))){
      if (is.null(mod$boot.results)){
      stop('\nresi function was not run with store.boot = TRUE option')}
    }
    if(is.null(mod$boot.results)){
      output = mod$anova[c(1:(which(colnames(mod$anova) == 'RESI')),
                              which(colnames(mod$anova)%in%
                                      c(paste(alpha/2*100, '%', sep=''),
                                        paste((1-alpha/2)*100, '%', sep=''))))]
    }
    else{
      output = mod$anova[,1:(which(colnames(mod$anova) == 'RESI'))]
      CIs = apply(mod$boot.results[,(ncol(mod$boot.results)-
                                          nrow(mod$anova[which(rownames(mod$anova)
                                                                  != "Residuals"),])+1):
                                        ncol(mod$boot.results)], 2,  quantile,
                  probs = sort(c(alpha/2, 1-alpha/2)), na.rm = TRUE)
      CIs = t(CIs)
      output[1:nrow(CIs), c(paste(alpha/2*100, '%', sep=''),
                            paste((1-rev(alpha)/2)*100, '%', sep=''))] = CIs
    }
  }
  return(output)
}

#' Anova method for resi objects
#'
#' After running the \code{\link{resi}} function on a fitted model, this function can be used to print the Anova-style table component. If the resi function was run with the `store.boot = TRUE` option to store the full matrix of bootstrapped estimates, the user can specify a different alpha level for this function's confidence intervals.
#' @param object an object resulting from resi function
#' @param alpha an optional new specification for the confidence level. Can be vector-valued
#' @param ... ignored
#' @return Returns an `anova` object containing the computed Anova-style table
#' @details The resi function uses the car::Anova function to compute the Anova table.
#'@examples
#' # fit a model
#' mod = lm(charges ~ bmi + sex, data = RESI::insurance)
#'
#' # run resi with the store.boot = TRUE option
#' resi.obj = resi(mod, nboot = 100, store.boot = TRUE, alpha = 0.01)
#'
#' # run anova, specifying a different alpha level if desired
#' anova(resi.obj, alpha = 0.05)
#' @export
anova.resi <- function(object, alpha = NULL, ...){
  if(is.null(object$anova)){
    stop('\nresi function was not run with anova = TRUE option')
  }

  if (is.null(alpha)){
    output = object$anova
  }
  else{
    if (!(all(alpha %in% object$alpha))){
      if (is.null(object$boot.results)){
        stop('\nresi function was not run with store.boot = TRUE option')}
    }
    if(is.null(object$boot.results)){
      output = object$anova[c(1:(which(colnames(object$anova) == 'RESI')),
                              which(colnames(object$anova)%in%
                                      c(paste(alpha/2*100, '%', sep=''),
                                        paste((1-alpha/2)*100, '%', sep=''))))]
    }
    else{
      output = object$anova[,1:(which(colnames(object$anova) == 'RESI'))]
      CIs = apply(object$boot.results[,(ncol(object$boot.results)-
                                          nrow(object$anova)+1):
                                        ncol(object$boot.results)], 2,
                  quantile, probs = sort(c(alpha/2, 1-alpha/2)), na.rm = TRUE)
      CIs = t(CIs)
      output[1:nrow(CIs), c(paste(alpha/2*100, '%', sep=''),
                            paste((1-rev(alpha)/2)*100, '%', sep=''))] = CIs
    }
  }
  return(output)
}
