
#' Analysis of Effect Sizes (ANOES) based on the Robust Effect Size index (RESI) for (generalized) linear regression models
#' @export
# anoes <- function(model.full, model.reduced = NULL,
#                      robust.var = TRUE,
#                      boot.type = 1, multi = 'none',
#                      nboot = 1000, alpha = 0.05){
#
#   output = boot.ci.glm(model.full = model.full, model.reduced = model.reduced,
#                     robust.var = robust.var,
#                     boot.type = boot.type, multi = multi,
#                     r = nboot,
#                     alpha = alpha)$ANOES
#
#   return(output)
# }
#
anoes <- function(x, ...){
  UseMethod("anoes")
}


#' Analysis of Effect Sizes (ANOES) based on the Robust Effect Size index (RESI) for (generalized) linear regression models
#' This function estiamtes RESI and it CIs in a fitted glm model object
#' The CI are calculated via non-parametric bootstraps
#' @param model.full the full model. It should be a `glm` object.
#' @param model.reduced the reduced `glm()` model to compare with the full model. By default `NULL`, it's the same model as the full model but only having intercept.
#' @param nboot numeric, the number of bootstrap replicates. By default, 1000 bootstraps will be implemented.
#' @param robust.var default to TRUE, whether to use the robust (sandwich) variance estimator when construct the Wald test statistic. If `TRUE`, the variance of the estimator will be obtained by using `sandwich::vcovHC()`` and the HC3 will be applied.
#' @param multi the distribution from which the multipliers will be drawn: 'none' = the multipliers equal constant 1 (default); 'rad' = rademacher; 'normal' = Std Normal distribution
#' @param boot.type which type of bootstrap to use. 1: resampling covariates along with residuals (default); 2: fixing covariates and only bootstrapping residulas; 3: resampling covariates and residuals independently w/ replacements; 4. no sampling, just multipliers
#' @param alpha significance level of the constructed CIs. By default, 0.05 will be used0
#' @param digits the number of decimal digits in the output ANOES table. By default, 3
#' @export
anoes.glm <- function(object = NULL, model.full = NULL, model.reduced = NULL,
                      robust.var = TRUE,
                      boot.type = 1, multi = 'none',
                      nboot = 1000, alpha = 0.05){
  if (is.null(model.reduced)){
    output = resi.glm(object)
    data = object$model
    # bootstrap
    output.boot = as.matrix(output[, 'RESI'])
    for (i in 1:nboot){
      boot.data = boot.samp(data)
      # re-fit the model
      boot.mod = update(object, data = boot.data)
      output.boot = cbind(output.boot, resi.glm(boot.mod)[, 'RESI'])
    }
    output.boot = output.boot[, -1]
    RESI.ci = apply(output.boot, 1, quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
    output = cbind(output, t(RESI.ci))
    return(output)
  } else{
   output = boot.ci(model.full = model.full, model.reduced = model.reduced,
                     robust.var = robust.var,
                     boot.type = boot.type, multi = multi,
                     r = nboot,
                     alpha = alpha)$ANOES
  }
   return(output)
}


#' Analysis of Effect Sizes (ANOES) based on the Robust Effect Size index (RESI) for GEE models
#' This function will estimate RESI and its CI for each factor in a fitted GEE model object.
#' The CIs are calculated via non-parametric bootstraps.
#' @param object the model object
#' @param alpha numeric, the type I error rate based on which CIs were constructed. By default, 0.05.
#' @param nboot numeric, the number of bootstraps used to construct CIs. By default, 1000.
#' @export
#' @return An ANOVA-type model summary output with RESI estimates and CIs added.
anoes.geeglm <- function(object, alpha = 0.05, nboot = 1000){
  output = resi.geeglm(object) # RESI point estimates
  data = object$data
  # id variable name
  id_var = as.character(object$call$id)
  # bootstrap
  output.boot = as.matrix(output[, 'RESI'])
  for (i in 1:nboot){
    boot.data = boot.samp(data, id.var = id_var)
    # re-fit the model
    boot.mod = update(object, data = boot.data)
    output.boot = cbind(output.boot, resi.geeglm(boot.mod)[, 'RESI'])
  }
  output.boot = output.boot[, -1]
  RESI.ci = apply(output.boot, 1, quantile, probs = c(alpha/2, 1-alpha/2))
  output = cbind(output, t(RESI.ci))
  return(output)
}

#' Analysis of Effect Sizes (ANOES) based on the Robust Effect Size index (RESI) for LME models
#' This function will estimate RESI and its CI for each factor in a fitted GEE model object.
#' The CIs are calculated via non-parametric bootstraps.
#' @param object the model object
#' @param alpha numeric, the type I error rate based on which CIs were constructed. By default, 0.05.
#' @param nboot numeric, the number of bootstraps used to construct CIs. By default, 1000.
#' @export
#' @return An ANOVA-type model summary output with RESI estimates and CIs added.
anoes.lme <- function(object, alpha = 0.05, nboot = 1000){
  output = resi.lme(object) # RESI point estimates
  data = object$data
  # id variable name
  id_var = attr(nlme::getGroups(object), "label")
  # bootstrap
  output.boot = as.matrix(output[, 'RESI'])
  for (i in 1:nboot){
    boot.data = boot.samp(data, id.var = id_var)
    # re-fit the model
    boot.mod = update(object, data = boot.data)
    output.boot = cbind(output.boot, resi.lme(boot.mod)[, 'RESI'])
  }
  output.boot = output.boot[, -1]
  RESI.ci = apply(output.boot, 1, quantile, probs = c(alpha/2, 1-alpha/2))
  output = cbind(output, t(RESI.ci))
  return(output)
}








