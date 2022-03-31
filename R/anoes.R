
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


#' Analysis of Effect Sizes (ANOES) based on the Robust Effect Size index (RESI) for (generalized) linear regression models\
#' @export
anoes.glm <- function(model.full, model.reduced = NULL,
                      robust.var = TRUE,
                      boot.type = 1, multi = 'none',
                      nboot = 1000, alpha = 0.05){

   output = boot.ci(model.full = model.full, model.reduced = model.reduced,
                     robust.var = robust.var,
                     boot.type = boot.type, multi = multi,
                     r = nboot,
                     alpha = alpha)$ANOES

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
  id_var = as.character(mod.gee$call$id)
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








