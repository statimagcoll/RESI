
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
resi <- function(x, ...){
  UseMethod("resi")
}


#' Analysis of Effect Sizes (ANOES) based on the Robust Effect Size index (RESI) for (generalized) linear regression models
#' This function estiamtes RESI and it CIs in a fitted glm model object
#' The CI are calculated via non-parametric bootstraps
#' @param model.full the full model. It should be a `glm` object.
#' @param model.reduced the reduced `glm()` model to compare with the full model. By default `NULL`, it's the same model as the full model but only having intercept.
#' @param nboot numeric, the number of bootstrap replicates. By default, 1000 bootstraps will be implemented.
#' @param robust.var default to TRUE, whether to use the robust (sandwich) variance estimator when construct the Wald test statistic. If `TRUE`, the variance of the estimator will be obtained by using `sandwich::vcovHC()`` and the HC3 will be applied.
#' @param boot.method which type of bootstrap to use: `nonparam` = non-parametric bootstrap (default); `bayes` = Bayesian bootstrap.
#' @param alpha significance level of the constructed CIs. By default, 0.05 will be used.
#' @export
resi.lm <- function(object, model.reduced = NULL,
                      robust.var = TRUE,
                      boot.method = "nonparam",
                      nboot = 1000, alpha = 0.05, ...){

  if (robust.var) {
    vcovfunc = sandwich::vcovHC
  } else {
    vcovfunc = vcov
  }

  if (! tolower(boot.method) %in% c("nonparam", "bayes")) stop("\n The bootstrap method should be either 'nonparam' for non-parametric bootstrap, or 'bayes' for Bayesian bootstrap")

  output = calc_resi(object, object.reduced = model.reduced, vcov. = vcovfunc, ...)$resi.tab
  data = object$model
  # bootstrap
  output.boot = as.matrix(output[, 'RESI'])
  ## non-param
  if (tolower(boot.method) == "nonparam"){
    for (i in 1:nboot){
      boot.data = boot.samp(data)
      # re-fit the model
      boot.mod = update(object, data = boot.data)
      if (is.null(model.reduced)){
        boot.mod.reduced = NULL
      } else {
        boot.mod.reduced = update(model.reduced, data = boot.data)
      }
      output.boot = cbind(output.boot, calc_resi(object = boot.mod, object.reduced = boot.mod.reduced, vcov. = vcovfunc, ...)$resi.tab[, 'RESI'])
    }
  }
  # bayesian boot
  if (tolower(boot.method)  == "bayes"){
    for (i in 1:nboot){
      boot.data = bayes.samp(data)
      # re-fit the model
      boot.mod = update(object, data = boot.data, weights = g)
      if (is.null(model.reduced)){
        boot.mod.reduced = NULL
      } else {
        boot.mod.reduced = update(model.reduced, data = boot.data, weights = g )
      }
      output.boot = cbind(output.boot, calc_resi(object = boot.mod, object.reduced = model.reduced, vcov. = vcovfunc, ...)$resi.tab[, 'RESI'])
    }
  } # end of `if (boot.method == "bayes")`

  output.boot = output.boot[, -1]
  RESI.ci = apply(output.boot, 1, quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
  output.tab = cbind(output, t(RESI.ci))
  if (robust.var) {
    colnames(output.tab)[which(colnames(output.tab) == 'Wald')] = 'Robust Wald'
    if (is.null(model.reduced)) colnames(output.tab)[which(colnames(output.tab) == 's.e.')] = 'Robust s.e.'
  }
  output = list(resi = output.tab,
                alpha = alpha,
                boot.value = output.boot, # bootstrapped values of RESI
                boot.method = tolower(boot.method),
                nboot = nboot,
                input.object = object,
                input.object.reduced = model.reduced,
                robust.var = robust.var
  )

  class(output) = "resi"
  return(output)
}


# resi.lm <- function(object, model.full = NULL, model.reduced = NULL,
#                     robust.var = TRUE,
#                     boot.method = "nonparam",
#                     nboot = 1000, alpha = 0.05){
#   if (is.null(model.reduced)){
#     output = calc_resi.glm(object)$resi.tab
#     data = object$model
#     # bootstrap
#     ## non-param
#     if (tolower(boot.method) == "nonparam"){
#       output.boot = as.matrix(output[, 'RESI'])
#       for (i in 1:nboot){
#         boot.data = boot.samp(data)
#         # re-fit the model
#         boot.mod = update(object, data = boot.data)
#         output.boot = cbind(output.boot, calc_resi.glm(boot.mod)$resi.tab[, 'RESI'])
#       }
#     }
#     # bayesian boot
#     if (tolower(boot.method)  == "bayes"){
#       output.boot = as.matrix(output[, 'RESI'])
#       for (i in 1:nboot){
#         boot.data = bayes.samp(data)
#         # re-fit the model
#         boot.mod = update(object, data = boot.data, weight = g)
#         output.boot = cbind(output.boot, calc_resi.glm(boot.mod)$resi.tab[, 'RESI'])
#       }
#     }
#
#     output.boot = output.boot[, -1]
#     RESI.ci = apply(output.boot, 1, quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
#     output.tab = cbind(output, t(RESI.ci))
#     output = list(resi = output.tab,
#                   alpha = alpha,
#                   boot.value = output.boot, # bootstrapped values of RESI
#                   boot.method = tolower(boot.method),
#                   nboot = nboot,
#                   input.object = object,
#                   robust.var = robust.var
#     )
#
#   } else{ # else: if there is a reduced model
#
#     output = boot.ci(model.full = model.full, model.reduced = model.reduced,
#                      robust.var = robust.var,
#                      boot.type = boot.type, multi = multi,
#                      r = nboot,
#                      alpha = alpha)$ANOES
#   }
#   class(output) = "resi"
#   return(output)
# }





#' Analysis of Effect Sizes (ANOES) based on the Robust Effect Size index (RESI) for GEE models
#' This function will estimate RESI and its CI for each factor in a fitted GEE model object.
#' The CIs are calculated via non-parametric bootstraps.
#' @param object the model object
#' @param alpha numeric, the type I error rate based on which CIs were constructed. By default, 0.05.
#' @param nboot numeric, the number of bootstraps used to construct CIs. By default, 1000.
#' @export
#' @return An ANOVA-type model summary output with RESI estimates and CIs added.
resi.geeglm <- function(object, robust.var = TRUE,
                        alpha = 0.05, nboot = 1000){
  output = calc_resi(object) # RESI point estimates
  data = object$data
  # id variable name
  id_var = as.character(object$call$id)
  # bootstrap
  output.boot = list(RESI = as.matrix(output[, 'RESI']),
                     pm_RESI = as.matrix(output[, 'pm-RESI']))
  corstr_spec = object$corstr
  for (i in 1:nboot){
    boot.data = boot.samp(data, id.var = id_var)
    # re-fit the model
    boot.mod = update(object, data = boot.data, corstr = corstr_spec)
    rv.boot = calc_resi(boot.mod, robust.var = robust.var)
    output.boot$RESI= cbind(output.boot$RESI, rv.boot[, 'RESI'])
    output.boot$pm_RESI = cbind(output.boot$pm_RESI, rv.boot[, 'pm-RESI'])
  }
  output.boot$RESI = output.boot$RESI[, -1]
  output.boot$pm_RESI = output.boot$pm_RESI[, -1]
  RESI.ci = apply(output.boot$RESI, 1, quantile, probs = c(alpha/2, 1-alpha/2))
  rownames(RESI.ci) = paste("RESI", rownames(RESI.ci))
  pm_RESI.ci = apply(output.boot$pm_RESI, 1, quantile, probs = c(alpha/2, 1-alpha/2))
  rownames(pm_RESI.ci) = paste("pm-RESI", rownames(pm_RESI.ci))
  output = cbind(output, t(RESI.ci), t(pm_RESI.ci))
  return(output)
}

# resi.gee <- function(object, data = NULL, alpha = 0.05, nboot = 1000){
#   output = calc_resi(object) # RESI point estimates
#
#   # NOTE: I couldn't find the data output from `gee` object.
#
#   # id variable name
#   id_var = as.character(object$call$id)
#   # bootstrap
#   output.boot = as.matrix(output[, 'RESI'])
#   for (i in 1:nboot){
#     boot.data = boot.samp(data, id.var = id_var)
#     # re-fit the model
#     boot.mod = update(object, data = boot.data)
#     output.boot = cbind(output.boot, calc_resi(boot.mod)[, 'RESI'])
#   }
#   output.boot = output.boot[, -1]
#   RESI.ci = apply(output.boot, 1, quantile, probs = c(alpha/2, 1-alpha/2))
#   output = cbind(output, t(RESI.ci))
#   return(output)
# }



#' Analysis of Effect Sizes (ANOES) based on the Robust Effect Size index (RESI) for LME models
#' This function will estimate RESI and its CI for each factor in a fitted GEE model object.
#' The CIs are calculated via non-parametric bootstraps.
#' @param object the model object
#' @param alpha numeric, the type I error rate based on which CIs were constructed. By default, 0.05.
#' @param nboot numeric, the number of bootstraps used to construct CIs. By default, 1000.
#' @export
#' @return An ANOVA-type model summary output with RESI estimates and CIs added.
resi.lme <- function(object, robust.var = TRUE,
                     alpha = 0.05, nboot = 1000){
  output = calc_resi(object) # RESI point estimates
  data = object$data
  # id variable name
  id_var = attr(nlme::getGroups(object), "label")
  # bootstrap (non-param for now)
  output.boot = list(RESI = as.matrix(output[, 'RESI']),
                     pm_RESI = as.matrix(output[, 'pm-RESI']))
  for (i in 1:nboot){
    boot.data = boot.samp(data, id.var = id_var)
    # re-fit the model
    boot.mod = update(object, data = boot.data)
    rv.boot = calc_resi(boot.mod, robust.var = robust.var)
    output.boot$RESI= cbind(output.boot$RESI, rv.boot[, 'RESI'])
    output.boot$pm_RESI = cbind(output.boot$pm_RESI, rv.boot[, 'pm-RESI'])
  }
  output.boot$RESI = output.boot$RESI[, -1]
  output.boot$pm_RESI = output.boot$pm_RESI[, -1]
  RESI.ci = apply(output.boot$RESI, 1, quantile, probs = c(alpha/2, 1-alpha/2))
  rownames(RESI.ci) = paste("RESI", rownames(RESI.ci))
  pm_RESI.ci = apply(output.boot$pm_RESI, 1, quantile, probs = c(alpha/2, 1-alpha/2))
  rownames(pm_RESI.ci) = paste("pm-RESI", rownames(pm_RESI.ci))
  output = cbind(output, t(RESI.ci), t(pm_RESI.ci))
  return(output)
}


#' Analysis of Effect Sizes (ANOES) based on the Robust Effect Size index (RESI) for LME models
#' This function will estimate RESI and its CI for each factor in a fitted LME model object using `lme4::lmer`.
#' The CIs are calculated via non-parametric bootstraps.
#' @param object the model object
#' @param alpha numeric, the type I error rate based on which CIs were constructed. By default, 0.05.
#' @param nboot numeric, the number of bootstraps used to construct CIs. By default, 1000.
#' @export
#' @return An ANOVA-type model summary output with RESI estimates and CIs added.
resi.lmerMod <- function(object, robust.var = TRUE, alpha = 0.05, nboot = 1000){
  output = calc_resi(object, robust.var = robust.var) # RESI point estimates
  data = object@frame
  # id variable name
  id_var = names(object@flist)
  # bootstrap (non-param for now)
  output.boot = list(RESI = as.matrix(output[, 'RESI']),
                     pm_RESI = as.matrix(output[, 'pm-RESI']))
  for (i in 1:nboot){
    boot.data = boot.samp(data, id.var = id_var)
    # re-fit the model
    boot.mod = update(object, data = boot.data)
    rv.boot = calc_resi(boot.mod, robust.var = robust.var)
    output.boot$RESI= cbind(output.boot$RESI, rv.boot[, 'RESI'])
    output.boot$pm_RESI = cbind(output.boot$pm_RESI, rv.boot[, 'pm-RESI'])
  }
  output.boot$RESI = output.boot$RESI[, -1]
  output.boot$pm_RESI = output.boot$pm_RESI[, -1]
  RESI.ci = apply(output.boot$RESI, 1, quantile, probs = c(alpha/2, 1-alpha/2))
  rownames(RESI.ci) = paste("RESI", rownames(RESI.ci))
  pm_RESI.ci = apply(output.boot$pm_RESI, 1, quantile, probs = c(alpha/2, 1-alpha/2))
  rownames(pm_RESI.ci) = paste("pm-RESI", rownames(pm_RESI.ci))
  output = cbind(output, t(RESI.ci), t(pm_RESI.ci))

  return(output)
}







