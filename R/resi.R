
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
resi.lm <- function(model.full, model.reduced = NULL, data, anova = TRUE,
                    coefficients = TRUE, nboot = 1000, boot.method = 'nonparam',
                    vcovfunc = sandwich::vcovHC, alpha = 0.05, store.boot = FALSE,
                    Anova.args = list(), vcov.args = list(), unbiased = TRUE, ...){
  boot.method = match.arg(tolower(boot.method), choices = c("nonparam", "bayes"))

  if (missing(data)){
    data = model.full$model
    tryCatch(update(model.full, data = data), error = function(e){
      message("Updating model fit failed. Try rerunning with providing data argument")})
  } else{
    data = as.data.frame(data)
  }

  # point estimation
  output <- list(alpha = alpha, nboot = nboot, boot.method = tolower(boot.method))
  output = c(output, resi_pe(model.full = model.full, model.reduced = model.reduced,
                             data = data, anova = anova, coefficients = coefficients,
                             vcovfunc = vcovfunc, Anova.args = Anova.args,
                             vcov.args = vcov.args, unbiased = unbiased, ...))

  # bootstrapping
  boot.results = data.frame(matrix(nrow = nboot, ncol = length(output$estimates)))
  colnames(boot.results) = names(output$estimates)
  # non-parametric bootstrap
  if (tolower(boot.method) == "nonparam"){
    for (i in 1:nboot){
      boot.data = boot.samp(data)
      boot.model.full <- update(model.full, data = boot.data)
      if (is.null(model.reduced)){
        boot.model.reduced = NULL
      }
      else{
        boot.model.reduced = update(model.reduced, data = boot.data)
      }
      boot.results[i,] = suppressWarnings(resi_pe(model.full = boot.model.full, model.reduced = boot.model.reduced,
                                                  data = boot.data, anova = anova, coefficients = coefficients,
                                                  vcovfunc = vcovfunc, Anova.args = Anova.args, vcov.args = vcov.args,
                                                  unbiased = unbiased, ...)$estimates)
    }}

  # bayesian bootstrap
  if (tolower(boot.method)  == "bayes"){
    `(weights)` = NULL
    for (i in 1:nboot){
      boot.data = bayes.samp(data)
      boot.model.full <- suppressWarnings(update(model.full, data = boot.data, weights = boot.data[,'g']))
      if (is.null(model.reduced)){
        boot.model.reduced = suppressWarnings(update(model.full, formula = as.formula(paste(format(formula(model.full)[[2]]), "~ 1")), data = boot.model.full$model, weights = `(weights)`))
      }
      else{
        boot.model.reduced = suppressWarnings(update(model.reduced, data = boot.data, weights = boot.data[,'g']))
      }
      boot.results[i,] = suppressWarnings(resi_pe(model.full = boot.model.full, model.reduced = boot.model.reduced,
                                                  data = boot.data, anova = anova, coefficients = coefficients,
                                                  vcovfunc = vcovfunc, Anova.args = Anova.args, vcov.args = vcov.args,
                                                  unbiased = unbiased, ...)$estimates)
    }}

  alpha.order = sort(c(alpha/2, 1-alpha/2))
  output$overall[nrow(output$overall),paste(alpha.order*100, '%', sep='')] = quantile(boot.results[,1], probs = alpha.order, na.rm = TRUE)

  if (coefficients){
    CIs = apply(boot.results[,2:(1+nrow(output$coefficients))], 2,  quantile, probs = alpha.order, na.rm = TRUE)
    CIs = t(CIs)
    output$coefficients[1:nrow(CIs), paste(alpha.order*100, '%', sep='')] = CIs
  }

  if (anova){
    CIs = apply(boot.results[,(ncol(boot.results)-length(which(rownames(output$anova) != "Residuals"))+1):ncol(boot.results)], 2,  quantile, probs = alpha.order, na.rm = TRUE)
    CIs = t(CIs)
    output$anova[1:nrow(CIs), paste(alpha.order*100, '%', sep='')] = CIs
    class(output$anova) = c("anova_resi", class(output$anova))
  }

  if(store.boot){
    output$boot.results = boot.results
  }
  class(output) = "resi"
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
resi.geeglm <- function(object,
                        alpha = 0.05, nboot = 1000, anova = TRUE){
  output = resi_pe(object, anova = anova) # RESI point estimates
  data = object$data
  # id variable name
  id_var = as.character(object$call$id)
  # bootstrap
  output.boot = list(RESI = as.matrix(output$resi[, 'RESI']),
                     pm_RESI = as.matrix(output$resi[, 'pm-RESI']))
  corstr_spec = object$corstr
  for (i in 1:nboot){
    skip_to_next <- FALSE
    boot.data = boot.samp(data, id.var = id_var)
    # re-fit the model
    boot.mod = update(object, data = boot.data, corstr = corstr_spec)
    rv.boot = tryCatch(resi_pe(boot.mod, robust.var = robust.var, anova = anova), error = function(e) {skip_to_next <<- TRUE})
    if (skip_to_next) next
    output.boot$RESI= cbind(output.boot$RESI, rv.boot$resi[, 'RESI'])
    output.boot$pm_RESI = cbind(output.boot$pm_RESI, rv.boot$resi[, 'pm-RESI'])
  }
  output.boot$RESI = output.boot$RESI[, -1]
  output.boot$pm_RESI = output.boot$pm_RESI[, -1]
  RESI.ci = apply(output.boot$RESI, 1, quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
  RESI_se = apply(output.boot$RESI, 1, sd, na.rm = TRUE)
  rownames(RESI.ci) = paste("RESI", rownames(RESI.ci))
  pm_RESI.ci = apply(output.boot$pm_RESI, 1, quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
  pm_RESI_se = apply(output.boot$pm_RESI, 1, sd, na.rm = TRUE)
  rownames(pm_RESI.ci) = paste("pm-RESI", rownames(pm_RESI.ci))
  output = cbind(output$resi, t(RESI.ci), t(pm_RESI.ci), RESI_se, pm_RESI_se)
  cat("Note: the CI is actually based on", ncol(output.boot$RESI), "bootstraps. \n")
  cat("Function modified Jan 29 10:52pm. \n")
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







