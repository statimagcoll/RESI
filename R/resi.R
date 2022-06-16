#' Robust Effect Size index (RESI) point and interval estimation for models
#'
#' This function will estimate the robust effect size (RESI) from Vandekar, Rao, & Blume (2020) and its confidence interval in various ways for a fitted model object. The overall RESI is estimated via a Wald test. RESI is (optionally) estimated for each factor in summary-style table. RESI is (optionally) estimated for each variable/interaction in an Anova-style table for models with existing Anova methods. CIs can be calculated using either non-parametric or Bayesian bootstrapping.
#' @param model.full \code{lm, glm, nls, survreg, coxph, hurdle, zeroinfl, gee, geeglm} or \code{lme} model object.
#' @param model.reduced Fitted model object of same type as model.full. By default `NULL`; the same model as the full model but only having intercept.
#' @param data Data.frame or object coercible to data.frame of model.full data (required for some model types).
#' @param vcovfunc The variance estimator function for constructing the Wald test statistic. By default, sandwich::vcovHC (the robust (sandwich) variance estimator).
#' @param summary Logical, whether to produce a summary (coefficients) table with the RESI columns added. By default = `TRUE`.
#' @param anova Logical, whether to produce an Anova table with the RESI columns added. By default = `TRUE`.
#' @param nboot Numeric, the number of bootstrap replicates. By default, 1000.
#' @param boot.method String, which type of bootstrap to use: `nonparam` = non-parametric bootstrap (default); `bayes` = Bayesian bootstrap.
#' @param alpha Numeric, significance level of the constructed CIs. By default, 0.05.
#' @param store.boot Logical, whether to store all the bootstrapped estimates. By default, `FALSE`.
#' @param ... Other arguments to be passed to Anova function.
#' @importFrom aod wald.test
#' @importFrom car Anova
#' @importFrom lmtest waldtest
#' @importFrom regtools nlshc
#' @importFrom sandwich vcovHC
#' @importFrom stats coef formula glm hatvalues pf predict quantile residuals update vcov
#' @export
#' @return Returns a list that includes function arguments, RESI point estimates, and confidence intervals in summary/anova-style tables

resi <- function(model.full, ...){
  UseMethod("resi")
}

#' @describeIn resi RESI point and interval estimation for models
#' @export
resi.default <- function(model.full, model.reduced = NULL, data, anova = TRUE, summary = TRUE,
                 nboot = 1000, boot.method = 'nonparam', vcovfunc = sandwich::vcovHC, alpha = 0.05, store.boot = FALSE, ...){
  if (! tolower(boot.method) %in% c("nonparam", "bayes")) stop("\n The bootstrap method should be either 'nonparam' for non-parametric bootstrap, or 'bayes' for Bayesian bootstrap")

  if (missing(data)){
    data = model.full$model
  }

  # point estimation
  output <- list(alpha = alpha, nboot = nboot, boot.method = tolower(boot.method))
  output = c(output, resi_pe(model.full = model.full, model.reduced = model.reduced, data = data, anova = anova, summary = summary, vcovfunc = vcovfunc, ...))

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
                                                  data = boot.data, anova = anova, summary = summary,
                                                  vcovfunc = vcovfunc, ...)$estimates)
    }}

  # bayesian bootstrap
  if (tolower(boot.method)  == "bayes"){
    for (i in 1:nboot){
      boot.data = bayes.samp(data)
      boot.model.full <- update(model.full, data = boot.data, weights = g)
      boot.model.reduced = update(model.reduced,  data = boot.data, weights = g)
      boot.results[i,] = suppressWarnings(resi_pe(model.full = boot.model.full, model.reduced = boot.model.reduced,
                                                  data = boot.data, anova = anova, summary = summary,
                                                  vcovfunc = vcovfunc, ...)$estimates)
    }}


  output$overall[nrow(output$overall),c(paste(alpha/2*100, '%', sep=''), paste((1-alpha/2)*100, '%', sep=''))] = quantile(boot.results[,1], probs = c(alpha/2, 1-alpha/2))

  if (summary){
    CIs = apply(boot.results[,2:(1+nrow(output$coefficients))], 2,  quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
    CIs = t(CIs)
    output$coefficients[1:nrow(CIs), c(paste(alpha/2*100, '%', sep=''), paste((1-alpha/2)*100, '%', sep=''))] = CIs
  }

  if (anova){
    CIs = apply(boot.results[,(ncol(boot.results)-nrow(output$anova)+1):ncol(boot.results)], 2,  quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
    CIs = t(CIs)
    output$anova[1:nrow(CIs), c(paste(alpha/2*100, '%', sep=''), paste((1-alpha/2)*100, '%', sep=''))] = CIs
  }

  if(store.boot){
    output$boot.results = boot.results
  }
  class(output) = "resi"
  return(output)
}

#' @describeIn resi RESI point and interval estimation for nls models
#' @export

resi.nls <- function(model.full, model.reduced = NULL, data, summary = TRUE,
                     nboot = 1000, boot.method = 'nonparam', vcovfunc = sandwich::vcovHC, alpha = 0.05, store.boot = FALSE, ...){
  if (! tolower(boot.method) %in% c("nonparam", "bayes")) stop("\n The bootstrap method should be either 'nonparam' for non-parametric bootstrap, or 'bayes' for Bayesian bootstrap")

  if (missing(data)){
    data = model.full$model
  }

  # point estimation
  output <- list(alpha = alpha, nboot = nboot, boot.method = tolower(boot.method))
  output = c(output, resi_pe(model.full = model.full, model.reduced = model.reduced, data = data, anova = anova, summary = summary, vcovfunc = vcovfunc, ...))

  # bootstrapping
  boot.results = data.frame(matrix(nrow = nboot, ncol = length(output$estimates)))
  colnames(boot.results) = names(output$estimates)
  # non-parametric bootstrap
  if (tolower(boot.method) == "nonparam"){
    for (i in 1:nboot){
      boot.data = boot.samp(data)
      boot.model.full <- update(model.full, data = boot.data)
      boot.results[i,] = suppressWarnings(resi_pe(model.full = boot.model.full, model.reduced = NULL,
                                                  data = boot.data, summary = summary,
                                                  vcovfunc = vcovfunc, ...)$estimates)

    }}

  # bayesian bootstrap
  if (tolower(boot.method)  == "bayes"){
    for (i in 1:nboot){
      boot.data = bayes.samp(data)
      boot.model.full <- update(model.full, data = boot.data, weights = g)
      boot.results[i,] = suppressWarnings(resi_pe(model.full = boot.model.full, model.reduced = NULL,
                                                  data = boot.data, summary = summary,
                                                  vcovfunc = vcovfunc, ...)$estimates)

    }}


  output$overall[nrow(output$overall),c(paste(alpha/2*100, '%', sep=''), paste((1-alpha/2)*100, '%', sep=''))] = quantile(boot.results[,1], probs = c(alpha/2, 1-alpha/2))

  if (summary){
    CIs = apply(boot.results[,2:(1+nrow(output$coefficients))], 2,  quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
    CIs = t(CIs)
    output$coefficients[1:nrow(CIs), c(paste(alpha/2*100, '%', sep=''), paste((1-alpha/2)*100, '%', sep=''))] = CIs
  }

  if(store.boot){
    output$boot.results = boot.results
  }
  class(output) = "resi"
  return(output)
}

#' @describeIn resi RESI point and interval estimation for hurdle models
#' @export
resi.hurdle <- function(model.full, model.reduced = NULL, data, summary = TRUE,
                     nboot = 1000, boot.method = 'nonparam', vcovfunc = sandwich::sandwich, alpha = 0.05, store.boot = FALSE, ...){
  if (! tolower(boot.method) %in% c("nonparam", "bayes")) stop("\n The bootstrap method should be either 'nonparam' for non-parametric bootstrap, or 'bayes' for Bayesian bootstrap")

  if (missing(data)){
    data = model.full$model
  }

  if (is.null(model.reduced)){
    model.reduced = update(model.full, formula = as.formula(paste(format(formula(model.full)[[2]]), "~ 1")),  data = model.full$model)
  }

  # point estimation
  output <- list(alpha = alpha, nboot = nboot, boot.method = tolower(boot.method))
  output = c(output, resi_pe(model.full = model.full, model.reduced = model.reduced, data = data, anova = anova, summary = summary, vcovfunc = vcovfunc, ...))

  # bootstrapping
  boot.results = data.frame(matrix(nrow = nboot, ncol = length(output$estimates)))
  colnames(boot.results) = names(output$estimates)
  # non-parametric bootstrap
  if (tolower(boot.method) == "nonparam"){
    for (i in 1:nboot){
      boot.data = boot.samp(data)
      boot.model.full <- update(model.full, data = boot.data)
      boot.model.reduced <- update(model.full, data = boot.data)
      boot.results[i,] = suppressWarnings(resi_pe(model.full = boot.model.full, model.reduced = NULL,
                                                  data = boot.data, summary = summary,
                                                  vcovfunc = vcovfunc, ...)$estimates)

    }}

  # bayesian bootstrap
  if (tolower(boot.method)  == "bayes"){
    for (i in 1:nboot){
      boot.data = suppressMessages(bayes.samp(data))
      boot.model.full <- suppressMessages(update(model.full, data = boot.data, weights = g))
      boot.model.reduced <- suppressMessages(update(model.reduced, data = boot.data, weights = g))
      boot.results[i,] = suppressMessages(resi_pe(model.full = boot.model.full, model.reduced = NULL,
                                                  data = boot.data, summary = summary,
                                                  vcovfunc = vcovfunc, ...)$estimates)
    }}


  output$overall[nrow(output$overall),c(paste(alpha/2*100, '%', sep=''), paste((1-alpha/2)*100, '%', sep=''))] = quantile(boot.results[,1], probs = c(alpha/2, 1-alpha/2))

  if (summary){
    CIs = apply(boot.results[,2:(1+nrow(output$coefficients))], 2,  quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
    CIs = t(CIs)
    output$coefficients[1:nrow(CIs), c(paste(alpha/2*100, '%', sep=''), paste((1-alpha/2)*100, '%', sep=''))] = CIs
  }

  if(store.boot){
    output$boot.results = boot.results
  }
  class(output) = "resi"
  return(output)
}

#' @describeIn resi RESI point and interval estimation for zeroinfl models
#' @export
resi.zeroinfl <- resi.hurdle


# need to finish making output format consistent
#' @describeIn resi RESI point and interval estimation for GEE models
#' @export
resi.geeglm <- function(model.full, alpha = 0.05, nboot = 1000){
  #browser()
  output <- list(alpha = alpha, nboot = nboot)
  # RESI point estimates
  output = c(output, resi_pe(model.full))
  data = model.full$data
  # id variable name
  id_var = as.character(model.full$call$id)
  # bootstrap
  output.boot = as.matrix(output$coefficients[, 'RESI'])
  for (i in 1:nboot){
    boot.data = boot.samp(data, id.var = id_var)
    # re-fit the model
    boot.mod = update(model.full, data = boot.data)
    output.boot = cbind(output.boot, resi_pe(boot.mod)$coefficients[, 'RESI'])
  }
  output.boot = output.boot[, -1]
  RESI.ci = apply(output.boot, 1, quantile, probs = c(alpha/2, 1-alpha/2))
  output = cbind(output$coefficients, t(RESI.ci))
  class(output)= 'resi'
  return(output)
}

#' @describeIn resi RESI point and interval estimation for LME models
#' @importFrom nlme getGroups
#' @export
resi.lme <- function(model.full, alpha = 0.05, nboot = 1000){
  output = resi_pe(model.full) # RESI point estimates
  data = model.full$data
  # id variable name
  id_var = attr(nlme::getGroups(model.full), "label")
  # bootstrap
  output.boot = as.matrix(output[, 'RESI'])
  for (i in 1:nboot){
    boot.data = boot.samp(data, id.var = id_var)
    # re-fit the model
    boot.mod = update(model.full, data = boot.data)
    output.boot = cbind(output.boot, resi_pe(boot.mod)[, 'RESI'])
  }
  output.boot = output.boot[, -1]
  RESI.ci = apply(output.boot, 1, quantile, probs = c(alpha/2, 1-alpha/2))
  output = cbind(output, t(RESI.ci))
  class(output) = 'resi'
  return(output)
}
