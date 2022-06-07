#' Robust Effect Size index (RESI) point and interval estimation for models
#' @param model.full the full model.
#' @param model.reduced the reduced model to compare with the full model. By default `NULL`, it's the same model as the full model but only having intercept.
#' @param data optional data frame or object coercible to data frame of model.full data (required for some model types)
#' @param anova whether to produce an Anova table with the RESI columns added. By default = `TRUE`
#' @param summary whether to produce a summary table with the RESI columns added. By default = `TRUE`
#' @param nboot numeric, the number of bootstrap replicates. By default, 1000 bootstraps will be implemented.
#' @param boot.method which type of bootstrap to use: `nonparam` = non-parametric bootstrap (default); `bayes` = Bayesian bootstrap.
#' @param vcovfunc the variance estimator function for constructing the Wald test statistic. By default, sandwich::vcovHC (the robust (sandwich) variance estimator)
#' @param alpha significance level of the constructed CIs. By default, 0.05 will be used.
#' @param ... Other arguments to be passed to Anova function
#' @importFrom aod wald.test
#' @importFrom car Anova
#' @importFrom lmtest waldtest
#' @importFrom regtools nlshc
#' @importFrom sandwich vcovHC
#' @importFrom stats coef formula glm hatvalues pf predict quantile residuals update vcov
#' @export
#' @return Returns a list that includes function arguments, RESI point estimates, and confidence intervals

resi <- function(model.full, ...){
  UseMethod("resi")
}

#' @describeIn resi Robust Effect Size index (RESI) point and interval estimation for models
#' @param model.full the full model.
#' @param model.reduced the reduced model to compare with the full model. By default `NULL`, it's the same model as the full model but only having intercept.
#' @param data optional data frame or object coercible to data frame of model.full data (required for some model types)
#' @param anova whether to produce an Anova table with the RESI columns added. By default = `TRUE`
#' @param summary whether to produce a summary table with the RESI columns added. By default = `TRUE`
#' @param nboot numeric, the number of bootstrap replicates. By default, 1000 bootstraps will be implemented.
#' @param boot.method which type of bootstrap to use: `nonparam` = non-parametric bootstrap (default); `bayes` = Bayesian bootstrap.
#' @param vcovfunc the variance estimator function for constructing the Wald test statistic. By default, sandwich::vcovHC (the robust (sandwich) variance estimator)
#' @param alpha significance level of the constructed CIs. By default, 0.05 will be used.
#' @param ... Other arguments to be passed to Anova function
#' @export
resi.default <- function(model.full, model.reduced = NULL, data, anova = TRUE, summary = TRUE,
                 nboot = 1000, boot.method = 'nonparam', vcovfunc = sandwich::vcovHC, alpha = 0.05, ...){
  if (! tolower(boot.method) %in% c("nonparam", "bayes")) stop("\n The bootstrap method should be either 'nonparam' for non-parametric bootstrap, or 'bayes' for Bayesian bootstrap")

  if (missing(data)){
    data = model.full$model
  }

  # point estimation
  output <- list(model.full = formula(model.full), model.reduced = NULL, alpha = alpha,
                 `number of bootstraps` = nboot)
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
        boot.model.reduced = update(model.full, formula = as.formula(paste(format(formula(model.full)[[2]]), "~ 1")),  data = boot.data)
      }
      else{
        boot.model.reduced = update(model.reduced,  data = boot.data)
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
      if (is.null(model.reduced)){
        boot.model.reduced = update(model.full, formula = as.formula(paste(format(formula(model.full)[[2]]), "~ 1")),  data = boot.data, weights = g)
      }
      else{
        boot.model.reduced = update(model.reduced,  data = boot.data, weights = g)
      }
      boot.results[i,] = suppressWarnings(resi_pe(model.full = boot.model.full, model.reduced = boot.model.reduced,
                                                  data = boot.data, anova = anova, summary = summary,
                                                  vcovfunc = vcovfunc, ...)$estimates)

    }}


  output$overall[nrow(output$overall),c("LL", "UL")] = quantile(boot.results[,1], probs = c(alpha/2, 1-alpha/2))

  if (summary){
    CIs = apply(boot.results[,2:(1+nrow(output$summary))], 2,  quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
    CIs = t(CIs)
    output$summary[1:nrow(CIs), c("LL", "UL")] = CIs
  }

  if (anova){
    CIs = apply(boot.results[,(ncol(boot.results)-nrow(output$anova)+1):ncol(boot.results)], 2,  quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
    CIs = t(CIs)
    output$anova[1:nrow(CIs), c("LL", "UL")] = CIs
  }



  output$boot.results = boot.results
  output$model.reduced = formula(boot.model.reduced)
  class(output) = "resi"
  print(output[which(!(names(output) %in% c("boot.results", "estimates")))])
  return(invisible(output))
}

#' @describeIn resi Robust Effect Size index (RESI) point and interval estimation for nls models
#' @param model.full the full model.
#' @param model.reduced the reduced model to compare with the full model. By default `NULL`, it's the same model as the full model but only having intercept.
#' @param data optional data frame or object coercible to data frame of model.full data (required for some model types)
#' @param summary whether to produce a summary table with the RESI columns added. By default = `TRUE`
#' @param nboot numeric, the number of bootstrap replicates. By default, 1000 bootstraps will be implemented.
#' @param boot.method which type of bootstrap to use: `nonparam` = non-parametric bootstrap (default); `bayes` = Bayesian bootstrap.
#' @param vcovfunc the variance estimator function for constructing the Wald test statistic. By default, sandwich::vcovHC (the robust (sandwich) variance estimator)
#' @param alpha significance level of the constructed CIs. By default, 0.05 will be used.
#' @param ... Other arguments to be passed to Anova function
#' @export

resi.nls <- function(model.full, model.reduced = NULL, data, summary = TRUE,
                     nboot = 1000, boot.method = 'nonparam', vcovfunc = sandwich::vcovHC, alpha = 0.05, ...){
  if (! tolower(boot.method) %in% c("nonparam", "bayes")) stop("\n The bootstrap method should be either 'nonparam' for non-parametric bootstrap, or 'bayes' for Bayesian bootstrap")

  if (missing(data)){
    data = model.full$model
  }

  # point estimation
  output <- list(model.full = formula(model.full), model.reduced = NULL, alpha = alpha,
                 `number of bootstraps` = nboot)
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


  output$overall[nrow(output$overall),c("LL", "UL")] = quantile(boot.results[,1], probs = c(alpha/2, 1-alpha/2))

  if (summary){
    CIs = apply(boot.results[,2:(1+nrow(output$summary))], 2,  quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
    CIs = t(CIs)
    output$summary[1:nrow(CIs), c("LL", "UL")] = CIs
  }



  output$boot.results = boot.results
  class(output) = "resi"
  print(output[which(!(names(output) %in% c("boot.results", "estimates")))])
  return(invisible(output))
}
