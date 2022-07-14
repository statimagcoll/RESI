#' Robust Effect Size Index (RESI) point and interval estimation for models
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
#' @importFrom stats anova as.formula coef formula glm hatvalues nobs pchisq pf predict quantile rbinom residuals rnorm runif update vcov
#' @importFrom utils capture.output
#' @export
#' @details The RESI, denoted as S, is applicable across many model types. It is a unitless
#' index and can be easily be compared across models. The RESI can also be
#' converted to Cohen's \emph{d} (\code{\link{S2d}}) under  model homoskedasticity.
#'
#' This function computes the RESI point estimates and bootstrapped confidence
#' intervals based on Chi-square, F, T, or Z statistics. The robust (sandwich)
#' variance is used by default, allowing for consistency under
#' model-misspecification. The RESI is related to the non-centrality parameter
#' of the test statistic. The RESI estimate is consistent for all four
#' (Chi-square, F, T, and Z) types of statistics used. The Chi-square and F-based
#' calculations rely on asymptotic theory, so they may be biased in small samples.
#' When possible, the T and Z statistics are used, as these are unbiased and can
#' be positive or negative. The RESI based on the Chi-Square and F statistics is
#' always greater than or equal to 0. The type of statistic used is listed with
#' the output.
#'
#' For most model types supported by this package, two bootstrap options are
#' available. The default is the standard non-parametric bootstrap. Bayesian
#' bootstrapping is currently available for models other than \code{survreg} and
#' \code{coxph}.
#'
#' Certain model types require the data used for the model be entered as an argument.
#' These are: \code{nls, survreg,} and \code{coxph}. Additionally, if a model
#' includes certain functions (splines, factor, I), the data needs to be provided.
#'
#' If running into convergence issues with nls models, it is advised to refit the
#' nls model with starting values equal to the estimates provided by the model
#' and then try rerunning \code{resi}.
#'
#' @return Returns a list that includes function arguments, RESI point estimates,
#' and confidence intervals in summary/anova-style tables
#' @family RESI functions
#' @seealso \code{\link{resi_pe}}, \code{\link{vcovHC}}, \code{\link{nlshc}},
#' \code{\link{f2S}}, \code{\link{chisq2S}}, \code{\link{z2S}}, \code{\link{t2S}}
#' @references Vandekar S, Tao R, Blume J. A Robust Effect Size Index. \emph{Psychometrika}. 2020 Mar;85(1):232-246. doi: 10.1007/s11336-020-09698-2.
#'
#' Kang, K., Armstrong, K., Avery, S., McHugo, M., Heckers, S., & Vandekar, S. (2021). Accurate confidence interval estimation for non-centrality parameters and effect size indices. \emph{arXiv preprint arXiv:2111.05966}.

resi <- function(model.full, ...){
  UseMethod("resi")
}

#' @describeIn resi RESI point and interval estimation for models
#' @export
resi.default <- function(model.full, model.reduced = NULL, data, anova = TRUE, summary = TRUE,
                 nboot = 1000, boot.method = 'nonparam', vcovfunc = sandwich::vcovHC, alpha = 0.05, store.boot = FALSE, ...){
  boot.method = match.arg(tolower(boot.method), choices = c("nonparam", "bayes"))

  if (missing(data)){
    data = model.full$model
    tryCatch(update(model.full, data = data), error = function(e){
      message("Updating model fit failed. Try rerunning with providing data argument")})
  }
  # else{
  #   data = as.data.frame(data)
  # }

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
      boot.model.full <- suppressWarnings(update(model.full, data = boot.data, weights = boot.data[,'g']))
      if (is.null(model.reduced)){
        boot.model.reduced = suppressWarnings(update(model.full, formula = as.formula(paste(format(formula(model.full)[[2]]), "~ 1")), data = boot.model.full$model, weights = `(weights)`))
      }
      else{
        boot.model.reduced = suppressWarnings(update(model.reduced, data = boot.data, weights = boot.data[,'g']))
      }
      boot.results[i,] = suppressWarnings(resi_pe(model.full = boot.model.full, model.reduced = boot.model.reduced,
                                                  data = boot.data, anova = anova, summary = summary,
                                                  vcovfunc = vcovfunc, ...)$estimates)
    }}

  alpha.order = sort(c(alpha/2, 1-alpha/2))
  output$overall[nrow(output$overall),paste(alpha.order*100, '%', sep='')] = quantile(boot.results[,1], probs = alpha.order, na.rm = TRUE)

  if (summary){
    CIs = apply(boot.results[,2:(1+nrow(output$coefficients))], 2,  quantile, probs = alpha.order, na.rm = TRUE)
    CIs = t(CIs)
    output$coefficients[1:nrow(CIs), paste(alpha.order*100, '%', sep='')] = CIs
  }

  if (anova){
    CIs = apply(boot.results[,(ncol(boot.results)-length(which(rownames(output$anova) != "Residuals"))+1):ncol(boot.results)], 2,  quantile, probs = alpha.order, na.rm = TRUE)
    CIs = t(CIs)
    output$anova[1:nrow(CIs), paste(alpha.order*100, '%', sep='')] = CIs
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
                     nboot = 1000, boot.method = 'nonparam', vcovfunc = regtools::nlshc, alpha = 0.05, store.boot = FALSE, ...){
  boot.method = match.arg(tolower(boot.method), choices = c("nonparam", "bayes"))

  if (missing(data)){
    stop('\nData argument is required for nls model')
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

  alpha.order = sort(c(alpha/2, 1-alpha/2))
  output$overall[nrow(output$overall),paste(alpha.order*100, '%', sep='')] = quantile(boot.results[,1], probs = alpha.order, na.rm = TRUE)

  if (summary){
    CIs = apply(boot.results[,2:(1+nrow(output$coefficients))], 2,  quantile, probs = alpha.order, na.rm = TRUE)
    CIs = t(CIs)
    output$coefficients[1:nrow(CIs), paste(alpha.order*100, '%', sep='')] = CIs
  }

  if(store.boot){
    output$boot.results = boot.results
  }
  class(output) = "resi"
  return(output)
}

#' @describeIn resi RESI point and interval estimation for survreg models
#' @export
resi.survreg <- function(model.full, model.reduced = NULL, data, anova = TRUE, summary = TRUE,
                        nboot = 1000, boot.method = "nonparam", vcovfunc = vcov, alpha = 0.05, store.boot = FALSE, ...){
  if (missing(data)){
    stop('\nData argument is required for survreg model')
  }
  boot.method = match.arg(tolower(boot.method), choices = c("nonparam", "bayes"))
  if (boot.method == "bayes"){
    warning("Bayesian bootstrap not currently supported for survreg models, using non-parametric bootstrap")
  }

  resi.default(model.full = model.full, model.reduced = model.reduced, data = data,
               anova = anova, summary = summary, nboot = nboot, vcovfunc = vcovfunc,
               boot.method = "nonparam", store.boot = store.boot, ...)
}

#' @describeIn resi RESI point and interval estimation for coxph models
#' @export
resi.coxph <- function(model.full, model.reduced = NULL, data, anova = TRUE, summary = TRUE,
                         nboot = 1000, boot.method = "nonparam", vcovfunc = vcov, alpha = 0.05, store.boot = FALSE, ...){
  if (missing(data)){
    stop('\nData argument is required for coxph model')
  }
  boot.method = match.arg(tolower(boot.method), choices = c("nonparam", "bayes"))
  if (boot.method == "bayes"){
    warning("Bayesian bootstrap not currently supported for coxph models, using non-parametric bootstrap")
  }

  resi.default(model.full = model.full, model.reduced = model.reduced, data = data,
               anova = anova, summary = summary, nboot = nboot, vcovfunc = vcovfunc,
               boot.method = "nonparam", store.boot = store.boot, ...)
}

#' @describeIn resi RESI point and interval estimation for hurdle models
#' @export
resi.hurdle <- function(model.full, model.reduced = NULL, data, summary = TRUE,
                     nboot = 1000, boot.method = 'nonparam', vcovfunc = sandwich::sandwich, alpha = 0.05, store.boot = FALSE, ...){
  boot.method = match.arg(tolower(boot.method), choices = c("nonparam", "bayes"))

  if (missing(data)){
    data = model.full$model
    tryCatch(update(model.full, data = data), error = function(e){
      message("Updating model fit failed. Try rerunning with providing data argument")})
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


  output$overall[nrow(output$overall),c(paste(alpha/2*100, '%', sep=''), paste((1-alpha/2)*100, '%', sep=''))] = quantile(boot.results[,1], probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)

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

#' @describeIn resi RESI point and interval estimation for GEE models
#' @export
resi.geeglm <- function(model.full, alpha = 0.05, nboot = 1000, ...){
  warning("\nInterval performance not yet evaluated for geeglm")
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
  RESI.ci = apply(output.boot, 1, quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
  output$coefficients = cbind(output$coefficients, t(RESI.ci))
  output$boot.method = "nonparam"
  class(output)= 'resi'
  return(output)
}

#' @describeIn resi RESI point and interval estimation for GEE models
#' @export
resi.gee <- function(model.full, data, alpha = 0.05, nboot = 1000, ...){
  if (missing(data)){
    stop('\nData argument is required for GEE models from gee package')
  }
  warning("\nInterval performance not yet evaluated for gee")
  output <- list(alpha = alpha, nboot = nboot)
  # RESI point estimates
  output = c(output, resi_pe(model.full))
  # id variable name
  id_var = as.character(model.full$call$id)
  # bootstrap
  output.boot = as.matrix(output$coefficients[, 'RESI'])
  for (i in 1:nboot){
    boot.data = boot.samp(data, id.var = id_var)
    # re-fit the model
    suppressMessages(capture.output(boot.mod <- update(model.full, data = boot.data), file =  nullfile()))
    output.boot = cbind(output.boot, resi_pe(boot.mod)$coefficients[, 'RESI'])
  }
  output.boot = output.boot[, -1]
  RESI.ci = apply(output.boot, 1, quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
  output$coefficients = cbind(output$coefficients, t(RESI.ci))
  output$boot.method = "nonparam"
  class(output)= 'resi'
  return(output)
}

# convergence issues - no fix yet

#' @describeIn resi RESI point and interval estimation for LME (nlme) models
#' @importFrom nlme getGroups
#' @export
resi.lme <- function(model.full, alpha = 0.05, nboot = 1000, vcovfunc = clubSandwich::vcovCR, ...){
  warning("\nInterval performance not yet evaluated for lme")
  output <- list(alpha = alpha, nboot = nboot)
  # RESI point estimates
  output = c(output, resi_pe(model.full = model.full, vcovfunc = vcovfunc))
  data = model.full$data
  # id variable name
  id_var = attr(nlme::getGroups(model.full), "label")
  # bootstrap
  output.boot = as.matrix(output$coefficients[, 'RESI'])
  tryCatch(update(model.full, data = data), error = function(e){
    message("Need to run `library(nlme)`")})
  for (i in 1:nboot){
    boot.data = boot.samp(data, id.var = id_var)
    # re-fit the model
    boot.mod = update(model.full, data = boot.data, fixed = as.formula(model.full$call$fixed), random = as.formula(model.full$call$random))
    output.boot = cbind(output.boot, resi_pe(model.full = boot.mod, vcovfunc = vcovfunc)$coefficients[, 'RESI'])
  }
  output.boot = output.boot[, -1]
  RESI.ci = apply(output.boot, 1, quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
  output$coefficients = cbind(output$coefficients, t(RESI.ci))
  output$boot.method = "nonparam"
  class(output) = 'resi'
  return(output)
}

#' @export
resi.lmerMod <- function(model.full, alpha = 0.05, nboot = 1000, vcovfunc = clubSandwich::vcovCR, ...){
  warning("\nInterval performance not yet evaluated for lmerMod")
  output = resi_pe(model.full, vcovfunc = vcovfunc) # RESI point estimates
  data = model.full@frame
  # id variable name
  id_var = names(model.full@flist)
  # bootstrap (non-param for now)
  output.boot = as.matrix(output$coefficients[, 'RESI'])
  for (i in 1:nboot){
    boot.data = boot.samp(data, id.var = id_var)
    # re-fit the model
    boot.mod = update(model.full, data = boot.data)
    rv.boot = resi_pe(boot.mod, vcovfunc = vcovfunc)
    output.boot = cbind(output.boot, rv.boot$coefficients[, 'RESI'])
  }
  output.boot = output.boot[, -1]
  RESI.ci = apply(output.boot, 1, quantile, probs = c(alpha/2, 1-alpha/2))
  output = list(coefficients = cbind(output$coefficients, t(RESI.ci)), naive.var = output$naive.var)

  return(output)
}
