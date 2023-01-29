#' Robust Effect Size Index (RESI) point and interval estimation for models
#'
#' This function will estimate the robust effect size (RESI) from Vandekar, Rao, & Blume (2020) and its confidence interval in various ways for a fitted model object. The overall RESI is estimated via a Wald test. RESI is (optionally) estimated for each factor in coefficients-style table. RESI is (optionally) estimated for each variable/interaction in an Anova-style table for models with existing Anova methods. CIs can be calculated using either non-parametric or Bayesian bootstrapping.
#' @param model.full \code{lm, glm, nls, survreg, coxph, hurdle, zeroinfl, gee, geeglm} or \code{lme} model object.
#' @param model.reduced Fitted model object of same type as model.full. By default `NULL`; the same model as the full model but only having intercept.
#' @param data Data.frame or object coercible to data.frame of model.full data (required for some model types).
#' @param vcovfunc The variance estimator function for constructing the Wald test statistic. By default, sandwich::vcovHC (the robust (sandwich) variance estimator).
#' @param coefficients Logical, whether to produce a coefficients (summary) table with the RESI columns added. By default = `TRUE`.
#' @param anova Logical, whether to produce an Anova table with the RESI columns added. By default = `TRUE`.
#' @param nboot Numeric, the number of bootstrap replicates. By default, 1000.
#' @param boot.method String, which type of bootstrap to use: `nonparam` = non-parametric bootstrap (default); `bayes` = Bayesian bootstrap.
#' @param alpha Numeric, significance level of the constructed CIs. By default, 0.05.
#' @param store.boot Logical, whether to store all the bootstrapped estimates. By default, `FALSE`.
#' @param Anova.args List, additional arguments to be passed to \link[car]{Anova} function.
#' @param vcov.args List, additional arguments to be passed to vcovfunc.
#' @param unbiased Logical, whether to use the unbiased or alternative T/Z statistic to RESI conversion. By default, `TRUE`. See details.
#' @param ... Ignored.
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
#' When possible, the T and Z statistics are used. There are two formulas for both
#' the T and Z statistic conversion. The first (default, unbiased = TRUE)
#' are based on solving the expected value of the T or Z statistic for the RESI.
#' The alternative is based on squaring the T or Z statistic and using the
#' F or Chi-square statistic conversion. Both of these methods are consistent, but
#' the alternative exhibits a notable amount of finite sample bias. The alternative
#' may be appealing because its absolute value will be equal to the RESI based on
#' the F or Chi-square statistic. The RESI based on the Chi-Square and F statistics
#' is always greater than or equal to 0. The type of statistic
#' used is listed with the output. See \code{\link{f2S}}, \code{\link{chisq2S}},
#' \code{\link{t2S}}, \code{\link{z2S}}, \code{\link{t2S_alt}}, and
#' \code{\link{z2S_alt}} for more details on the formulas.
#'
#' For most \code{lm} and \code{nls} model types, there is a Bayesian bootstrap
#' option available as an alternative to the default, standard non-parametric
#' bootstrap. The interpretation of a Bayesian bootstrapped interval is similar to
#' that of a credible interval.
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
#' and confidence intervals in coefficients/anova-style tables
#' @family RESI functions
#' @seealso \code{\link{resi_pe}}, \code{\link{vcovHC}}, \code{\link{nlshc}},
#' \code{\link{f2S}}, \code{\link{chisq2S}}, \code{\link{z2S}}, \code{\link{t2S}}
#'
#' @examples
#' ## for timing purposes, a small number of bootstrap replicates is used in the
#' ## examples. Run them with a higher or default `nboot` argument for better performance
#'
#' ## RESI on a linear model
#' # fit linear model
#' mod = lm(charges ~ region * age + bmi + sex, data = RESI::insurance)
#'
#' # run resi on fitted model with desired number of bootstrap replicates
#' resi(mod, nboot = 50)
#'
#' # fit a reduced model for comparison
#' mod.red = lm(charges ~ bmi, data = RESI::insurance)
#'
#' # running resi and including the reduced model will provide almost the exact same
#' # output as not including a reduced model. The difference is that the "overall"
#' # portion of the output will compare the full model to the reduced model.
#' # The "summary" and "anova" RESI estimates will be the same. (The bootstrapped
#' # confidence intervals may differ.)
#' resi(model.full = mod, model.reduced = mod.red, nboot = 10)
#'
#' # for some model types and formula structures, data argument is required
#' library(splines)
#' # fit logistic regression model with splines
#' mod = glm(smoker ~ ns(age, df = 3) + region, data = RESI::insurance, family = "binomial")
#'
#' # store bootstrap results for calculating different CIs later
#' # specify additional arguments to the variance-covariance function via vcov.args
#' resi.obj = resi(mod, data = RESI::insurance, store.boot = TRUE, alpha = 0.01,
#' vcov.args = list(type = "HC0"), nboot = 25)
#' summary(resi.obj, alpha = c(0.1))
#' car::Anova(resi.obj, alpha = 0.1)
#'
#' # the result of resi, as well as the summary or Anova of a `resi` object can be plotted
#' # if the resi object was created with the store.boot = `TRUE` option, any alpha
#' # can be specified
#' plot(resi.obj, alpha = 0.05)
#' # if the variable names on the y-axis are too long, you can reduce their size with
#' # the ycex.axis argument (or use regular common solutions like changing the margins)
#' plot(resi.obj, alpha = 0.05, ycex.axis = 0.5)
#'
#' ## RESI on a survival model with alternate Z2S
#' # load survival library
#' library(survival)
#'
#' # fit coxph model on example data from survival package
#' # Note: for survival models, you need to specify robust variance in the model
#' # creation. resi will ignore the vcovfunc argument for this reason.
#' mod.coxph =  coxph(Surv(time, status) ~ age + sex + wt.loss, data=lung, robust = TRUE)
#'
#' # run resi on the model
#' # to use the alternative Z to RESI formula (which is equal in absolute value to the
#' # chi-square to RESI (S) formula), specify unbiased = FALSE.
#' resi(mod.coxph, data = lung, unbiased = FALSE, nboot = 10)
#'
#' @references Vandekar S, Tao R, Blume J. A Robust Effect Size Index. \emph{Psychometrika}. 2020 Mar;85(1):232-246. doi: 10.1007/s11336-020-09698-2.
#'
#' Kang, K., Armstrong, K., Avery, S., McHugo, M., Heckers, S., & Vandekar, S. (2021). Accurate confidence interval estimation for non-centrality parameters and effect size indices. \emph{arXiv preprint arXiv:2111.05966}.

resi <- function(model.full, ...){
  UseMethod("resi")
}

# take out bayes option for glm

#' @describeIn resi RESI point and interval estimation for models
#' @export
resi.default <- function(model.full, model.reduced = NULL, data, anova = TRUE,
                         coefficients = TRUE, nboot = 1000,
                         vcovfunc = sandwich::vcovHC, alpha = 0.05, store.boot = FALSE,
                         Anova.args = list(), vcov.args = list(), unbiased = TRUE, ...){
  dots = list(...)
  if ("boot.method" %in% names(dots)){
    message("Only nonparametric bootstrap supported for model type")
  }

  if (missing(data)){
    data = model.full$model
    tryCatch(update(model.full, data = data), error = function(e){
      message("Updating model fit failed. Try rerunning with providing data argument")})
  }
  else{
    data = as.data.frame(data)
  }

  # point estimation
  output <- list(alpha = alpha, nboot = nboot, boot.method = "nonparam")
  output = c(output, resi_pe(model.full = model.full, model.reduced = model.reduced,
                             data = data, anova = anova, coefficients = coefficients,
                             vcovfunc = vcovfunc, Anova.args = Anova.args,
                             vcov.args = vcov.args, unbiased = unbiased, ...))

  # bootstrapping
  boot.results = data.frame(matrix(nrow = nboot, ncol = length(output$estimates)))
  colnames(boot.results) = names(output$estimates)
  # non-parametric bootstrap
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
  }


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

#' @describeIn resi RESI point and interval estimation for models
#' @export
resi.glm <- resi.default

#' @describeIn resi RESI point and interval estimation for lm models
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
  }
  else{
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

#' @describeIn resi RESI point and interval estimation for nls models
#' @export
resi.nls <- function(model.full, model.reduced = NULL, data, coefficients = TRUE,
                     nboot = 1000, boot.method = 'nonparam',
                     vcovfunc = regtools::nlshc, alpha = 0.05, store.boot = FALSE,
                     vcov.args = list(), unbiased = TRUE, ...){

  boot.method = match.arg(tolower(boot.method), choices = c("nonparam", "bayes"))

  if (missing(data)){
    stop('\nData argument is required for nls model')
  }
  data = as.data.frame(data)

  # point estimation
  output <- list(alpha = alpha, nboot = nboot, boot.method = tolower(boot.method))
  output = c(output, resi_pe(model.full = model.full, model.reduced = model.reduced,
                             data = data, anova = anova, coefficients = coefficients,
                             vcovfunc = vcovfunc, vcov.args = vcov.args, unbiased = unbiased, ...))

  # bootstrapping
  boot.results = data.frame(matrix(nrow = nboot, ncol = length(output$estimates)))
  colnames(boot.results) = names(output$estimates)
  fail = 0
  # non-parametric bootstrap
  if (tolower(boot.method) == "nonparam"){
    for (i in 1:nboot){

      boot.data = boot.samp(data)
      t <- try(update(model.full, data = boot.data), silent = T)
      if (inherits(t, "try-error")){
        fail = fail + 1
        boot.results[i,] = NA
      }
      else{
        boot.model.full <- update(model.full, data = boot.data)

        # for catching overall error
        t1 <- try(resi_pe(model.full = boot.model.full, model.reduced = NULL,
                          data = boot.data, coefficients = coefficients,
                          vcovfunc = vcovfunc, vcov.args = vcov.args,
                          unbiased = unbiased, ...)$estimates, silent = T)
        if (inherits(t1, "try-error")){
          fail = fail + 1
        }
        else{
          boot.results[i,] = suppressWarnings(resi_pe(model.full = boot.model.full, model.reduced = NULL,
                                                      data = boot.data, coefficients = coefficients,
                                                      vcovfunc = vcovfunc, vcov.args = vcov.args,
                                                      unbiased = unbiased, ...)$estimates)
        }

      }

    }}

  # bayesian bootstrap
  if (tolower(boot.method)  == "bayes"){
    g = NULL
    for (i in 1:nboot){
      boot.data = bayes.samp(data)

      t <- try(update(model.full, data = boot.data, weights = g), silent = T)
      if (inherits(t, "try-error")){
        fail = fail + 1
        boot.results[i,] = NA
      }
      else{
        boot.model.full <- update(model.full, data = boot.data, weights = g)

        # for catching overall error
        t1 <- try(resi_pe(model.full = boot.model.full, model.reduced = NULL,
                          data = boot.data, coefficients = coefficients,
                          vcovfunc = vcovfunc, vcov.args = vcov.args,
                          unbiased = unbiased, ...)$estimates, silent = T)
        if (inherits(t1, "try-error")){
          fail = fail + 1
        }
        else{
          boot.results[i,] = suppressWarnings(resi_pe(model.full = boot.model.full, model.reduced = NULL,
                                                      data = boot.data, coefficients = coefficients,
                                                      vcovfunc = vcovfunc, vcov.args = vcov.args,
                                                      unbiased = unbiased, ...)$estimates)
        }}}}

  alpha.order = sort(c(alpha/2, 1-alpha/2))
  output$overall[nrow(output$overall),paste(alpha.order*100, '%', sep='')] = quantile(boot.results[,1], probs = alpha.order, na.rm = TRUE)

  if (coefficients){
    CIs = apply(boot.results[,2:(1+nrow(output$coefficients))], 2,  quantile, probs = alpha.order, na.rm = TRUE)
    CIs = t(CIs)
    output$coefficients[1:nrow(CIs), paste(alpha.order*100, '%', sep='')] = CIs
  }

  if(store.boot){
    output$boot.results = boot.results
  }
  output$nfail = fail
  class(output) = "resi"
  return(output)
}

#' @describeIn resi RESI point and interval estimation for survreg models
#' @export
resi.survreg <- function(model.full, model.reduced = NULL, data, anova = TRUE,
                         coefficients = TRUE, nboot = 1000,
                         vcovfunc = vcov, alpha = 0.05, store.boot = FALSE,
                         Anova.args = list(), unbiased = TRUE, ...){

  if (missing(data)){
    stop('\nData argument is required for survreg model')
  }
  else{
    data = as.data.frame(data)
  }

  dots = list(...)
  if ("vcov.args" %in% names(dots)){
    warning("vcov.args ignored for survreg objects")
  }
  if ("boot.method" %in% names(dots)){
    message("Only nonparametric bootstrap supported for model type")
  }

  resi.default(model.full = model.full, model.reduced = model.reduced, data = data,
               anova = anova, coefficients = coefficients, nboot = nboot, vcovfunc = vcovfunc,
               store.boot = store.boot, Anova.args = Anova.args,
               unbiased = unbiased, alpha = alpha, ...)
}

#' @describeIn resi RESI point and interval estimation for coxph models
#' @export
resi.coxph <- function(model.full, model.reduced = NULL, data, anova = TRUE,
                       coefficients = TRUE, nboot = 1000,
                       vcovfunc = vcov, alpha = 0.05, store.boot = FALSE,
                       Anova.args = list(), unbiased = TRUE, ...){
  if (missing(data)){
    stop('\nData argument is required for coxph model')
  }
  else{
    data = as.data.frame(data)
  }

  dots = list(...)
  if ("vcov.args" %in% names(dots)){
    warning("vcov.args ignored for survreg objects")
  }
  if ("boot.method" %in% names(dots)){
    message("Only nonparametric bootstrap supported for model type")
  }

  resi.default(model.full = model.full, model.reduced = model.reduced, data = data,
               anova = anova, coefficients = coefficients, nboot = nboot, vcovfunc = vcovfunc,
               store.boot = store.boot, Anova.args = Anova.args,
               unbiased = unbiased, alpha = alpha, ...)
}

#' @describeIn resi RESI point and interval estimation for hurdle models
#' @export
resi.hurdle <- function(model.full, model.reduced = NULL, data, coefficients = TRUE,
                     nboot = 1000, vcovfunc = sandwich::sandwich,
                     alpha = 0.05, store.boot = FALSE, vcov.args = list(), unbiased = TRUE, ...){

  dots = list(...)
  if ("boot.method" %in% names(dots)){
    message("Only nonparametric bootstrap supported for model type")
  }

  if (missing(data)){
    data = model.full$model
    tryCatch(update(model.full, data = data), error = function(e){
      message("Updating model fit failed. Try rerunning with providing data argument")})
  }
  else{
    data = as.data.frame(data)
  }

  if (is.null(model.reduced)){
    model.reduced = update(model.full, formula = as.formula(paste(format(formula(model.full)[[2]]), "~ 1")),  data = model.full$model)
  }

  # point estimation
  output <- list(alpha = alpha, nboot = nboot, boot.method = "nonparam")
  output = c(output, resi_pe(model.full = model.full, model.reduced = model.reduced,
                             data = data, anova = anova, coefficients = coefficients,
                             vcovfunc = vcovfunc, vcov.args = vcov.args, unbiased = unbiased, ...))

  # bootstrapping
  boot.results = data.frame(matrix(nrow = nboot, ncol = length(output$estimates)))
  colnames(boot.results) = names(output$estimates)
  # non-parametric bootstrap
  for (i in 1:nboot){
    boot.data = boot.samp(data)
    boot.model.full <- update(model.full, data = boot.data)
    boot.model.reduced <- update(model.full, data = boot.data)
    boot.results[i,] = suppressWarnings(resi_pe(model.full = boot.model.full, model.reduced = NULL,
                                                data = boot.data, coefficients = coefficients,
                                                vcovfunc = vcovfunc, vcov.args = vcov.args,
                                                unbiased = unbiased, ...)$estimates)

  }


  output$overall[nrow(output$overall),c(paste(alpha/2*100, '%', sep=''), paste((1-alpha/2)*100, '%', sep=''))] = quantile(boot.results[,1], probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)

  if (coefficients){
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
resi.geeglm <- function(model.full, data, anova = TRUE,
                        coefficients = TRUE, nboot = 1000,
                        alpha = 0.05, store.boot = FALSE,
                        unbiased = TRUE, ...){
  dots = list(...)
  if ("boot.method" %in% names(dots)){
    message("Only nonparametric bootstrap supported for model type")
  }

  if (missing(data)){
    data = model.full$data
    tryCatch(update(model.full, data = data), error = function(e){
      message("Updating model fit failed. Try rerunning with providing data argument")})
  }
  else{
    data = as.data.frame(data)
  }

  # point estimation
  output <- list(alpha = alpha, nboot = nboot, boot.method = "nonparam")
  output = c(output, resi_pe(object = model.full, data = data, anova = anova,
                             coefficients = coefficients, unbiased = unbiased, ...))

  # id variable name
  id_var = as.character(model.full$call$id)
  corstr_spec = model.full$corstr
  # bootstrapping
  boot.results = data.frame(matrix(nrow = nboot, ncol = length(output$estimates)))
  colnames(boot.results) = names(output$estimates)
  fail = 0
  for (i in 1:nboot){
    skip_to_next <- FALSE
    boot.data = boot.samp(data, id.var = id_var)
    # re-fit the model
    boot.model.full = update(model.full, data = boot.data, corstr = corstr_spec,
                             id = bootid)
    rv.boot = tryCatch(resi_pe(boot.model.full, data = boot.data, anova = anova,
                               coefficients = coefficients, unbiased = unbiased, ...),
                       error = function(e){skip_to_next <<- TRUE})
    if (skip_to_next) {
      fail = fail + 1
      next
      }
    # output.boot$RESI= cbind(output.boot$RESI, rv.boot$resi[, 'RESI'])
    boot.results[i,] = resi_pe(object = boot.model.full,
                              data = boot.data, anova = anova,
                              coefficients = coefficients,
                              unbiased = unbiased)$resi
  }

  alpha.order = sort(c(alpha/2, 1-alpha/2))

  if (coefficients){
    CIs = apply(boot.results[,1:nrow(output$coefficients)], 2,  quantile,
                probs = alpha.order, na.rm = TRUE)
    CIs = t(CIs)
    output$coefficients[1:nrow(CIs), paste(alpha.order*100, '%', sep='')] = CIs
  }

  if (anova){
    CIs = apply(boot.results[,(ncol(boot.results)-
                                 length(rownames(output$anova))+1):ncol(boot.results)],
                2,  quantile, probs = alpha.order, na.rm = TRUE)
    CIs = t(CIs)
    output$anova[1:nrow(CIs), paste(alpha.order*100, '%', sep='')] = CIs
    class(output$anova) = c("anova_resi", class(output$anova))
  }

  if(store.boot){
    output$boot.results = boot.results
  }

  output$nfail = fail

  class(output) = "resi"
  return(output)

  # RESI_se = apply(output.boot$RESI, 1, sd, na.rm = TRUE)

  return(output)
}


#' @describeIn resi RESI point and interval estimation for GEE models
#' @export
resi.gee <- function(model.full, data, nboot = 1000, alpha = 0.05,
                     store.boot = FALSE, unbiased = TRUE, ...){
  if (missing(data)){
    stop('\nData argument is required for GEE models from gee package')
  }
  else{
    data = as.data.frame(data)
  }

  output <- list(alpha = alpha, nboot = nboot, boot.method = "nonparam")
  # RESI point estimates
  output = c(output, resi_pe(model.full, unbiased, ...))
  # id variable name
  id_var = as.character(model.full$call$id)

  # bootstrapping
  boot.results = data.frame(matrix(nrow = nboot, ncol = length(output$estimates)))
  colnames(boot.results) = names(output$estimates)
  for (i in 1:nboot){
    boot.data = boot.samp(data, id.var = id_var)
    # re-fit the model
    suppressMessages(capture.output(boot.model.full <-
                                      update(model.full, data = boot.data, id = bootid),
                                    file =  nullfile()))
    boot.results[i,] = suppressWarnings(resi_pe(model.full = boot.model.full,
                                                unbiased = unbiased, ...)$estimates)
  }

  alpha.order = sort(c(alpha/2, 1-alpha/2))

  CIs = apply(boot.results[,1:nrow(output$coefficients)], 2,  quantile,
                probs = alpha.order, na.rm = TRUE)
  CIs = t(CIs)
  output$coefficients[1:nrow(CIs), paste(alpha.order*100, '%', sep='')] = CIs

  if(store.boot){
    output$boot.results = boot.results
  }

  class(output)= 'resi'
  return(output)
}

# convergence issues - no fix yet

#' @describeIn resi RESI point and interval estimation for LME (nlme) models
#' @importFrom nlme getGroups
#' @export
resi.lme <- function(model.full, alpha = 0.05, nboot = 1000, vcovfunc = clubSandwich::vcovCR,
                     vcov.args = list(), ...){
  warning("\nConfidence Interval procedure not developed for lme, returning point estimates only")
  output <- list(alpha = alpha, nboot = 0)
  # RESI point estimates
  output = c(output, resi_pe(model.full = model.full, vcovfunc = vcovfunc, vcov.args = vcov.args))
  data = model.full$data
  # id variable name
  id_var = attr(nlme::getGroups(model.full), "label")
  # bootstrap
  #output.boot = as.matrix(output$coefficients[, 'RESI'])
  #fun <- utils::getFromNamespace("update.lme", "nlme")
  #tryCatch(fun(model.full, data = data), error = function(e){
    #message("Try running `library(nlme)`")})
  # for (i in 1:nboot){
  #   boot.data = boot.samp(data, id.var = id_var)
  #   # re-fit the model
  #   boot.mod = fun(model.full, data = boot.data, fixed = as.formula(model.full$call$fixed),
  #                     random = as.formula(model.full$call$random))
  #   output.boot = cbind(output.boot, resi_pe(model.full = boot.mod, vcovfunc = vcovfunc,
  #                                            vcov.args = vcov.args)$coefficients[, 'RESI'])
  # }
  # output.boot = output.boot[, -1]
  # RESI.ci = apply(output.boot, 1, quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
  # output$coefficients = cbind(output$coefficients, t(RESI.ci))
  output$boot.method = "nonparam"
  class(output) = 'resi'
  return(output)
}

#' @describeIn resi RESI point and interval estimation for lmerMod models
#' @export
resi.lmerMod<- function(model.full, alpha = 0.05, nboot = 1000,
                        vcovfunc = clubSandwich::vcovCR, vcov.args = list(), ...){
  # warning("\nConfidence Interval procedure not developed for lmerMod, returning point estimates only")
  output <- list(alpha = alpha, nboot = nboot)
  output = c(output, resi_pe(model.full, vcovfunc = vcovfunc, vcov.args = vcov.args)) # RESI point estimates
  data = model.full@frame
  # id variable name
  id_var = names(model.full@flist)
  form_full = formula(model.full) %>% format
  form_full_boot = gsub(id_var, "bootid", form_full) %>% as.formula
  # bootstrap (non-param for now)
  output.boot = as.matrix(output$coefficients[, 'RESI'])
  for (i in 1:nboot){
    boot.data = boot.samp(data, id.var = id_var)
    # re-fit the model
    boot.mod = lme4::lmer(form_full_boot, data = boot.data) # using `bootid` as the new ID
    # update(model.full, formula = form_full_boot, data = boot.data)
    rv.boot = resi_pe(boot.mod, vcovfunc = vcovfunc, vcov.args = vcov.args)
    output.boot = cbind(output.boot, rv.boot$coefficients[, 'RESI'])
  }
  output.boot = output.boot[, -1]
  RESI.ci = apply(output.boot, 1, quantile, probs = c(alpha/2, 1-alpha/2))
  output$coefficients = cbind(output$coefficients, t(RESI.ci))
  output$boot.method = "nonparam"
  class(output) = 'resi'
  return(output)
}
