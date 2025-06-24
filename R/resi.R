#' Robust Effect Size Index (RESI) point and interval estimation for models
#'
#' This function will estimate the robust effect size (RESI) from Vandekar, Tao, & Blume (2020) and its confidence interval in various ways for a fitted model object. The overall RESI is estimated via a Wald test. RESI is (optionally) estimated for each factor in coefficients-style table. RESI is (optionally) estimated for each variable/interaction in an Anova-style table for models with existing Anova methods. CIs can be calculated using either non-parametric or Bayesian bootstrapping.
#' @param model.full \code{lm, glm, nls, survreg, coxph, hurdle, zeroinfl, gee, geeglm} or \code{lme} model object.
#' @param model.reduced Fitted model object of same type as model.full. By default `NULL`; the same model as the full model but only having intercept.
#' @param data Data.frame or object coercible to data.frame of model.full data (required for some model types).
#' @param vcovfunc The variance estimator function for constructing the Wald test statistic. By default, \link[sandwich]{vcovHC} (the robust (sandwich) variance estimator).
#' @param coefficients Logical, whether to produce a coefficients (summary) table with the RESI columns added. By default = `TRUE`.
#' @param anova Logical, whether to produce an Anova table with the RESI columns added. By default = `TRUE`.
#' @param overall Logical, whether to produce an overall Wald test comparing full to reduced model with RESI columns added. By default = `TRUE`.
#' @param nboot Numeric, the number of bootstrap replicates. By default, 1000.
#' @param boot.method String, which type of bootstrap to use: `nonparam` = non-parametric bootstrap (default); `bayes` = Bayesian bootstrap.
#' @param alpha Numeric, significance level of the constructed CIs. By default, 0.05.
#' @param store.boot Logical, whether to store all the bootstrapped estimates. By default, `FALSE`.
#' @param Anova.args List, additional arguments to be passed to \link[car]{Anova} function.
#' @param vcov.args List, additional arguments to be passed to vcovfunc.
#' @param unbiased Logical, whether to use the unbiased or alternative T/Z statistic to RESI conversion. By default, `TRUE`. See details.
#' @param parallel See documentation for \link[boot]{boot}.
#' @param ncpus See documentation for \link[boot]{boot}.
#' @param long Logical, whether the data is longitudinal/clustered. By default, `FALSE`.
#' @param clvar Character, the name of the cluster/id variable if data is clustered. By default, `NULL`.
#' @param ... Ignored.
#' @importFrom aod wald.test
#' @importFrom boot boot
#' @importFrom car Anova
#' @importFrom lmtest waldtest
#' @importFrom sandwich vcovHC
#' @importFrom stats anova as.formula coef formula glm hatvalues lm nobs pchisq pf predict quantile rbinom residuals rnorm runif update vcov
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
#' \code{\link{t2S}}, and \code{\link{z2S}} for more details on the formulas.
#'
#' For GEE (\link[geepack]{geeglm}, \link[glmtoolbox]{glmgee}) models, a longitudinal RESI (L-RESI) and a cross-sectional,
#' per-measurement RESI (CS-RESI) is estimated. The longitudinal RESI takes the
#' specified clustering into account, while the cross-sectional RESI is designed to
#' estimate the effect size if a random observation for each participant were collected cross-sectionally.
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
#' @seealso \code{\link{resi_pe}}, \link[sandwich]{vcovHC},
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
#' # store bootstrap results for calculating different CIs later
#' resi_obj = resi(mod, nboot = 50, store.boot = TRUE)
#' # print output
#' resi_obj
#'
#' # fit a reduced model for comparison
#' mod_red = lm(charges ~ bmi, data = RESI::insurance)
#'
#' # running resi and including the reduced model will provide almost the exact same
#' # output as not including a reduced model. The difference is that the "overall"
#' # portion of the output will compare the full model to the reduced model.
#' # The "summary" and "anova" RESI estimates will be the same. (The bootstrapped
#' # confidence intervals may differ.)
#' resi(model.full = mod, model.reduced = mod_red, nboot = 10)
#'
#' # used stored bootstrap results to get a different alpha-level confidence interval
#' summary(resi_obj, alpha = c(0.01, 0.1))
#' car::Anova(resi_obj, alpha = c(0.01, 0.1))
#'
#' # the result of resi, as well as the summary or Anova of a `resi` object can be plotted
#' # if the resi object was created with the store.boot = `TRUE` option, any alpha
#' # can be specified
#' plot(resi_obj, alpha = 0.01)
#' # if the variable names on the y-axis are too long, you can reduce their size with
#' # the ycex.axis argument (or use regular common solutions like changing the margins)
#' plot(resi_obj, alpha = 0.01, ycex.axis = 0.5)
#'
#' # for some model types and formula structures, data argument is required
#' if(requireNamespace("splines")){
#'   # fit logistic regression model with splines
#'   mod = glm(smoker ~ splines::ns(age, df = 3) + region, data = RESI::insurance,
#'     family = "binomial")
#'
#'   # specify additional arguments to the variance-covariance function via vcov.args
#'   resi_obj = resi(mod, data = RESI::insurance, alpha = 0.01,
#'     vcov.args = list(type = "HC0"), nboot = 25)
#'   summary(resi_obj)
#'   car::Anova(resi_obj)}
#'
#'
#' ## RESI on a survival model with alternate Z2S
#' if(requireNamespace("survival")){
#'   # fit coxph model on example data from survival package
#'   # Note: for survival models, you need to specify robust variance in the model
#'   # creation. resi will ignore the vcovfunc argument for this reason.
#'   mod.coxph =  survival::coxph(survival::Surv(time, status) ~ age + sex + wt.loss,
#'    data=survival::lung, robust = TRUE)
#'
#'   # run resi on the model
#'   # to use the alternative Z to RESI formula (which is equal in absolute value to the
#'   # chi-square to RESI (S) formula), specify unbiased = FALSE.
#'   resi(mod.coxph, data = survival::lung, unbiased = FALSE, nboot = 10)}
#'
#' @references Vandekar S, Tao R, Blume J. A Robust Effect Size Index. \emph{Psychometrika}. 2020 Mar;85(1):232-246. doi: 10.1007/s11336-020-09698-2.
#'
#' Kang, K., Armstrong, K., Avery, S., McHugo, M., Heckers, S., & Vandekar, S. (2021). Accurate confidence interval estimation for non-centrality parameters and effect size indices. \emph{arXiv preprint arXiv:2111.05966}.
#'
#' Jones, M., Kang, K., Vandekar, S. (2025). \emph{Journal of Statistical Software}. RESI: An R Package for Robust Effect Sizes.<doi:10.18637/jss.v112.i03>

resi = function(model.full, ...){
  UseMethod("resi")
}

#' @describeIn resi RESI point and interval estimation for models
#' @export
resi.default = function(model.full, model.reduced = NULL, data, anova = TRUE,
                         coefficients = TRUE, overall = TRUE, nboot = 1000,
                         boot.method = "nonparam", vcovfunc = sandwich::vcovHC,
                         alpha = 0.05, store.boot = FALSE, Anova.args = list(),
                         vcov.args = list(), unbiased = TRUE,
                         parallel = c("no", "multicore", "snow"),
                         ncpus = getOption("boot.ncpus", 1L), long = FALSE,
                         clvar = NULL, ...){

  # check for supported model type
  if(!(any(class(model.full) %in%
    sapply(as.character(utils::methods(RESI::resi)), function(x) substr(x, 6, nchar(x)))))){
    warning("Model type not implemented in RESI package, attempting default")
  }
  dots = list(...)

  if (missing(data)){
    data = model.full$model
    tryCatch(update(model.full, data = data), error = function(e){
      message("Updating model fit failed. Try running with providing data argument")})
  }
  else{
    if (!(is.null(model.full$na.action))){
      data = data[which(!(1:nrow(data)%in% model.full$na.action)),]
    }
    data = as.data.frame(data)
  }

  if (is.null(model.reduced) & overall){
    if(!("skip.red" %in% names(dots))){
      form.reduced = as.formula(paste(format(formula(model.full)[[2]]), "~ 1"))
      if (!(form.reduced == formula(model.full))){
        if(!(is.null(model.full$model)) & !long){
          model.reduced = try(update(model.full, formula = form.reduced,
                                      data = model.full$model), silent = T)
        } else{
          model.reduced = try(update(model.full, formula = form.reduced,
                                      data = data), silent = T)
        }
        if(inherits(model.reduced, "try-error")){
          warning("Fitting intercept-only model failed. No overall test computed")
          overall = FALSE
        }
    }}}

  # point estimation
  output = list(alpha = alpha, nboot = nboot, boot.method = boot.method)
  output = c(output, resi_pe(model.full = model.full, model.reduced = model.reduced,
                             data = data, anova = anova, coefficients = coefficients,
                             overall = overall, vcovfunc = vcovfunc,
                             Anova.args = Anova.args, vcov.args = vcov.args,
                             unbiased = unbiased, ...))

  if (long){
    data = unique(data[, clvar])
  }

  # bootstrapping
  suppressMessages(capture.output(boot_out <- boot::boot(data = data, R = nboot, statistic = resi_stat, parallel = parallel,
                  ncpus = ncpus, mod.full = model.full, mod.reduced = model.reduced,
                  anova = anova, coefficients = coefficients, overall = overall, vcovfunc = vcovfunc,
                  Anova.args = Anova.args, vcov.args = vcov.args, unbiased = unbiased,
                  boot.method = boot.method, clvar = clvar, nest = length(output$estimates), ...), file = nullfile()))

  # bootstrapped estimates
  boot.results = boot_out$t
  colnames(boot.results) = names(boot_out$t0)

  if (all(is.na(boot.results))){
    stop("\nBootstrapping failed. Run resi_pe for point estimates")
  }

  alpha.order = sort(c(alpha/2, 1-alpha/2))

  if (!long){
    if (overall & !is.null(output$overall)){
      if(nrow(output$overall) == 1){
        names.overall = colnames(output$overall)
        for (i in 1:length(alpha.order)){
          output$overall = cbind(output$overall, quantile(boot.results[,1],
                                                          probs = alpha.order[i], na.rm = TRUE))
        }
        colnames(output$overall) = c(names.overall, paste(alpha.order*100, "%", sep=""))
      } else{
        output$overall[nrow(output$overall),paste(alpha.order*100, "%", sep="")] =
          quantile(boot.results[,1], probs = alpha.order, na.rm = TRUE)}
    }

    if (coefficients){
      co = boot.results[,(1+!(is.null(output$overall))):((!is.null(output$overall))+nrow(output$coefficients))]
      if (is.null(dim(co))){
        CIs = quantile(co, probs = alpha.order, na.rm = TRUE)
      } else{
        CIs = apply(co, 2,  quantile, probs = alpha.order, na.rm = TRUE)
      }
      CIs = t(CIs)
      output$coefficients[1:nrow(CIs), paste(alpha.order*100, "%", sep="")] = CIs
    }

    if (anova){
      an = boot.results[,(ncol(boot.results)-length(which(rownames(output$anova) != "Residuals"))+1):ncol(boot.results)]
      if (is.null(dim(an))){
        CIs = quantile(an, probs = alpha.order, na.rm = TRUE)
      } else{
        CIs = apply(an, 2,  quantile, probs = alpha.order, na.rm = TRUE)
      }
      CIs = t(CIs)
      output$anova[1:nrow(CIs), paste(alpha.order*100, "%", sep="")] = CIs
    }
  } else {
    r = 1
    if (!is.null(output$overall)){
      # L-RESI for geeglm model overall Wald test
      output$overall[nrow(output$overall),paste("L ",alpha.order*100, "%", sep="")] =
        quantile(boot.results[,1], probs = alpha.order, na.rm = TRUE)
      r = 2
    }
    # L-RESI and CS-RESI for longitudinal models
    if (coefficients){
      lco = boot.results[,r:(nrow(output$coefficients) + r -1)]
      if (is.null(dim(lco))){
        lCIs = quantile(lco, probs = alpha.order, na.rm = TRUE)
      } else{
        lCIs = apply(lco, 2,  quantile,
                     probs = alpha.order, na.rm = TRUE)
      }
      lCIs = t(lCIs)
      output$coefficients[1:nrow(lCIs), paste("L ",alpha.order*100, "%", sep="")] = lCIs
      cco = boot.results[,(nrow(output$coefficients)+r):(2*nrow(output$coefficients) + r - 1)]
      if (is.null(dim(cco))){
        cCIs = quantile(cco, probs = alpha.order, na.rm = TRUE)
      } else{
        cCIs = apply(cco, 2,  quantile, probs = alpha.order, na.rm = TRUE)
      }
      cCIs = t(cCIs)
      output$coefficients[1:nrow(cCIs), paste("CS ",alpha.order*100, "%", sep="")] = cCIs
    }

    if (anova){
      lan = boot.results[,(ncol(boot.results)-
                             2*length(rownames(output$anova))+1):
                           (ncol(boot.results)-length(rownames(output$anova)))]
      if (is.null(dim(lan))){
        lCIs = quantile(lan, probs = alpha.order, na.rm = TRUE)
      } else{
        lCIs = apply(lan, 2,  quantile, probs = alpha.order, na.rm = TRUE)
      }
      lCIs = t(lCIs)
      output$anova[1:nrow(lCIs), paste("L ", alpha.order*100, "%", sep="")] = lCIs
      can = boot.results[,(ncol(boot.results)-length(rownames(output$anova))+1):
                           ncol(boot.results)]
      if (is.null(dim(can))){
        cCIs = quantile(can, probs = alpha.order, na.rm = TRUE)
      } else{
        cCIs = apply(can, 2,  quantile, probs = alpha.order, na.rm = TRUE)
      }
      cCIs = t(cCIs)
      output$anova[1:nrow(cCIs), paste("CS ", alpha.order*100, "%", sep="")] = cCIs
    }
  }

  if(store.boot){
    #output$boot.results = boot.results
    output$boot.results = boot_out
  }
  class(output) = "resi"
  return(output)
}

#' @describeIn resi RESI point and interval estimation for models
#' @export
resi.glm = function(model.full, model.reduced = NULL, data, anova = TRUE,
                    coefficients = TRUE, overall = TRUE, nboot = 1000,
                    vcovfunc = sandwich::vcovHC, alpha = 0.05, store.boot = FALSE,
                    Anova.args = list(), vcov.args = list(), unbiased = TRUE,
                    parallel = c("no", "multicore", "snow"), ncpus = getOption("boot.ncpus", 1L),
                    ...){
  dots = list(...)
  if ("boot.method" %in% names(dots)){
    stop("\nOnly nonparametric bootstrap supported for model type")
  }
  resi.default(model.full = model.full, model.reduced = model.reduced, data = data,
               anova = anova, coefficients = coefficients, overall = overall,
               nboot = nboot, vcovfunc = vcovfunc, store.boot = store.boot,
               Anova.args = Anova.args, vcov.args = vcov.args,
               unbiased = unbiased, alpha = alpha, parallel = parallel,
               ncpus = ncpus, ...)

}

#' @describeIn resi RESI point and interval estimation for lm models
#' @export
resi.lm = function(model.full, model.reduced = NULL, data, anova = TRUE,
                    coefficients = TRUE, overall = TRUE, nboot = 1000,
                    boot.method = "nonparam", vcovfunc = sandwich::vcovHC,
                    alpha = 0.05, store.boot = FALSE, Anova.args = list(),
                    vcov.args = list(), unbiased = TRUE,
                    parallel = c("no", "multicore", "snow"),
                    ncpus = getOption("boot.ncpus", 1L), ...){
  boot.method = match.arg(tolower(boot.method), choices = c("nonparam", "bayes"))

  resi.default(model.full = model.full, model.reduced = model.reduced, data = data,
               anova = anova, coefficients = coefficients, overall = overall,
               nboot = nboot, vcovfunc = vcovfunc, store.boot = store.boot,
               Anova.args = Anova.args, vcov.args = vcov.args, unbiased = unbiased,
               alpha = alpha, parallel = parallel, ncpus = ncpus,
               boot.method = boot.method, ...)

}

#' @describeIn resi RESI point and interval estimation for nls models
#' @export
resi.nls = function(model.full, model.reduced = NULL, data, coefficients = TRUE,
                     overall = TRUE, nboot = 1000, boot.method = "nonparam",
                     anova = FALSE, vcovfunc = r_nlshc, alpha = 0.05,
                     store.boot = FALSE, vcov.args = list(), unbiased = TRUE,
                     parallel = c("no", "multicore", "snow"), ncpus = getOption("boot.ncpus", 1L),
                     ...){

  boot.method = match.arg(tolower(boot.method), choices = c("nonparam", "bayes"))

  if (missing(data)){
    stop("\nData argument is required for nls model")
  }
  if (identical(vcovfunc, sandwich::vcovHC)){
    vcovfunc = r_nlshc
    warning("Sandwich vcov function not applicable for nls model type, vcovfunc set to regtools::nlshc")
  }

  output = resi.default(model.full = model.full, model.reduced = model.reduced,
                        data = data, anova = FALSE, coefficients = coefficients,
                        overall = overall, nboot = nboot, vcovfunc = vcovfunc,
                        store.boot = TRUE, unbiased = unbiased, alpha = alpha,
                        vcov.args = vcov.args, parallel = parallel,
                        ncpus = ncpus, boot.method = boot.method, skip.red = TRUE, ...)


  # number of bootstrap replicates with failed updating
  output$nfail = length(which(is.na(output$boot.results$t[,1])))

  if (!store.boot){
    output = output[names(output) != "boot.results"]
  }
  class(output) = "resi"
  return(output)
}

#' @describeIn resi RESI point and interval estimation for survreg models
#' @export
resi.survreg = function(model.full, model.reduced = NULL, data, anova = TRUE,
                         coefficients = TRUE, overall = TRUE, nboot = 1000,
                         vcovfunc = vcov, alpha = 0.05, store.boot = FALSE,
                         Anova.args = list(), unbiased = TRUE,
                         parallel = c("no", "multicore", "snow"),
                         ncpus = getOption("boot.ncpus", 1L),...){
  dots = list(...)
  if (missing(data)){
    stop("\nData argument is required for model type")
  }
  if (!identical(vcovfunc, vcov)){
    warning("vcovfunc argument ignored for survreg objects")
  }
  if ("vcov.args" %in% names(dots)){
    warning("vcov.args ignored for model type")
  }
  if ("boot.method" %in% names(dots)){
    stop("\nOnly nonparametric bootstrap supported for model type")
  }

  resi.default(model.full = model.full, model.reduced = model.reduced, data = data,
               anova = anova, coefficients = coefficients, overall = overall,
               nboot = nboot, vcovfunc = vcov, store.boot = store.boot,
               Anova.args = Anova.args, unbiased = unbiased, alpha = alpha,
               boot.method = "nonparam", parallel = parallel, ncpus = ncpus, ...)
}

#' @describeIn resi RESI point and interval estimation for coxph models
#' @export
resi.coxph = function(model.full, model.reduced = NULL, data, anova = TRUE,
                       coefficients = TRUE, overall = TRUE, nboot = 1000,
                       vcovfunc = vcov, alpha = 0.05, store.boot = FALSE,
                       Anova.args = list(), unbiased = TRUE,
                       parallel = c("no", "multicore", "snow"),
                       ncpus = getOption("boot.ncpus", 1L),...){
  dots = list(...)
  if (missing(data)){
    stop("\nData argument is required for model type")
  }
  if (!identical(vcovfunc, vcov)){
    warning("vcovfunc argument ignored for coxph objects")
  }
  if ("vcov.args" %in% names(dots)){
    warning("vcov.args ignored for model type")
  }
  if ("boot.method" %in% names(dots)){
    stop("\nOnly nonparametric bootstrap supported for model type")
  }
  if (!is.null(model.reduced)){
    warning("Reduced model argument ignored for coxph model")
  }
  resi.default(model.full = model.full, model.reduced = NULL, data = data,
               anova = anova, coefficients = coefficients, overall = overall,
               nboot = nboot, vcovfunc = vcov, store.boot = store.boot,
               Anova.args = Anova.args, unbiased = unbiased, alpha = alpha,
               boot.method = "nonparam", parallel = parallel, ncpus = ncpus,
               skip.red = TRUE, ...)
}


#' @describeIn resi RESI point and interval estimation for hurdle models
#' @export
resi.hurdle = function(model.full, model.reduced = NULL, data, coefficients = TRUE,
                        overall = TRUE, nboot = 1000, vcovfunc = sandwich::sandwich,
                        anova = FALSE, alpha = 0.05, store.boot = FALSE,
                        vcov.args = list(), unbiased = TRUE,
                        parallel = c("no", "multicore", "snow"),
                        ncpus = getOption("boot.ncpus", 1L), ...){

  dots = list(...)
  if ("boot.method" %in% names(dots)){
    stop("\nOnly nonparametric bootstrap supported for model type")
  }

  resi.default(model.full = model.full, model.reduced = model.reduced, data = data,
               coefficients = coefficients, overall = overall, nboot = nboot,
               anova = FALSE, vcovfunc = vcovfunc, store.boot = store.boot,
               vcov.args = vcov.args, unbiased = unbiased, alpha = alpha,
               boot.method = "nonparam", parallel = parallel, ncpus = ncpus, ...)

}

#' @describeIn resi RESI point and interval estimation for zeroinfl models
#' @export
resi.zeroinfl = resi.hurdle

#' @describeIn resi RESI point and interval estimation for GEE models
#' @export
resi.geeglm = function(model.full, model.reduced = NULL, data, anova = TRUE,
                        coefficients = TRUE, overall = TRUE, nboot = 1000,
                        alpha = 0.05, store.boot = FALSE,
                        unbiased = TRUE, parallel = c("no", "multicore", "snow"),
                        ncpus = getOption("boot.ncpus", 1L), ...){
  dots = list(...)
  if ("boot.method" %in% names(dots)){
    stop("\nOnly nonparametric bootstrap supported for model type")
  }

  if (missing(data)){
    data = model.full$data
    tryCatch(update(model.full, data = data), error = function(e){
      message("Updating model fit failed. Try rerunning with providing data argument")})
  }
  else{
    data = as.data.frame(data)
  }

  # id variable name
  id_var = as.character(model.full$call$id)

  # bootstrapping
  output = resi.default(model.full = model.full, model.reduced = model.reduced,
                        data = data, anova = anova, coefficients = coefficients,
                        overall = overall, nboot = nboot,
                        store.boot = TRUE,
                        unbiased = unbiased, alpha = alpha,
                        parallel = parallel, ncpus = ncpus, boot.method = "nonparam",
                        cluster = TRUE, clvar = id_var, mod.dat = data, long = TRUE, ...)

  # number of bootstrap replicates with failed updating
  output$nfail = length(which(is.na(output$boot.results$t[,1])))


  if (!store.boot){
    output = output[names(output) != "boot.results"]
  }

  class(output) = "resi"

  return(output)
}

#' @describeIn resi RESI point and interval estimation for GEE models in glmtoolbox
#' @export
resi.glmgee = function(model.full, model.reduced = NULL, data, anova = FALSE,
                                     coefficients = TRUE, overall = TRUE, nboot = 1000,
                                     alpha = 0.05, store.boot = FALSE,
                                     unbiased = TRUE, parallel = c("no", "multicore", "snow"),
                                     ncpus = getOption("boot.ncpus", 1L), ...){
  # Capture the call to my_function
  call <- match.call()

  # Get the argument names and their default values
  formal_args <- formals(resi.glmgee)

  # Create a list of all arguments including defaults
  args_list <- as.list(formal_args)
  # remove ellipses
  args_list = args_list[-length(args_list)]

  # Update the list with the user's inputs from the call (if provided)
  for (name in names(call)[-1]) {  # Exclude the function name
    args_list[[name]] <- call[[name]]
  }
  do.call(resi.geeglm, args_list)
}

 #' @describeIn resi RESI point and interval estimation for GEE models
#' @export
resi.gee = function(model.full, data, nboot = 1000, alpha = 0.05,
                     store.boot = FALSE, unbiased = TRUE,
                     parallel = c("no", "multicore", "snow"),
                     ncpus = getOption("boot.ncpus", 1L), ...){
  dots = list(...)
  if ("boot.method" %in% names(dots)){
    stop("\nOnly nonparametric bootstrap supported for model type")
  }
  if (missing(data)){
    stop("\nData argument is required for GEE models from gee package")
  }
  else{
    data = as.data.frame(data)
  }

  id_var = as.character(model.full$call$id)

  # bootstrapping
  output = resi.default(model.full = model.full, model.reduced = NULL, data = data,
                        coefficients = TRUE, nboot = nboot,
                        store.boot = TRUE, anova = FALSE,
                        unbiased = unbiased, alpha = alpha, parallel = parallel,
                        ncpus = ncpus, boot.method = "nonparam", cluster = TRUE,
                        clvar = id_var, mod.dat = data, skip.red = TRUE, long = TRUE,
                        overall = FALSE, ...)

  # number of bootstrap replicates with failed updating
  output$nfail = length(which(is.na(output$boot.results$t[,1])))

  if (!store.boot){
    output = output[names(output) != "boot.results"]
  }

  class(output) = "resi"

  return(output)
}

#' @describeIn resi RESI point and interval estimation for LME (nlme) models
#' @importFrom nlme getGroups
#' @export
resi.lme = function(model.full, alpha = 0.05, nboot = 1000, anova = TRUE, vcovfunc = clubSandwich::vcovCR,
                     vcov.args = list(), ...){
  warning("\nConfidence Interval procedure not developed for lme, returning point estimates only")
  output = list(alpha = alpha, nboot = 0)
  # RESI point estimates
  output = c(output, resi_pe(model.full = model.full, vcovfunc = vcovfunc, vcov.args = vcov.args,
                             anova = anova, ...))
  output$boot.method = "nonparam"
  class(output) = "resi"
  return(output)
}

#' @describeIn resi RESI point and interval estimation for lmerMod models
#' @export
resi.lmerMod = function(model.full, alpha = 0.05, nboot = 1000, anova = TRUE,
                         vcovfunc = clubSandwich::vcovCR, vcov.args = list(), ...){
  warning("\nConfidence Interval procedure not developed for lmerMod, returning point estimates only")
  output = list(alpha = alpha, nboot = 0)
  output = c(output, resi_pe(model.full, vcovfunc = vcovfunc, vcov.args = vcov.args, anova = anova, ...)) # RESI point estimates
  output$boot.method = "nonparam"
  class(output) = "resi"
  return(output)
}
