#' Robust Effect Size Index (RESI) Point Estimation
#'
#' This function will estimate the robust effect size (RESI) from Vandekar, Tao, & Blume (2020).
#' The overall RESI is estimated via a Wald test. RESI is (optionally) estimated for each factor in coefficients-style table.
#' RESI is (optionally) estimated for each variable/interaction in an Anova-style table
#' for models with existing Anova methods. This function is the building block for the \code{\link{resi}} function.
#' @param model.full \code{lm, glm, nls, survreg, coxph, hurdle, zeroinfl, gee, geeglm} or \code{lme} model object.
#' @param model.reduced Fitted model object of same type as model.full. By default `NULL`; the same model as the full model but only having intercept.
#' @param data Data.frame or object coercible to data.frame of model.full data (required for some model types).
#' @param vcovfunc The variance estimator function for constructing the Wald test statistic. By default, sandwich::vcovHC (the robust (sandwich) variance estimator).
#' @param coefficients Logical, whether to produce a coefficients (summary) table with the RESI columns added. By default = `TRUE`.
#' @param anova Logical, whether to produce an Anova table with the RESI columns added. By default = `TRUE`.
#' @param overall Logical, whether to produce an overall Wald test comparing full to reduced model with RESI columns added. By default = `TRUE`.
#' @param Anova.args List, additional arguments to be passed to Anova function.
#' @param vcov.args List, additional arguments to be passed to vcovfunc.
#' @param unbiased Logical, whether to use the unbiased or alternative T/Z statistic to RESI conversion. By default, `TRUE`. See details.
#' @param waldtype Numeric, indicates which function to use for overall Wald test. 0 (default) = lmtest::waldtest Chi-square, 1 = lmtest::waldtest F, 2 = aod::wald.test
#' @param ... Ignored.
#' @importFrom aod wald.test
#' @importFrom car Anova
#' @importFrom lmtest waldtest
#' @importFrom sandwich vcovHC
#' @importFrom stats coef df.residual formula glm hatvalues pf predict quantile residuals update vcov
#' @importFrom utils methods
#' @export
#' @details The Robust Effect Size Index (RESI) is an effect size measure based on M-estimators.
#' This function is called by \code{\link{resi}} a specified number of times to
#' form bootstrapped confidence intervals. Called by itself, this function will
#' only calculate point estimates.
#'
#' The RESI, denoted as S, is applicable across many model types. It is a unitless
#' index and can be easily be compared across models. The RESI can also be
#' converted to Cohen's \emph{d} (\code{\link{S2d}}) under model homoskedasticity.
#'
#' The RESI is related to the non-centrality parameter
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
#' For GEE (\code{geeglm}) models, a longitudinal RESI (L-RESI) and a cross-sectional,
#' per-measurement RESI (CS-RESI) is estimated. The longitudinal RESI takes the
#' specified clustering into account, while the cross-sectional RESI is estimated
#' using a model where each measurement is its own cluster.
#'
#' @return Returns a list containing RESI point estimates
#' @examples
#' # This function produces point estimates for the RESI. The resi function will
#' # provide the same point estimates but adds confidence intervals. See resi for
#' # more detailed examples.
#'
#' ## resi_pe for a linear model
#' # fit linear model
#' mod <- lm(charges ~ region * age + bmi + sex, data = RESI::insurance)
#' # run resi_pe on the model
#' resi_pe(mod)
#'
#' # if you want to have RESI estimates in the coefficient table that are equal in absolute
#' # value to those in the Anova table (except for those with >1 df and/or included in other
#' # interaction terms), you can specify unbiased = FALSE to use the alternate conversion.
#' resi_pe(mod, unbiased = FALSE)
#' @references Vandekar S, Tao R, Blume J. A Robust Effect Size Index. \emph{Psychometrika}. 2020 Mar;85(1):232-246. doi: 10.1007/s11336-020-09698-2.

resi_pe = function(model.full, ...){
  UseMethod("resi_pe")
}

#' @describeIn resi_pe RESI point estimation
#' @export
resi_pe.default = function(model.full, model.reduced = NULL, data, anova = TRUE,
                           coefficients = TRUE, overall = TRUE, vcovfunc = sandwich::vcovHC,
                           Anova.args = list(), vcov.args = list(), unbiased = TRUE,
                           waldtype = 0, ...){

  # check for supported model type
  if(!(any(class(model.full) %in%
           sapply(as.character(utils::methods(RESI::resi_pe)), function(x) substr(x, 9, nchar(x)))))){
    warning("Model type not implemented in RESI package, attempting default")
  }

  dots = list(...)
  if (missing(data)){
    if (!is.null(model.full$model)){
      data = model.full$model
    } else{
      stop("\nData argument needed for model type")
    }
  }
  else{
    data = as.data.frame(data)
  }

  if (is.null(model.reduced) & overall){
    if(!("skip.red.pe" %in% names(dots))){
      form.reduced = as.formula(paste(format(formula(model.full)[[2]]), "~ 1"))
      if ((form.reduced == formula(model.full))){
        overall = FALSE} else{
          if (is.null(model.full$na.action)){
            model.reduced <- try(update(model.full, formula = form.reduced,
                                      data = data), silent = T)
          }
          else{
            if (!is.null(model.full$model)){
              model.reduced <- try(update(model.full, formula = form.reduced,
                                      data = model.full$model), silent = T)
            } else{
              model.reduced <- try(update(model.full, formula = form.reduced,
                                      data = data[which(!(1:nrow(data)%in% model.full$na.action)),]),
                                   silent = T)
            }}
          if(inherits(model.reduced, "try-error")){
            warning("Fitting intercept-only model failed. No overall test computed")
            overall = FALSE}
          }}}

  # dealing with additional vcov args
  if (length(vcov.args) == 0){
    vcovfunc2 <- vcovfunc
  } else{
    vcov.args = c(formals(vcovfunc)[1], vcov.args)
    vcovfunc2 <- function(x) {
      vcov.args[[1]] = x
      do.call(vcovfunc, vcov.args)}
  }

  # covariance matrix
  var.names = names(coef(model.full))
  tryCatch(vcovmat <- vcovfunc2(model.full), error = function(e){
    stop("\nComputing covariance failed. Trying running with different function for vcovfunc")
  })
  tryCatch(rownames(vcovmat) <- colnames(vcovmat) <- var.names, error = function(e){})


  # set up output
  if ((is.null(model.reduced) & !("skip.red.pe" %in% names(dots))) |
      inherits(model.reduced, "try-error")){
    output = list(model.full = list(call = model.full$call, formula = formula(model.full)),
                  model.reduced = NULL)
    overall = FALSE
  } else{
  output = list(model.full = list(call = model.full$call, formula = formula(model.full)),
                 model.reduced = list(call = model.reduced$call, formula = unlist(formula(model.reduced))))
  }
  names.est = c()
  # overall
  if (overall){
    # check if full and reduced model are the same
    if (!(is.null(model.reduced))){
      if (model.full$call == model.reduced$call){
        stop("\nFull and reduced model are identical")
      }}

    if (waldtype %in% c(0,1)){
      # 0 is Chisq with lmtest
      # 1 is F with lmtest
      wtype = ifelse(waldtype == 0, "Chisq", "F")
      tryCatch(overall.tab <- lmtest::waldtest(model.reduced, model.full, vcov = vcovfunc2,
                                               test = wtype), error = function(e){
                                                 stop("Overall Wald test failed, try running with overall = FALSE")})

      stats = overall.tab[,wtype][2]
      overall.df = overall.tab$Df[2]
      res.df = overall.tab$Res.Df[2]
    }

    if (waldtype == 2){
      # 2 is aod (nls, coxph)
      tryCatch(overall.tab <- aod::wald.test(vcovmat, coef(model.full),
                                             Terms = 1:length(coef(model.full)))$result$chi2,
               error = function(e){
                 stop("Overall Wald test failed, try running with overall = FALSE")})
      stats = overall.tab["chi2"]
      overall.df = overall.tab["df"]
      res.df = nrow(data) - overall.df
      overall.tab = t(as.data.frame(overall.tab))
    }
    overall.resi.hat = ifelse(waldtype == 1, f2S(stats, overall.df, res.df, nrow(data)),
                              chisq2S(stats, overall.df, nrow(data)))
    if (nrow(overall.tab) == 1){
      overall.tab = cbind(overall.tab, RESI = overall.resi.hat)
    } else{
      overall.tab[nrow(overall.tab), "RESI"] = overall.resi.hat
    }
    output$estimates = overall.resi.hat
    output$overall = overall.tab
    names.est = "Overall"
    names(output$estimates) = names.est
  }

  # coefficients table (z or t statistics)
  if (coefficients){
    if ("coeftest_df" %in% names(dots)){
      tryCatch(coefficients.tab <- lmtest::coeftest(model.full, vcov. = vcovmat, df = dots$coeftest_df),
               error = function(e){
                 stop("Coefficients table failed, try running with coefficients = FALSE")})
      } else{
    tryCatch(coefficients.tab <- lmtest::coeftest(model.full, vcov. = vcovmat), error = function(e){
      stop("Overall Wald test failed, try running with overall = FALSE")})}

    torZ = ifelse("z value" %in% colnames(coefficients.tab), "z", "t")
    coefficients.df = data.frame(coefficients.tab[,"Estimate"], coefficients.tab[,"Std. Error"],
                                 coefficients.tab[,paste(torZ, "value")],
                                 coefficients.tab[,paste("Pr(>|", torZ, "|)", sep = "")],
                                 row.names = rownames(coefficients.tab))
    colnames(coefficients.df) = colnames(coefficients.tab)

    if (torZ == "z"){
      coefficients.df[,"RESI"] = suppressWarnings(z2S(coefficients.df[,3], nrow(data), unbiased))
    } else{
      if(!overall){
        res.df = model.full$df.residual
      }
      coefficients.df[,"RESI"] = suppressWarnings(t2S(coefficients.df[,"t value"],
                                                      res.df, nrow(data), unbiased))
    }
    output$coefficients = coefficients.df
    output$estimates = c(output$estimates, coefficients.df$RESI)
    names.est = c(names.est, rownames(coefficients.df))
    names(output$estimates) = names.est
  }

  if (anova){
    if ("anovaF" %in% names(dots)){
      tryCatch(suppressMessages(anova.tab <- do.call(car::Anova, c(list(mod = model.full,
                                              vcov. = vcovmat), Anova.args))),
               error = function(e){
                 stop("car::Anova failed. Try rerunning with anova = FALSE")})
      anova.tab = anova.tab[which(rownames(anova.tab) != "Residuals"),]
      anova.tab[,"RESI"] = f2S(anova.tab[,"F"], anova.tab[,"Df"], res.df, nrow(data))
    } else {
      tryCatch(suppressMessages(anova.tab <- do.call(car::Anova, c(list(mod = model.full,
                                             test.statistic = "Wald",
                                             vcov. = vcovmat), Anova.args))),
               error = function(e){
          stop("car::Anova failed. Try rerunning with anova = FALSE")})
      anova.tab = anova.tab[which(rownames(anova.tab) != "Residuals"),]
      anova.tab[,"RESI"] = chisq2S(anova.tab[,"Chisq"], anova.tab[,"Df"], nrow(data))
    }

    output$anova = anova.tab
    names.est = names(output$estimates)
    output$estimates = c(output$estimates, anova.tab$RESI)
    names.est = c(names.est, rownames(anova.tab))
    names(output$estimates) = names.est
    class(output$anova) = c("anova_resi", class(output$anova))
  }

  if(identical(vcovfunc, stats::vcov)){
    output$naive.var = TRUE
  }
  else{
    output$naive.var = FALSE
  }

  class(output) = c("resi", "list")

  return(output)
}


#' @describeIn resi_pe RESI point estimation for generalized linear models
#' @export
resi_pe.glm = resi_pe.default

#' @describeIn resi_pe RESI point estimation for linear models
#' @export
resi_pe.lm = function(model.full, model.reduced = NULL, data, anova = TRUE,
                      coefficients = TRUE, vcovfunc = sandwich::vcovHC, Anova.args = list(),
                      vcov.args = list(), unbiased = TRUE,
                      overall = TRUE, ...){

  output = resi_pe.default(model.full = model.full, model.reduced = model.reduced,
                           data = data, anova = anova, coefficients = coefficients,
                           vcovfunc = vcovfunc, Anova.args = Anova.args,
                           vcov.args = vcov.args, unbiased = unbiased,
                           overall = overall, waldtype = 1, anovaF = TRUE, ...)
  return(output)
}



#' @describeIn resi_pe RESI point estimation for nonlinear least squares models
#' @export
resi_pe.nls = function(model.full, model.reduced = NULL, data, coefficients = TRUE,
                       anova = FALSE, vcovfunc = r_nlshc, vcov.args = list(),
                       unbiased = TRUE, overall = TRUE, ...){
  if (missing(data) | is.null(data)){
    stop("\nData argument is required for nls model")
  }
  data = as.data.frame(data)

  # currently not accepting other reduced models for nls
  if (!is.null(model.reduced)){
    warning("Reduced model argument ignored for nls model")
  }

  if (identical(vcovfunc, sandwich::vcovHC)){
    vcovfunc = r_nlshc
    warning("Sandwich vcov function not applicable for nls model type, vcovfunc set to regtools::nlshc")
  }

  output = resi_pe.default(model.full = model.full, model.reduced = NULL, data = data,
                           anova = FALSE, coefficients = coefficients, vcovfunc = vcovfunc,
                           vcov.args = vcov.args, unbiased = unbiased, overall = overall,
                           waldtype = 2, skip.red.pe = TRUE, ...)

  return(output)
}

#' @describeIn resi_pe RESI point estimation for survreg
#' @export
resi_pe.survreg = function(model.full, model.reduced = NULL, data, anova = TRUE,
                           coefficients = TRUE, vcovfunc = vcov, Anova.args = list(),
                           unbiased = TRUE, overall = TRUE, ...){

  if (missing(data)){
    stop("\nData argument is required for survreg model")
  }
  else{
    data = as.data.frame(data)
  }

  if (!identical(vcovfunc, vcov)){
    warning("vcovfunc argument ignored for survreg objects")
  }

  output = resi_pe.default(model.full = model.full, model.reduced = model.reduced,
                           data = data, anova = anova, coefficients = coefficients,
                           vcovfunc = vcov, Anova.args = Anova.args, unbiased = unbiased,
                           overall = overall, waldtype = 0, ...)

  if(is.null(model.full$naive.var)){
    output$naive.var = TRUE
  }
  else{
    output$naive.var = FALSE
  }

  return(output)
}

#' @describeIn resi_pe RESI point estimation for coxph models
#' @export
resi_pe.coxph = function(model.full, model.reduced = NULL, data, anova = TRUE,
                         coefficients = TRUE, vcovfunc = vcov, Anova.args = list(),
                         unbiased = TRUE, overall = TRUE, ...){
  if (missing(data)){
    stop("\nData argument is required for coxph model")
  }
  else{
    data = as.data.frame(data)
  }

  if (!identical(vcovfunc, vcov)){
    warning("vcovfunc argument ignored for coxph objects")
  }

  # currently not accepting other reduced models for coxph
  if (!is.null(model.reduced)){
    warning("Reduced model argument ignored for coxph model")
  }

  output = resi_pe.default(model.full = model.full, model.reduced = NULL,
                           data = data, anova = anova, coefficients = coefficients,
                           vcovfunc = vcov, Anova.args = Anova.args, unbiased = unbiased,
                           overall = overall, waldtype = 2, ...)

  if(is.null(model.full$naive.var)){
    output$naive.var = TRUE
  }
  else{
    output$naive.var = FALSE
  }

  return(output)
}

#' @describeIn resi_pe RESI point estimation for hurdle models
#' @export
resi_pe.hurdle = function(model.full, model.reduced = NULL, data, coefficients = TRUE,
                          anova = TRUE, vcovfunc = sandwich::sandwich, vcov.args = list(),
                          unbiased = TRUE, overall = TRUE, ...){

  if (identical(vcovfunc, sandwich::vcovHC)){
    vcovfunc = sandwich::sandwich
  }

  output = resi_pe.default(model.full = model.full, model.reduced = model.reduced,
                           data = data, anova = FALSE, coefficients = coefficients,
                           vcovfunc = vcovfunc, vcov.args = vcov.args, unbiased = unbiased,
                           overall = overall, waldtype = 0,
                           coeftest_df = 0, ...)

  return(output)
}

#' @describeIn resi_pe RESI point estimation for zeroinfl models
#' @export
resi_pe.zeroinfl = resi_pe.hurdle

#' @describeIn resi_pe RESI point estimation for geeglm object
#' @export
resi_pe.geeglm <- function(model.full, model.reduced = NULL, data, anova = TRUE,
                           coefficients = TRUE, overall = TRUE, unbiased = TRUE, ...){

  if (missing(data)){
    data = model.full$data
  }

  # sample size
  N = length(summary(model.full)$clusz)
  # total num of observations
  tot_obs = nrow(data)

  # num of observations from each subject
  n_i = table(model.full$id)
  n_i = rep(n_i, times = n_i)
  # weight in independent model
  w = 1 / n_i
  data$w = w
  # sample size
  N = length(unique(model.full$id))
  # total num of observations
  tot_obs = nrow(data)

  # model form
  form = formula(model.full)
  # independence model
  mod_indg = suppressWarnings(glm(formula = form, family = model.full$family, data = data,
                                  weights = w, contrasts = model.full$contrasts))
  mod_indg$residuals = mod_indg$residuals / sqrt(mod_indg$weights)
  # the var-cov matrix estimate from the independence model
  # Note: this is the estimate for Cov[(\hat{\beta}_ind - \beta_0)]
  cov_ind = sandwich::vcovHC(mod_indg, type = "HC0")

  # make copy of model.full, replace vbeta with independence
  mod_ind = model.full
  mod_ind$geese$vbeta = cov_ind

  output = list(model.full = list(call = model.full$call, formula = form),
                model.reduced = NULL)
  names.est = c()
  if (overall){
  # reduced model
  if (is.null(model.reduced)){
    form.reduced = as.formula(paste(format(formula(model.full)[[2]]), "~ 1"))
    model.reduced = update(model.full, formula = form.reduced, data = data)
  } else{
    form.reduced = formula(model.reduced)
  }

  # overall
  # longitudinal
  overall.tab = anova(model.full, model.reduced)
  overall.tab[,"L-RESI"] = chisq2S(overall.tab[,"X2"], overall.tab[,"Df"], N)

  output$model.reduced = list(call = model.reduced$call, formula = form.reduced)
  output$estimates = c(overall.tab[,"L-RESI"])
  output$overall = overall.tab
  names.est = c("Overall-L")
  names(output$estimates) = names.est
  }

  # longitudinal RESI
  # coefficients (z statistics)
  if (coefficients) {
    # longitudinal resi
    # geeglm uses robust estimate by default, do not need to specify vcov
    coefficients.tab <- lmtest::coeftest(model.full)
    coefficients.df = data.frame(coefficients.tab[,"Estimate"],
                                 coefficients.tab[,"Std. Error"],
                                 coefficients.tab[,"z value"],
                                 coefficients.tab[,"Pr(>|z|)"],
                                 row.names = rownames(coefficients.tab))
    colnames(coefficients.df) = colnames(coefficients.tab)
    coefficients.df[,"L-RESI"] = z2S(coefficients.df[,"z value"], N, unbiased)

    # CS RESI
    # M: use independence mod
    coefficients.tabcs <- lmtest::coeftest(mod_ind, vcov. = cov_ind)
    z_cs = coefficients.tabcs[,"z value"]
    coefficients.df[,"CS-RESI"] = z2S(z_cs, N, unbiased)

    output$coefficients = coefficients.df
    output$estimates = c(output$estimates, coefficients.df$`L-RESI`,
                         coefficients.df[,"CS-RESI"])
    names.est = c(names.est, rep(rownames(coefficients.df), 2))
    names(output$estimates) = names.est
  }

  # anova
  if (anova) {
    # longitudinal RESI
    anova.tab <- anova(model.full)
    anova.tab[,"L-RESI"] = chisq2S(anova.tab[,"X2"], anova.tab[,"Df"], N)

    # CS RESI
    # M: use anova() on new mod with substituted variance
    anova.tabcs <- suppressMessages(car::Anova(mod_indg, vcov. = cov_ind,
                                               test.statistic = "Wald"))
    anova.tab[,"CS-RESI"] = chisq2S(anova.tabcs[,"Chisq"], anova.tabcs[,"Df"], N)

    output$anova = anova.tab
    output$estimates = c(output$estimates, anova.tab$`L-RESI`,
                         anova.tab[,"CS-RESI"])
    names.est = c(names.est, rep(rownames(anova.tab), 2))
    names(output$estimates) = names.est
    class(output$anova) = c("anova_resi", class(output$anova))
  }

  output$naive.var = FALSE
  class(output) = c("resi", "list")
  return(output)
}

#' @describeIn resi_pe RESI point estimation for gee object
#' @export
resi_pe.gee <- function(model.full, data, unbiased = TRUE, ...){
  # sample size
  N = length(unique(model.full$id))
  # total num of observations
  tot_obs = nrow(data)

  # num of observations from each subject
  n_i = table(model.full$id)
  n_i = rep(n_i, times = n_i)
  # weight in independent model
  w = 1 / n_i
  data$w = w

  mod_indg = suppressWarnings(glm(formula = formula(model.full), family = model.full$family, data = data,
                                  contrasts = model.full$contrasts, weights = w))
  mod_indg$residuals = mod_indg$residuals / sqrt(mod_indg$weights)
  # the var-cov matrix estimate from the independence model
  # Note: this is the estimate for Cov[(\hat{\beta}_ind - \beta_0)]
  cov_ind = sandwich::vcovHC(mod_indg, type = "HC0")
  # cov_ind = cov_ind * tot_obs/N
  # make copy of full model and substitute independence variance
  mod_ind = model.full
  mod_ind$robust.variance = cov_ind

  output <- list(model.full = list(call = model.full$call, formula = formula(model.full)),
                 estimates = c())
  names.est = c()



  # longitudinal RESI
  # coefficients (z statistics)
  # uses different function to calculate than other classes
  coefficients.tab <- summary(model.full)$coefficients
  coefficients.df = data.frame(coefficients.tab[,'Estimate'],
                               coefficients.tab[,'Robust S.E.'],
                               coefficients.tab[,'Robust z'],
                               row.names = rownames(coefficients.tab))
  colnames(coefficients.df) = c('Estimate', 'Std. Error', 'z value')
  coefficients.df[,'L-RESI'] = z2S(coefficients.df[,'z value'], N, unbiased)

  # CS RESI
  coefficients.tabcs = summary(mod_ind)$coefficients
  z_cs = coefficients.tabcs[,'Robust z']
  coefficients.df[,'CS-RESI'] = z2S(z_cs, N, unbiased)

  output$coefficients = coefficients.df
  output$estimates = c(coefficients.df$`L-RESI`, coefficients.df[,"CS-RESI"])
  names.est = rep(rownames(coefficients.df),2)
  names(output$estimates) = names.est

  output$naive.var = FALSE
  class(output) = c("resi", "list")
  return(output)
}

#' @describeIn resi_pe RESI point estimation for lme object
#' @importFrom clubSandwich vcovCR
#' @export
resi_pe.lme <- function(model.full, anova = TRUE, vcovfunc = clubSandwich::vcovCR,
                        Anova.args = list(), vcov.args = list(), ...){
  x = as.matrix(summary(model.full)$tTable)
  #sample size
  N = summary(model.full)$dims$ngrps[1]
  # robust se
  if (length(vcov.args) == 0){
    vcov.args = c(formals(vcovfunc)[1], list(type = "CR3"))
  } else{
    vcov.args = c(formals(vcovfunc)[1], vcov.args)
    if (!("type" %in% names(vcov.args))){
      vcov.args = c(vcov.args, list(type = "CR3"))
    }}
  vcovfunc2 <- function(x) {
    vcov.args[[1]] = x
    do.call(vcovfunc, vcov.args)}

  robust.var = diag(vcovfunc2(model.full))
  robust.se = sqrt(robust.var)
  coefficients.df = as.data.frame(cbind(x, 'Robust.SE' = robust.se,
                                        'Robust Wald' = (x[, 'Value']^2/robust.var),
                                        RESI = RESI::chisq2S(x[, 'Value']^2/robust.var, 1, N)))
  output = list(model.full = list(call = model.full$call, formula = formula(model.full)),
                estimates = coefficients.df$RESI, coefficients = coefficients.df)
  names(output$estimates) = rownames(coefficients.df)

  # Anova table (Chi sq statistics)
  if (anova){
    suppressMessages(anova.tab <- do.call(car::Anova, c(list(mod = model.full, vcov. = vcovfunc2), Anova.args)))
    anova.tab[,'RESI'] = chisq2S(anova.tab[,'Chisq'], anova.tab[,'Df'], N)
    names.est = names(output$estimates)
    output$estimates = c(output$estimates, anova.tab$RESI)
    names.est = c(names.est, rownames(anova.tab))
    names(output$estimates) = names.est
    output$anova = anova.tab
    class(output$anova) = c("anova_resi", class(output$anova))
  }
  if(identical(vcovfunc, stats::vcov)){
    output$naive.var = TRUE
  }
  else{
    output$naive.var = FALSE
  }
  class(output) = c("resi", "list")
  return(output)
}

#' @describeIn resi_pe RESI point estimation for lmerMod object
#' @export
resi_pe.lmerMod <- function(model.full, anova = TRUE, vcovfunc = clubSandwich::vcovCR,
                            Anova.args = list(), vcov.args = list(), ...){
  x = as.matrix(summary(model.full)$coefficients)
  #sample size
  N = summary(model.full)$ngrps

  # robust se
  if (identical(vcovfunc, stats::vcov)) {
    coefficients.tab = cbind(x, 'Wald' = x[, 't value']^2, RESI = RESI::chisq2S(x[, 't value']^2, 1, N))
    output = list(estimates = coefficients.tab[,"RESI"], coefficients = coefficients.tab, naive.var = TRUE)
    vcovfunc2 = vcovfunc
  } else {
    if (length(vcov.args) == 0){
      if (identical(vcovfunc, clubSandwich::vcovCR)){
        vcov.args = c(formals(vcovfunc)[1], list(type = "CR3"))
      }} else{
        if (identical(vcovfunc, clubSandwich::vcovCR)){
          if (!("type" %in% names(vcov.args))){
            vcov.args = c(vcov.args, list(type = "CR3"))}}
        vcov.args = c(formals(vcovfunc)[1], vcov.args)}

    vcovfunc2 <- function(x) {
      vcov.args[[1]] = x
      do.call(vcovfunc, vcov.args)}
    robust_var = diag(vcovfunc2(model.full))
    robust_se = sqrt(robust_var)
    coefficients.tab = cbind(x, 'Robust.SE' = robust_se,
                             'Robust Wald' = (x[, 'Estimate']^2/robust_var),
                             RESI = RESI::chisq2S(x[, 'Estimate']^2/robust_var, 1, N))
    # model full is not the same process as others, so print method is missing the call
    output = list(estimates = coefficients.tab[,"RESI"], coefficients = coefficients.tab,
                  naive.var = FALSE)
  }
  names(output$estimates) = rownames(coefficients.tab)
  # Anova table (Chi sq statistics)
  if (anova){
    suppressMessages(anova.tab <- do.call(car::Anova,
                                          c(list(mod = model.full, vcov. = vcovfunc2), Anova.args)))
    anova.tab[,'RESI'] = chisq2S(anova.tab[,'Chisq'], anova.tab[,'Df'], N)
    names.est = names(output$estimates)
    output$estimates = c(output$estimates, anova.tab$RESI)
    names.est = c(names.est, rownames(anova.tab))
    names(output$estimates) = names.est
    output$anova = anova.tab
    class(output$anova) = c("anova_resi", class(output$anova))
  }
  output$model.full = list(call = model.full@call, formula = formula(model.full))
  class(output) = c("resi", "list")
  return(output)
}
