#' Robust Effect Size index (RESI) Point Estimation
#' @param model.full \code{lm, glm, nls, survreg, coxph, hurdle, zeroinfl, gee, geeglm} or \code{lme} model object.
#' @param model.reduced Fitted model object of same type as model.full. By default `NULL`; the same model as the full model but only having intercept.
#' @param data Data.frame or object coercible to data.frame of model.full data (required for some model types).
#' @param vcovfunc The variance estimator function for constructing the Wald test statistic. By default, sandwich::vcovHC (the robust (sandwich) variance estimator).
#' @param summary Logical, whether to produce a summary table with the RESI columns added. By default = `TRUE`.
#' @param anova Logical, whether to produce an Anova table with the RESI columns added. By default = `TRUE`.
#' @param ... Other arguments to be passed to Anova function.
#' @importFrom aod wald.test
#' @importFrom car Anova
#' @importFrom lmtest waldtest
#' @importFrom regtools nlshc
#' @importFrom sandwich vcovHC
#' @importFrom stats coef formula glm hatvalues pf predict quantile residuals update vcov
#' @export
#' @return Returns a list containing RESI point estimates

resi_pe <- function(model.full, ...){
  UseMethod("resi_pe")
}

#' @describeIn resi_pe RESI point estimation
#' @export
resi_pe.default <- function(model.full, model.reduced = NULL, data,
                    summary = TRUE, vcovfunc = sandwich::vcovHC){
  if (missing(data)){
      data = model.full$model
  }

  if (is.null(model.reduced)){
    form.reduced = as.formula(paste(format(formula(model.full)[[2]]), "~ 1"))
    model.reduced <- update(model.full, formula = form.reduced, data = data)
  }

  # overall
  wald.test = lmtest::waldtest(model.reduced, model.full, vcov = vcovfunc, test = 'Chisq')
  stats = wald.test$Chisq[2]
  overall.df = wald.test$Df[2]
  res.df = wald.test$Res.Df[2]
  overall.resi.hat = chisq2S(stats, overall.df, res.df)
  wald.test[2, 'RESI'] = overall.resi.hat
  output <- list(estimates = overall.resi.hat, overall = wald.test)
  names.est = "Overall"
  names(output$estimates) = names.est

  # summary table (z statistics)
  if (summary){
    summary.tab <- lmtest::coeftest(model.full, vcov. = vcovfunc)
    summary.df = data.frame(summary.tab[,'Estimate'], summary.tab[,'Std. Error'],
                             summary.tab[,'z value'], summary.tab[,'Pr(>|z|)'], row.names = rownames(summary.tab))
    colnames(summary.df) = colnames(summary.tab)
    summary.df[,'RESI'] = z2S(summary.df[,'z value'], nrow(data))
    output$summary = summary.df
    output$estimates = c(output$estimates, summary.df$RESI)
    names.est = c(names.est, rownames(summary.df))
    names(output$estimates) = names.est
  }
  return(output)
}

#' @describeIn resi_pe RESI point estimation for generalized linear models
#' @export
resi_pe.glm <- function(model.full, model.reduced = NULL, data, anova = TRUE,
                            summary = TRUE, vcovfunc = sandwich::vcovHC, ...){
  output <- resi_pe.default(model.full, model.reduced, data, summary, vcovfunc)

  # Anova table (Chi sq statistics)
  if (anova){
    suppressMessages(anova.tab <- car::Anova(model.full, test.statistic = 'Wald', vcov. = vcovfunc, ...))
    anova.tab[,'RESI'] = chisq2S(anova.tab[,'Chisq'], anova.tab[,'Df'], model.full$df.residual)
    output$anova = anova.tab
    names.est = names(output$estimates)
    output$estimates = c(output$estimates, anova.tab$RESI)
    names.est = c(names.est, rownames(anova.tab))
    names(output$estimates) = names.est
  }
  return(output)
}

#' @describeIn resi_pe RESI point estimation for linear models
#' @export
resi_pe.lm <- function(model.full, model.reduced = NULL, data, anova = TRUE,
                       summary = TRUE, vcovfunc = sandwich::vcovHC, ...){
  # data required with splines
  if (missing(data)){
    data = model.full$model
  }

  if (is.null(model.reduced)){
    form.reduced = as.formula(paste(format(formula(model.full)[[2]]), "~ 1"))
    model.reduced <- update(model.full, formula = form.reduced, data = data)
  }

  # overall
  wald.test = lmtest::waldtest(model.reduced, model.full, vcov = vcovfunc, test = 'Chisq')
  stats = wald.test$Chisq[2]
  overall.df = wald.test$Df[2]
  res.df = wald.test$Res.Df[2]
  overall.resi.hat = chisq2S(stats, overall.df, res.df)
  wald.test[2, 'RESI'] = overall.resi.hat
  output <- list(estimates = overall.resi.hat, overall = wald.test)
  names.est = "Overall"
  names(output$estimates) = names.est

  # summary table (t statistics)
  if (summary){
    ## KK: are we sure this function always returns the t-statistics?
    ## KK: will it change with the different model type passed by `model.full`?
    ## MJ: for lm it will be t unless you manually set the df argument for coeftest to Inf or 0
    summary.tab <- lmtest::coeftest(model.full, vcov. = vcovfunc)
    summary.df = data.frame(summary.tab[,'Estimate'], summary.tab[,'Std. Error'],
                            summary.tab[,'t value'], summary.tab[,'Pr(>|t|)'], row.names = rownames(summary.tab))
    colnames(summary.df) = colnames(summary.tab)
    summary.df[,'RESI'] = t2S(summary.df[,'t value'], model.full$df.residual)
    output$summary = summary.df
    output$estimates = c(output$estimates, summary.df$RESI)
    names.est = c(names.est, rownames(summary.df))
    names(output$estimates) = names.est
  }

  # Anova table (F statistics)
  if (anova){
    suppressMessages(anova.tab <- car::Anova(model.full, vcov. = vcovfunc, ...))
    anova.tab[,'RESI'] = f2S(anova.tab[,'F'], anova.tab[,'Df'], res.df)
    output$anova = anova.tab
    output$estimates = c(output$estimates, anova.tab$RESI)
    names.est = c(names.est, rownames(anova.tab))
    names(output$estimates) = names.est
  }
  return(output)
}

#' @describeIn resi_pe RESI point estimation for nonlinear least squares models
#' @export
resi_pe.nls <- function(model.full, model.reduced = NULL, data,
                       summary = TRUE, vcovfunc = regtools::nlshc, ...){
  if (missing(data) | is.null(data)){
    stop('\n Data argument is required for nls model')
  }

  # currently not accepting other reduced models for nls
  if (!is.null(model.reduced)){
    warning('Reduced model argument ignored for nls model')
  }

  if (identical(vcovfunc, sandwich::vcovHC)){
    vcovfunc = regtools::nlshc
    warning("Sandwich vcov function not applicable for nls model type, vcovfunc set to regtools::nlshc")
  }

  var.names = names(coef(model.full))
  vcovmat = vcovfunc(model.full)
  rownames(vcovmat) = var.names
  colnames(vcovmat) = var.names

  # overall
  overall.tab = aod::wald.test(vcov(model.full), coef(model.full), Terms = 1:length(coef(model.full)))$result$chi2
  stats = overall.tab["chi2"]
  overall.df = overall.tab["df"]
  res.df = nrow(data) - overall.df
  overall.resi.hat = chisq2S(stats, overall.df, res.df)
  overall.tab['RESI'] = overall.resi.hat
  overall.tab = as.data.frame(t(overall.tab))
  rownames(overall.tab) = "Wald Test"
  output <- list(estimates = overall.resi.hat, overall = overall.tab)
  names.est = "Overall"
  names(output$estimates) = names.est

  # summary table (t statistics)
  if (summary){
    summary.tab <- lmtest::coeftest(model.full, vcov. = vcovmat)
    summary.df = data.frame(summary.tab[,'Estimate'], summary.tab[,'Std. Error'],
                            summary.tab[,'t value'], summary.tab[,'Pr(>|t|)'], row.names = rownames(summary.tab))
    colnames(summary.df) = colnames(summary.tab)
    summary.df[,'RESI'] = t2S(summary.df[,'t value'], res.df)
    output$summary = summary.df
    output$estimates = c(output$estimates, summary.df$RESI)
    names.est = c(names.est, rownames(summary.df))
    names(output$estimates) = names.est
  }
  return(output)
}

#' @describeIn resi_pe RESI point estimation for survreg
#' @export
resi_pe.survreg <- function(model.full, model.reduced = NULL, data, anova = TRUE,
                        summary = TRUE, vcovfunc = vcov, ...){
  if (missing(data)){
    stop('\n Data argument is required for survreg model')
  }

  if (!identical(vcovfunc, vcov)){
    warning("vcov function ignored for survreg objects")
  }

  output = resi_pe.default(model.full, model.reduced, data, summary, vcovfunc = vcov)

  # Anova table (Chi sq statistics)
  if (anova){
    suppressMessages(anova.tab <- car::Anova(model.full, ...))
    anova.tab[,'RESI'] = chisq2S(anova.tab[,'LR Chisq'], anova.tab[,'Df'], model.full$df.residual)
    output$anova = anova.tab
    names.est = names(output$estimates)
    output$estimates = c(output$estimates, anova.tab$RESI)
    names.est = c(names.est, rownames(anova.tab))
    names(output$estimates) = names.est
  }

  return(output)
}

#' @describeIn resi_pe RESI point estimation for hurdle and zeroinfl models
#' @export
resi_pe.zeroinfl <- resi_pe.hurdle <- function(model.full, model.reduced = NULL, data,
                           summary = TRUE, vcovfunc = sandwich::sandwich, ...){
  if (missing(data)){
    data = model.full$model
  }

  if (is.null(model.reduced)){
    form.reduced = as.formula(paste(format(formula(model.full)[[2]]), "~ 1"))
    model.reduced <- update(model.full, formula = form.reduced, data = data)
  }

  if (identical(vcovfunc, sandwich::vcovHC)){
    vcovfunc = sandwich::sandwich
  }

  # overall
  wald.test = lmtest::waldtest(model.reduced, model.full, vcov = vcovfunc, test = 'Chisq')
  stats = wald.test$Chisq[2]
  overall.df = wald.test$Df[2]
  res.df = wald.test$Res.Df[2]
  overall.resi.hat = chisq2S(stats, overall.df, res.df)
  wald.test[2, 'RESI'] = overall.resi.hat
  output <- list(estimates = overall.resi.hat, overall = wald.test)
  names.est = "Overall"
  names(output$estimates) = names.est

  # summary table (z statistics)
  if (summary){
    summary.tab <- lmtest::coeftest(model.full, vcov. = vcovfunc, df = 0)
    summary.df = data.frame(summary.tab[,'Estimate'], summary.tab[,'Std. Error'],
                            summary.tab[,'z value'], summary.tab[,'Pr(>|z|)'], row.names = rownames(summary.tab))
    colnames(summary.df) = colnames(summary.tab)
    summary.df[,'RESI'] = z2S(summary.df[,'z value'], nrow(data))
    output$summary = summary.df
    output$estimates = c(output$estimates, summary.df$RESI)
    names.est = c(names.est, rownames(summary.df))
    names(output$estimates) = names.est
  }

  return(output)
}

#' @describeIn resi_pe RESI point estimation for coxph models
#' @export
resi_pe.coxph <- function(model.full, model.reduced = NULL, data, anova = TRUE,
                       summary = TRUE, vcovfunc = vcov, ...){
  if (missing(data)){
    stop('\n Data argument is required for coxph model')
  }

  if (!identical(vcovfunc, vcov)){
    warning("vcov function ignored for coxph objects")
  }

  # currently not accepting other reduced models for coxph
  if (!is.null(model.reduced)){
    warning('Reduced model argument ignored for coxph model')
  }

  # overall
  overall.tab = aod::wald.test(vcov(model.full), coef(model.full), Terms = 1:length(coef(model.full)))$result$chi2
  stats = overall.tab["chi2"]
  overall.df = overall.tab["df"]
  res.df = nrow(data) - overall.df
  overall.resi.hat = chisq2S(stats, overall.df, res.df)
  overall.tab['RESI'] = overall.resi.hat
  overall.tab = as.data.frame(t(overall.tab))
  rownames(overall.tab) = "Wald Test"
  ## the output object
  output <- list(estimates = overall.resi.hat, overall = overall.tab)
  names.est = "Overall"
  names(output$estimates) = names.est

  if (summary){
    summary.tab <- lmtest::coeftest(model.full, vcov. = vcov)
    summary.df = data.frame(summary.tab[,'Estimate'], summary.tab[,'Std. Error'],
                            summary.tab[,'z value'], summary.tab[,'Pr(>|z|)'], row.names = rownames(summary.tab))
    colnames(summary.df) = colnames(summary.tab)
    summary.df[,'RESI'] = z2S(summary.df[,'z value'], nrow(data))
    output$summary = summary.df
    output$estimates = c(output$estimates, summary.df$RESI)
    names.est = c(names.est, rownames(summary.df))
    names(output$estimates) = names.est
  }

  # Anova table (Chi sq statistics)
  if (anova){
    suppressMessages(anova.tab <- car::Anova(model.full, test.statistic = 'Wald', ...))
    anova.tab[,'RESI'] = chisq2S(anova.tab[,'Chisq'], anova.tab[,'Df'], res.df)
    output$anova = anova.tab
    output$estimates = c(output$estimates, anova.tab$RESI)
    names.est = c(names.est, rownames(anova.tab))
    names(output$estimates) = names.est
  }

  return(output)
}


#' @describeIn resi_pe RESI point estimation for geeglm object
#' @export
resi_pe.geeglm <- function(model.full, ...){
  x = as.matrix(summary(model.full)$coefficients)
  #sample size
  N = length(summary(model.full)$clusz)
  output = cbind(x, RESI = RESI::chisq2S(x[, 'Wald'], 1, N))
  return(output)
}

#' @describeIn resi_pe RESI point estimation for gee object
#' @export
resi_pe.gee <- function(model.full, ...){
  x = as.matrix(summary(model.full)$coefficients)
  #sample size
  N = length(unique(model.full$id))
  output = cbind(x, RESI = RESI::chisq2S(x[, 'Robust z']^2, 1, N))
  return(output)
}

#' @describeIn resi_pe RESI point estimation for lme object
#' @importFrom clubSandwich vcovCR
#' @export
resi_pe.lme <- function(model.full, ...){
  x = as.matrix(summary(model.full)$tTable)
  #sample size
  N = summary(model.full)$dims$ngrps[1]
  # robust se
  robust.var = diag(clubSandwich::vcovCR(model.full, type = "CR3"))
  robust.se = sqrt(robust.var)
  output = cbind(x, 'Robust.SE' = robust.se, 'Robust Wald' = (x[, 'Value']^2/robust.var), RESI = RESI::chisq2S(x[, 'Value']^2/robust.var, 1, N))
  return(output)
}


