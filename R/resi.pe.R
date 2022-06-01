#' Robust Effect Size index (RESI) point estimation
#' @param model.full the full model.
#' @param model.reduced the reduced model to compare with the full model. By default `NULL`, it's the same model as the full model but only having intercept.
#' @param data optional data frame or object coercible to data frame of model.full data (required for some model types)
#' @param anova whether to produce an Anova table with the RESI columns added. By default = `TRUE`
#' @param summary whether to produce a summary table with the RESI columns added. By default = `TRUE`
#' @param vcovfunc the variance estimator function for constructing the Wald test statistic. By default, sandwich::vcovHC (the robust (sandwich) variance estimator)
#' @param ... Other arguments to be passed to Anova function
#' @importFrom aod wald.test
#' @importFrom car Anova
#' @importFrom lmtest waldtest
#' @importFrom regtools nlshc
#' @importFrom sandwich vcovHC
#' @importFrom stats coef formula glm hatvalues pf predict quantile residuals update vcov
#' @export
#' @return Returns a list of RESI point estimates

resi.pe <- function(model.full, ...){
  UseMethod("resi.pe", model.full)
}

resi.pe.glm <- function(model.full, model.reduced = NULL, data, anova = TRUE,
                    summary = TRUE, vcovfunc = sandwich::vcovHC, ...){
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

  # Anova table (Chi sq statistics)
  if (anova){
    suppressMessages(anova.tab <- car::Anova(model.full, test.statistic = 'Wald', vcov. = vcovfunc, ...))
    anova.tab[,'RESI'] = chisq2S(anova.tab[,'Chisq'], anova.tab[,'Df'], model.full$df.residual)
    output$anova = anova.tab
    output$estimates = c(output$estimates, anova.tab$RESI)
    names.est = c(names.est, rownames(anova.tab))
    names(output$estimates) = names.est
  }

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


resi.pe.lm <- function(model.full, model.reduced = NULL, data, anova = TRUE,
                       summary = TRUE, vcovfunc = sandwich::vcovHC, ...){
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

  # Anova table (F statistics)
  if (anova){
    suppressMessages(anova.tab <- car::Anova(model.full, vcov. = vcovfunc, ...))
    anova.tab[,'RESI'] = f2S(anova.tab[,'F'], anova.tab[,'Df'], res.df)
    output$anova = anova.tab
    output$estimates = c(output$estimates, anova.tab$RESI)
    names.est = c(names.est, rownames(anova.tab))
    names(output$estimates) = names.est
  }

  # summary table (t statistics)
  if (summary){
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
  return(output)
}

# add nls, survreg, coxph methods
