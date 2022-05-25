#' Robust effect size add-on for model fit objects
#' @export
calc_resi <- function(x, ...){
  UseMethod("calc_resi")
}

#' Robust effect size add-on for lm and glm
#' Effect size summary that can compliment your summary output. All conversions are based on the robust effect size index estimator.
#' @param object The model object.
#' @param object.reduced The reduced model to compare with the object while calculating the effect size. By default, `model.reduced = NULL`; if a reduced model is specified, only the effect size for the different parameters will be shown.
#' @param vcov. The variance covariance matrix to use. Defaults to one of the heteroskedasticity consistent estimator. Note that this differs from the default in lmtest.
#' @param ... Arguments passed to `vcov.`
#' @keywords lm effect size
#' @return Returns a summary of tests and the robust effect size index for linear models.
#' @details Stuff and stuff.
#' @importFrom sandwich vcovHC
#' @importFrom lmtest coeftest
#' @export
calc_resi.lm = function(object, object.reduced = NULL, vcov.= sandwich::vcovHC, ...){
  if (is.null(object.reduced)) {
    x = as.matrix(summary(object)$coefficients)
    se = sqrt(diag(vcov.(object, ...)))
    x = cbind('Estimate' = x[, 1], 's.e.' = se)
    x = cbind(x, 'Wald' = (x[, "Estimate"]/x[, 's.e.'])^2)
    x = cbind(x, 'p-value'= pchisq(x[, "Wald"], df = 1, lower.tail = FALSE))
    resi.tab = cbind(x, RESI = RESI::chisq2S(x[, 'Wald'], df = 1, object$df.residual))
  } else {

    # Overall (Wald) test stat
    wald.test = lmtest::waldtest(object.reduced, object, vcov = vcov., test = 'Chisq')
    stats = wald.test$Chisq[2]
    overall.df = wald.test$Df[2]
    res.df = wald.test$Res.Df[2]
    # overall RESI
    overall.resi.hat = RESI::chisq2S(stats, overall.df, res.df)

    # set up ANOES table
    anoes.tab = matrix(rep(NA, 2*4), nrow = 2, ncol = 4)
    anoes.tab[1, c(1, 2, 4) ] = c(stats, overall.df, overall.resi.hat)
    anoes.tab[2, 2] = res.df
    ## rename
    rownames(anoes.tab) = c("Tested", "Residual")
    colnames(anoes.tab) = c("Wald", "df", "p-val", "RESI")
    ## calculate p-values form Chi-sq dist
    anoes.tab[, 'p-val'] = pchisq(anoes.tab[, "Wald"], df = anoes.tab[, "df"], lower.tail = FALSE)
    resi.tab = anoes.tab
  }
  output = list(resi.tab = resi.tab,
                vcov. = vcov.,
                object = object,
                object.reduced = object.reduced)
  return(output)
}

#' Robust effect size add-on for glm
#'
#' Effect size summary that can compliment your summary output. All conversions are based on the robust effect size index estimator.
# #' @param object The model object.
# #' @param vcov. The variance covariance matrix to use. Defaults to one of the heteroskedasticity consistent estimator. Note that this differs from the default in lmtest.
# #' @param ... Arguments passed to coeftest
# #' @keywords effect size glm
# #' @return Returns a summary of tests and the robust effect size index for generalized linear models.
# #' @details Stuff and more stuff.
# #' @importFrom sandwich vcovHC
# #' @importFrom lmtest coeftest
# calc_resi.glm = function(object, object.reduced = NULL, vcov.=sandwich::vcovHC, ...){
#   if (is.null(object.reduced)){
#     x = as.matrix(summary(object)$coefficients)
#     robust.se = sqrt(diag(vcov.(object, ...)))
#     x = cbind('Estimate' = x[, 1], 'Robust s.e.' = robust.se)
#     x = cbind(x, 'Robust Wald' = (x[, "Estimate"]/x[, 'Robust s.e.'])^2)
#     x = cbind(x, 'p-value'= pchisq(x[, "Robust Wald"], df = 1, lower.tail = FALSE))
#     resi.tab = cbind(x, RESI = RESI::chisq2S(x[,'Robust Wald'], 1, object$df.residual))
#   } else {
#
#   }
#   output = list(resi.tab = resi.tab,
#                 vcov. = vcov.)
#   return(output)
# }

# Robust effect size add-on for Wald tests and anova

# Effect size summary that can compliment your summary output. All conversions are based on the robust effect size index estimator. Object handling is handled by the 'lmtest' package.
# #' @param object The model object.
#c#' @param vcov. The variance covariance matrix to use. Defaults to one of the heteroskedasticity consistent estimator. Note that this differs from the default in lmtest.
# #' @param ... Arguments passed to waldtest
# #' @keywords wald test, anova
# #' @return Returns an anova-like test and the robust effect size index for comparing two or more models.
# #' @importFrom sandwich vcovHC
# #' @importFrom lmtest waldtest
# #' @importFrom stats pt qnorm
# calc_resi.waldtest = function(object, ..., vcov.=sandwich::vcovHC){
#   x = waldtest(object, ...=..., vcov=vcov.)
#   cbind(x, S = RESI::chisq2S(qnorm(pt(x[,'t value'], df = object$df.residual))^2, 1, object$df.residual))
# }


#' RESI from gee or geeglm object
#' This function calculate the RESI from geeglm model object
#' @param object The model object
#' @return returns the ANOVA-type summary table with RESI estimate for each factor
#' @export
calc_resi.geeglm <- function(object, ...){
  x = as.matrix(summary(object)$coefficients)
  #sample size
  N = length(summary(object)$clusz)
  output = cbind(x, RESI = RESI::chisq2S(x[, 'Wald'], 1, N))
  return(output)
}

calc_resi.gee <- function(object, ...){
  x = as.matrix(summary(object)$coefficients)
  #sample size
  N = length(unique(object$id))
  output = cbind(x, RESI = RESI::chisq2S(x[, 'Robust z']^2, 1, N))
  return(output)
}

#' RESI from lme object
#' This function calculate the RESI from lme model object
#' @param object The lme model object
#' @return returns the ANOVA-type summary table with RESI estimate for each factor
#' @export
calc_resi.lme <- function(object, ...){
  x = as.matrix(summary(object)$tTable)
  #sample size
  N = summary(object)$dims$ngrps[1]
  # robust se
  robust.var = diag(clubSandwich::vcovCR(object, type = "CR3"))
  robust.se = sqrt(robust.var)
  output = cbind(x, 'Robust.SE' = robust.se, 'Robust Wald' = (x[, 'Value']^2/robust.var), RESI = RESI::chisq2S(x[, 'Value']^2/robust.var, 1, N))
  return(output)
}




