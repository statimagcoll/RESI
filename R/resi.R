#' Robust effect size add-on for lm objects
#'
#' Effect size summary that can compliment your summary output. All conversions are based on the robust effect size index estimator.
#' @param object The model object.
#' @param vcov. The variance covariance matrix to use. Defaults to one of the heteroskedasticity consistent estimator. Note that this differs from the default in lmtest.
#' @param ... Arguments passed to coeftest
#' @keywords lm effect size
#' @return Returns a summary of tests and the robust effect size index for linear models.
#' @details Stuff and stuff.
#' @importFrom sandwich vcovHC
#' @importFrom lmtest coeftest
#' @importFrom stats pt qnorm
#' @export
resi.lm = function(object, vcov.=sandwich::vcovHC, ...){
  x = coeftest(object, vcov.=vcov., ...)
  cbind(x, S = RESI::chisq2S(qnorm(pt(x[,'t value'], df = object$df.residual))^2, 1, object$df.residual))
}

#' Robust effect size add-on for glm
#'
#' Effect size summary that can compliment your summary output. All conversions are based on the robust effect size index estimator.
#' @param object The model object.
#' @param vcov. The variance covariance matrix to use. Defaults to one of the heteroskedasticity consistent estimator. Note that this differs from the default in lmtest.
#' @param ... Arguments passed to coeftest
#' @keywords effect size glm
#' @return Returns a summary of tests and the robust effect size index for generalized linear models.
#' @details Stuff and more stuff.
#' @importFrom sandwich vcovHC
#' @importFrom lmtest coeftest
#' @export
resi.glm = function(object, vcov.=sandwich::vcovHC, ...){
  x = coeftest(object, vcov.=vcov., ...)
  cbind(x, S = RESI::chisq2S(x[,'z value']^2, 1, object$df.residual))
}

#' Robust effect size add-on for Wald tests and anova
#'
#' Effect size summary that can compliment your summary output. All conversions are based on the robust effect size index estimator. Object handling is handled by the 'lmtest' package.
#' @param object The model object.
#' @param vcov. The variance covariance matrix to use. Defaults to one of the heteroskedasticity consistent estimator. Note that this differs from the default in lmtest.
#' @param ... Arguments passed to waldtest
#' @keywords wald test, anova
#' @return Returns an anova-like test and the robust effect size index for comparing two or more models.
#' @importFrom sandwich vcovHC
#' @importFrom lmtest waldtest
#' @importFrom stats pt qnorm
#' @export
resi.waldtest = function(object, ..., vcov.=sandwich::vcovHC){
  x = waldtest(object, ...=..., vcov=vcov.)
  cbind(x, S = RESI::chisq2S(qnorm(pt(x[,'t value'], df = object$df.residual))^2, 1, object$df.residual))
}
