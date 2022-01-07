
#' Analysis of Effect Sizes (ANOES) based on the Robust Effect Size index (RESI) for (generalized) linear regression models
#' @param model.full model formula for full model. `glm` object
#' @param model.reduced the reduced model or comparison model. By dafult, `NULL`
#' @param roubst.var logic, indicating whether robust (sandwich) variance estimator should be used in the construction of test statistics. Default to `TRUE` and `sandwich::vcovHC` will be used to estimate HC3.
#' @param boot.type the type of bootstraps used to construct the CIs for the RESI estimates. By default, along with `multi = "none"`, non-parametric bootstraps will be used.
#' @param multi default to `"none"`. It indicates the distribution from with the multiplers for the wild bootstrapp will be drawn. `"rad"` = Rademacher distribution and `"normal"` = Std Normal distribution.
#' @param alpha The significance level of the constructed CI. By default, 0.05.
#' @param nboot the number of bootstraps that will be implemented to construct the CIs. By default, 1000 bootstraps will be applied.
#' @importFrom stats coefficients hatvalues pf quantile residuals update
#' @return
#' @export


anoes <- function(model.full, model.reduced = NULL,
                     robust.var = TRUE,
                     boot.type = 1, multi = 'none',
                     nboot = 1000, alpha = 0.05){

  output = boot.ci(model.full = model.full, model.reduced = model.reduced,
                    robust.var = robust.var,
                    boot.type = boot.type, multi = multi,
                    r = nboot,
                    alpha = alpha)$ANOES

  return(output)
}



