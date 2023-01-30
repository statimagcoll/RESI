#' Robust Effect Size Index (RESI) Point Estimation
#'
#' This function will estimate the robust effect size (RESI) from Vandekar, Rao, & Blume (2020).
#' The overall RESI is estimated via a Wald test. RESI is (optionally) estimated for each factor in coefficients-style table.
#' RESI is (optionally) estimated for each variable/interaction in an Anova-style table
#' for models with existing Anova methods. This function is the building block for the \code{\link{resi}} function.
#' @param model.full \code{lm, glm, nls, survreg, coxph, hurdle, zeroinfl, gee, geeglm} or \code{lme} model object.
#' @param model.reduced Fitted model object of same type as model.full. By default `NULL`; the same model as the full model but only having intercept.
#' @param data Data.frame or object coercible to data.frame of model.full data (required for some model types).
#' @param vcovfunc The variance estimator function for constructing the Wald test statistic. By default, sandwich::vcovHC (the robust (sandwich) variance estimator).
#' @param coefficients Logical, whether to produce a coefficients (summary) table with the RESI columns added. By default = `TRUE`.
#' @param anova Logical, whether to produce an Anova table with the RESI columns added. By default = `TRUE`.
#' @param Anova.args List, additional arguments to be passed to Anova function.
#' @param vcov.args List, additional arguments to be passed to vcovfunc.
#' @param unbiased Logical, whether to use the unbiased or alternative T/Z statistic to RESI conversion. By default, `TRUE`. See details.
#' @param ... Ignored.
#' @importFrom aod wald.test
#' @importFrom car Anova
#' @importFrom lmtest waldtest
#' @importFrom regtools nlshc
#' @importFrom sandwich vcovHC
#' @importFrom stats coef df.residual formula glm hatvalues pf predict quantile residuals update vcov
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
#' \code{\link{t2S}}, \code{\link{z2S}}, \code{\link{t2S_alt}}, and
#' \code{\link{z2S_alt}} for more details on the formulas.
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

resi_pe <- function(model.full, ...){
  UseMethod("resi_pe")
}

#' @describeIn resi_pe RESI point estimation
#' @export
resi_pe.default <- function(model.full, model.reduced = NULL, data,
                            coefficients = TRUE, vcovfunc = sandwich::vcovHC, Anova.args = list(),
                            vcov.args = list(), unbiased = TRUE, ...){
  if (missing(data)){
    data = model.full$model
  }
  else{
    data = as.data.frame(data)
  }

  if (is.null(model.reduced)){
    form.reduced = as.formula(paste(format(formula(model.full)[[2]]), "~ 1"))
    if (is.null(model.full$na.action)){
      model.reduced <- update(model.full, formula = form.reduced, data = data)
    }
    else{
      model.reduced <- update(model.full, formula = form.reduced, data = model.full$model)
    }
  }

  # dealing with additional vcov args
  if (length(vcov.args) == 0){
    vcovfunc2 <- vcovfunc
  } else{
    vcov.args = c(formals(vcovfunc)[1], vcov.args)
    vcovfunc2 <- function(x) {
      vcov.args[[1]] = x
      do.call(vcovfunc, vcov.args)}
  }

  # overall
  wald.test = lmtest::waldtest(model.reduced, model.full, vcov = vcovfunc2, test = 'Chisq')
  stats = wald.test$Chisq[2]
  overall.df = wald.test$Df[2]
  res.df = wald.test$Res.Df[2]
  overall.resi.hat = chisq2S(stats, overall.df, nrow(data))
  wald.test[2, 'RESI'] = overall.resi.hat
  output <- list(model.full = list(call = model.full$call, formula = formula(model.full)),
                 model.reduced = list(call = model.reduced$call, formula = formula(model.reduced)),
                 estimates = overall.resi.hat, overall = wald.test)
  names.est = "Overall"
  names(output$estimates) = names.est

  # coefficients table (z statistics)
  if (coefficients){
    coefficients.tab <- lmtest::coeftest(model.full, vcov. = vcovfunc2)
    coefficients.df = data.frame(coefficients.tab[,'Estimate'], coefficients.tab[,'Std. Error'],
                                 coefficients.tab[,'z value'], coefficients.tab[,'Pr(>|z|)'], row.names = rownames(coefficients.tab))
    colnames(coefficients.df) = colnames(coefficients.tab)
    if (unbiased){
      coefficients.df[,'RESI'] = z2S(coefficients.df[,'z value'], nrow(data))
    }
    else{
      coefficients.df[,'RESI'] = suppressWarnings(z2S_alt(coefficients.df[,'z value'], nrow(data)))
    }
    output$coefficients = coefficients.df
    output$estimates = c(output$estimates, coefficients.df$RESI)
    names.est = c(names.est, rownames(coefficients.df))
    names(output$estimates) = names.est
  }

  return(output)
}

#' @describeIn resi_pe RESI point estimation for generalized linear models
#' @export
resi_pe.glm <- function(model.full, model.reduced = NULL, data, anova = TRUE,
                        coefficients = TRUE, vcovfunc = sandwich::vcovHC,
                        Anova.args = list(), vcov.args = list(), unbiased = TRUE, ...){
  output <- resi_pe.default(model.full = model.full, model.reduced = model.reduced,
                            data = data, coefficients = coefficients, vcovfunc = vcovfunc,
                            Anova.args = Anova.args, vcov.args = vcov.args, unbiased = unbiased, ...)

  if (length(vcov.args) == 0){
    vcovfunc2 <- vcovfunc
  } else{
    vcov.args = c(formals(vcovfunc)[1], vcov.args)
    vcovfunc2 <- function(x) {
      vcov.args[[1]] = x
      do.call(vcovfunc, vcov.args)}
  }

  # Anova table (Chi sq statistics)
  if (anova){
    suppressMessages(anova.tab <- do.call(car::Anova, c(list(mod = model.full, test.statistic = 'Wald', vcov. = vcovfunc2), Anova.args)))
    anova.tab[,'RESI'] = chisq2S(anova.tab[,'Chisq'], anova.tab[,'Df'], nrow(data))
    output$anova = anova.tab[which(rownames(anova.tab) != "Residuals"),]
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
  return(output)
}

#' @describeIn resi_pe RESI point estimation for linear models
#' @export
resi_pe.lm <- function(model.full, model.reduced = NULL, data, anova = TRUE,
                       coefficients = TRUE, vcovfunc = sandwich::vcovHC,
                       Anova.args = list(), vcov.args = list(), unbiased = TRUE, ...){
  if (missing(data)){
    data = model.full$model
  }
  else{
    data = as.data.frame(data)
  }

  if (is.null(model.reduced)){
    form.reduced = as.formula(paste(format(formula(model.full)[[2]]), "~ 1"))
    if (is.null(model.full$na.action)){
      model.reduced <- update(model.full, formula = form.reduced, data = data)
    }
    else{
      model.reduced <- update(model.full, formula = form.reduced, data = model.full$model)
    }
  }

  if (length(vcov.args) == 0){
    vcovfunc2 <- vcovfunc
  } else{
    vcov.args = c(formals(vcovfunc)[1], vcov.args)
    vcovfunc2 <- function(x) {
      vcov.args[[1]] = x
      do.call(vcovfunc, vcov.args)}
  }

  # overall
  wald.test = lmtest::waldtest(model.reduced, model.full, vcov = vcovfunc2, test = 'Chisq')
  stats = wald.test$Chisq[2]
  overall.df = wald.test$Df[2]
  res.df = wald.test$Res.Df[2]
  overall.resi.hat = chisq2S(stats, overall.df, nrow(data))
  wald.test[2, 'RESI'] = overall.resi.hat
  output <- list(model.full = list(call = model.full$call, formula = formula(model.full)),
                 model.reduced = list(call = model.reduced$call, formula = formula(model.reduced)),
                 estimates = overall.resi.hat, overall = wald.test)
  names.est = "Overall"
  names(output$estimates) = names.est

  # coefficients table (t statistics)
  if (coefficients){
    coefficients.tab <- lmtest::coeftest(model.full, vcov. = vcovfunc2)
    coefficients.df = data.frame(coefficients.tab[,'Estimate'], coefficients.tab[,'Std. Error'],
                                 coefficients.tab[,'t value'], coefficients.tab[,'Pr(>|t|)'], row.names = rownames(coefficients.tab))
    colnames(coefficients.df) = colnames(coefficients.tab)
    if (unbiased){
      coefficients.df[,'RESI'] = t2S(coefficients.df[,'t value'], nrow(data), model.full$df.residual)
    } else{
      coefficients.df[,'RESI'] = t2S_alt(coefficients.df[,'t value'], model.full$df.residual)
    }

    output$coefficients = coefficients.df
    output$estimates = c(output$estimates, coefficients.df$RESI)
    names.est = c(names.est, rownames(coefficients.df))
    names(output$estimates) = names.est
  }

  # Anova table (F statistics)
  if (anova){
    suppressMessages(anova.tab <- do.call(car::Anova, c(list(mod = model.full,
                                                             vcov. = vcovfunc2), Anova.args)))
    anova.tab[,'RESI'] = f2S(anova.tab[,'F'], anova.tab[,'Df'], res.df)
    output$anova = anova.tab[which(rownames(anova.tab) != "Residuals"),]
    output$estimates = c(output$estimates, anova.tab$RESI)
    names.est = c(names.est, rownames(anova.tab))
    names(output$estimates) = names.est
    output$estimates = output$estimates[which(output$estimates != "Residuals")]
    class(output$anova) = c("anova_resi", class(output$anova))
  }

  if(identical(vcovfunc, stats::vcov)){
    output$naive.var = TRUE
  }
  else{
    output$naive.var = FALSE
  }
  return(output)
}



#' @describeIn resi_pe RESI point estimation for hurdle models
#' @export
#' @describeIn resi_pe RESI point estimation for geeglm object
#' @export
resi_pe.geeglm <- function(object, anova = TRUE, ...){

  data = object$data
  # sample size
  N = length(summary(object)$clusz)
  # total num of observations
  tot_obs = nrow(data)
  # working corstr of the object
  corstr_spec = object$corstr

  # longitudinal RESI
  if (anova) {
    mod_tab = anova(object)
    mod_tab$RESI = RESI::chisq2S(mod_tab[,'X2'], mod_tab$Df, N)
  } else {
    mod_tab = as.data.frame(summary(object)$coefficients)
    mod_tab = cbind(mod_tab, RESI = RESI::chisq2S(mod_tab[,'Wald'], 1, N))
  }

  # model form
  form = formula(object)
  # independence model
  mod_ind = update(object, id = c(1:nrow(data)), corstr = corstr_spec, data = data)

  # the var-cov matrix estimate from the independence model
  # Note: this is the estimate for Cov[(\hat{\beta}_ind - \beta_0)]
  cov_ind = vcov(mod_ind)
  # Convert it to the estimate for \Sigma_ind = Cov(\sqrt{N}(\hat{\beta}_ind - \beta_0))
  cov_ind = cov_ind * N
  # # Convert it to the estimate for \Sigma_cs = \Sigma_ind * m
  # cov_ind = cov_ind * tot_obs/N

  # The var-cov matrix estimate from the analysis model (the one considering correlation )
  # Note: this is the estimate for Cov(\hat{\beta}_long)
  cov_long = vcov(object)
  # convert it to the estimate for \Sigma_long = Cov(\sqrt{N}(\hat{\beta}_long - \beta_0))
  cov_long = cov_long * N

  # The pm-RESI estimates
  # mod_tab_ind = anova(mod_ind)
  coef_mod = coef(mod_ind)
  n.vars = length(coef_mod ) # num of coefficients in the model
  var_names = names(coef_mod) # coefficient names ...
  term_names = labels(terms(mod_ind)) # terms ...
  n.terms = length(term_names) # num of terms ...
  # Cmat: the matrix identifying the "relative" coefficients
  # it may be problematic if a variable's name is contained in another one's...
  if (anova){
    Cmat = matrix(NA, ncol = n.vars, nrow = n.terms)
    colnames(Cmat) = var_names
    rownames(Cmat) = term_names
    for (term in 1:n.terms){
      Cmat[term, ] = grepl(term_names[term], var_names, fixed = TRUE)
    }
  } else {
    Cmat = diag(n.vars) == 1
    colnames(Cmat) = rownames(Cmat) = var_names
  }

  # the pm-RESI for each variable..
  Wald_ind = NULL
  for (term in rownames(Cmat)){
    index = Cmat[term, ]
    sub_coef = coef(object)[index] # the coefficient(s) from the input object
    sub_vcov_cor = cov_long[index, index] # the vcov from the gee model
    sub_vcov_ind = cov_ind[index, index] # the vcov from the independence model
    # # ESS
    # w = sub_coef %*% solve(sub_vcov_cor) %*% sub_coef / sub_coef %*% solve(sub_vcov_ind) %*% sub_coef
    # ess = w * tot_obs
    # mod_tab[term, "ess"] = ess
    # mod_tab[term, "pm-RESI"] = sqrt(N / ess * mod_tab[term, "RESI"]^2)

    # Wald test statistics using vcov from the independence model
    Wald_ind = N * sub_coef %*% solve(sub_vcov_ind) %*%  sub_coef
    df = sum(index) # df
    S_ind =  chisq2S(Wald_ind, df = df, n = N)
    mod_tab[term, "pm-RESI"] = sqrt(N/tot_obs * S_ind^2)
  }


  output <- list(model.full = list(call = object$call, formula = formula(object)),
                 resi = mod_tab)

  output$robust.var = "The GEE models always use robust (sandwich) variance estiamtor."
  output$info = "Function modified Jan 29 10:32pm"
  return(output)
}

#' @describeIn resi_pe RESI point estimation for gee object
#' @export
resi_pe.gee <- function(model.full, ...){
  x = as.matrix(summary(model.full)$coefficients)
  #sample size
  N = length(unique(model.full$id))
  output <- list(model.full = list(call = model.full$call, formula = formula(model.full)),
                 coefficients =  as.data.frame(cbind(x, RESI = RESI::chisq2S(x[, 'Robust z']^2, 1, N))))
  output$naive.var = FALSE
  return(output)
}

#' @describeIn resi_pe RESI point estimation for lme object
#' @importFrom clubSandwich vcovCR
#' @export
resi_pe.lme <- function(model.full, vcovfunc = clubSandwich::vcovCR, vcov.args = list(), ...){
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
  output = list(model.full = list(call = model.full$call, formula = formula(model.full)),
                coefficients =  as.data.frame(cbind(x, 'Robust.SE' = robust.se,
                                                    'Robust Wald' = (x[, 'Value']^2/robust.var), RESI = RESI::chisq2S(x[, 'Value']^2/robust.var, 1, N))))
  if(identical(vcovfunc, stats::vcov)){
    output$naive.var = TRUE
  }
  else{
    output$naive.var = FALSE
  }
  return(output)
}

#' @describeIn resi_pe RESI point estimation for lmerMod object
#' @export
resi_pe.lmerMod <- function(model.full, vcovfunc = clubSandwich::vcovCR, vcov.args = list(), ...){
  x = as.matrix(summary(model.full)$coefficients)
  #sample size
  N = summary(model.full)$ngrps

  # robust se
  if (identical(vcovfunc, stats::vcov)) {
    output = list(coefficients = cbind(x, 'Wald' = x[, 't value']^2, RESI = RESI::chisq2S(x[, 't value']^2, 1, N)), naive.var = TRUE)
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
    output = list(coefficients = cbind(x, 'Robust.SE' = robust_se, 'Robust Wald' = (x[, 'Estimate']^2/robust_var), RESI = RESI::chisq2S(x[, 'Estimate']^2/robust_var, 1, N)), naive.var = FALSE)
  }

  return(output)
}


