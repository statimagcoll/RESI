#' Calculating Effective Sample Size
#'
#' This function calculates the effective sample size (ESS) from Kang et. al (2022)
#' @param obj `lme` or `gee` object
#' @param total whether to calculate RESI for each covariate
#' @param robust.var whether robust covariance was used
#' @param xTab xTab
#' @param ... ignored
#' @export
ess <- function(obj, ...) {
  UseMethod("ess")
}

#' @describeIn ess Calculates the effective sample size from a fitted linear mixed effect model object (lme) or a GEE model object (gee)
#' @export
ess.lme <- function(obj, total = TRUE, ...){
  # dataset
  data = obj$data
  # sample size
  N = length(unique(getGroups(obj)))
  # model form
  form = formula(obj)

  # 1. RESI (from the original model)
  # WALD STATISTICS
  ## point estimates
  beta = as.matrix(fixed.effects(obj))
  ## sandwich covariance
  cov.beta = clubSandwich::vcovCR(obj, type = "CR3")
  ## wald statistics from the full model
  if (total){
    chisq.ori = beta^2 / diag(cov.beta)
    # observed RESI from original model
    # RESI.sq.ori = chisq2Ssq(chisq, 1, N)
  } else {
    C = diag(rep(1, length(beta)))
    chisq.ori = t(C %*% beta) %*% solve(C %*% cov.beta %*% t(C)) %*% (C %*% beta)
    # observed RESI from original model
    # RESI.sq.ori = chisq2Ssq(chisq, sum(C), N)
  }

  # 2. RESI (from the independence model)
  mod.ind = glm(form, data = data)
  beta.ind = coef(mod.ind)
  cov.beta.ind = sandwich::vcovHC(mod.ind)
  if (total){
    chisq.ind = beta.ind^2 / diag(cov.beta.ind)
    # observed RESI from original model
    # RESI.sq.ind = chisq2Ssq(chisq, 1, N)
  } else {
    C = diag(rep(1, length(beta)))
    chisq.ind = t(C %*% beta) %*% solve(C %*% cov.beta %*% t(C)) %*% (C %*% beta)
    # observed RESI from original model
    # RESI.sq.ind = chisq2Ssq(chisq, sum(C), N)
  }

  # 3. EFFECTIVE SAMPLE SIZE (ESS)
  ## total num of observations
  tot_obs = nrow(data)
  w = chisq.ori/chisq.ind
  ess = tot_obs * w
  ess = cbind(tot_obs, chisq.ori, chisq.ind, w, ess)
  colnames(ess) = c("Total obs", "Original Wald", "Ind Wald", "weight", "ESS")
  return(ess)
}

#' @describeIn ess Calculates the effective sample size from a fitted LME model.
#' @importFrom buildmer remove.terms
#' @importFrom lme4 findbars
#' @export
ess.lmerMod <- function(obj, xTab, robust.var = TRUE, ...){
  # dataset
  data = obj@frame
  # id variable
  id.var = names(obj@flist)
  # sample size
  N = length(unique(data[, id.var]))
  # total num of observations
  tot_obs = nrow(data)
  # model form
  form =  obj@call$formula
  ## random effect term
  re = paste0("(", lme4::findbars(form), ")")
  ## remove random effects
  # not sure if this is working as intended: form_lm = buildmer::remove.terms(form, remove = re)
  # using structure from remove.terms example
  step1 <- buildmer::remove.terms(form, re)
  re2 = paste0("(", lme4::findbars(step1), ")")
  form_lm <- buildmer::remove.terms(step1, re2)

  # The Covariance Matrix from the independence model
  mod_ind = glm(form_lm, data = data)

  if (robust.var) {
    cov_ind = sandwich::vcovHC(mod_ind)
  } else {
    cov_ind = vcov(mod_ind)
  }

  # the (wald) statistics under independence
  stat_ind = xTab[, "Estimate"]^2/diag(cov_ind) #using the same estimates
  # the RESI estimates under the
  resi_ind = RESI::chisq2Ssq(stat_ind, 1, N)
  # the pm-RESI
  resi_pm = sqrt(N/tot_obs * resi_ind^2)
  # the ESS
  ess = xTab[, "RESI"]^2/resi_ind^2 * tot_obs
  # weights
  w = xTab[, "RESI"]^2/resi_ind^2
  output = cbind(xTab, 'pm-RESI' = resi_pm, 'ESS' = ess, "weights" = w)
  return(output)
}


#' @describeIn ess Calculates the effective sample size from a fitted GEE model.
#' @export
ess.geeglm <- function(obj, total = TRUE, ...){
  # dataset
  data = obj$data
  # sample size
  N = length(summary(obj)$clusz)
  # total num of observations
  tot_obs = nrow(data)
  # model form
  form =  formula(obj)

  # The Covariance Matrix from the independence model
  mod_ind = glm(form, data = data)

  if (robust.var) {
    cov_ind = sandwich::vcovHC(mod_ind)
  } else {
    cov_ind = vcov(mod_ind)
  }

  # the (wald) statistics under independence
  stat_ind = xTab[, "Estimate"]^2/diag(cov_ind) #using the same estimates
  # the RESI estimates under the
  resi_ind = RESI::chisq2Ssq(stat_ind, 1, N)
  # the pm-RESI
  resi_pm = sqrt(N/tot_obs * resi_ind^2)
  # the ESS
  ess = xTab[, "RESI"]^2/resi_ind^2 * tot_obs
  # weights
  w = xTab[, "RESI"]^2/resi_ind^2
  output = cbind(xTab, 'pm-RESI' = resi_pm, 'ESS' = ess, "weights" = w)
  return(output)
  # # dataset
  # data = obj$data
  # # sample size
  # N = length(unique(obj$id))
  # # model form
  # form = formula(obj)
  #
  # # 1. RESI (from original model)
  # # WALD STATISTICS
  # ## point estimates
  # beta = as.matrix(obj$coefficients)
  # ## sandwich covariance
  # cov.beta = vcov(obj)
  # ## wald statistics from the full model
  # if (total){
  #   chisq.ori = beta^2 / diag(cov.beta)
  #   # observed RESI from original model
  #   # RESI.sq.ori = chisq2Ssq(chisq, 1, N)
  # } else {
  #   C = diag(rep(1, length(beta)))
  #   chisq.ori = t(C %*% beta) %*% solve(C %*% cov.beta %*% t(C)) %*% (C %*% beta)
  #   # observed RESI from original model
  #   # RESI.sq.ori = chisq2Ssq(chisq, sum(C), N)
  # }
  #
  # # 2. RESI (from independence model)
  # mod.ind = glm(form, data = data)
  #   # geeglm(form, data = data, id = obj$id, constr = "independence")
  # beta.ind = coef(mod.ind)
  # cov.beta.ind = sandwich::vcovHC(mod.ind)
  # if (total){
  #   chisq.ind = beta.ind^2 / diag(cov.beta.ind)
  #   # observed RESI from independence model
  #   # RESI.sq.ind = chisq2Ssq(chisq.ind, 1, N)
  # } else {
  #   C = diag(rep(1, length(beta.ind)))
  #   chisq.ind = t(C %*% beta.ind) %*% solve(C %*% cov.beta.ind %*% t(C)) %*% (C %*% beta.ind)
  #   # observed RESI from independence model
  #   # RESI.sq.ind = chisq2Ssq(chisq.ind, sum(C), N)
  #   }
  #
  # # 3. EFFECTIVE SAMPLE SIZE (ESS)
  # ## total num of observations
  # tot_obs = nobs(obj)
  # w = chisq.ori/chisq.ind
  # ess = tot_obs * w
  # ess = cbind(tot_obs, chisq.ori, chisq.ind, w, ess)
  # colnames(ess) = c("Total obs", "Original Wald", "Ind Wald", "weight", "ESS")
  # return(ess)
}
