#' Calculating Effective Sample Size
#'
#' This function calcualte the effective sample size (ESS) from Kang et. al (2022)

ess <- function(x, ...) {
UseMethod("ess")
}

#' This function calculates the effective sample size from a fitted linear mixed effect model object (lme) or a GEE model object (gee)
#' @param obj `lme` or `gee` object
#' @param constr matrix. The linear constraint used to form an hypothesis regading which the RESI is built. For example, for parameters `beta`, `constr = diag(c(1, 1, 1))` means testing `constr %*% beta = 0`
#' By default, `constr = NULL` where the ESS for each of the parameters will be computed.
#'
#' @export
ess.lme <- function(obj, constr = NULL){
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
    chisq = beta^2 / diag(cov.beta)
    # observed RESI from original model
    RESI.sq.ori = chisq2Ssq(chisq, 1, N)
  } else {
    C = diag(rep(1, length(beta)))
    chisq = t(C %*% beta) %*% solve(C %*% cov.beta %*% t(C)) %*% (C %*% beta)
    # observed RESI from original model
    RESI.sq.ori = chisq2Ssq(chisq, sum(C), N)
  }

  # 2. RESI (from the independence model)
  mod.ind = glm(form, data = data)
  beta.ind = coef(mod.ind)
  cov.beta.ind = sandwich::vcovHC(mod.ind)
  if (total){
    chisq = beta.ind^2 / diag(cov.beta.ind)
    # observed RESI from original model
    RESI.sq.ind = chisq2Ssq(chisq, 1, N)
  } else {
    C = diag(rep(1, length(beta)))
    chisq = t(C %*% beta) %*% solve(C %*% cov.beta %*% t(C)) %*% (C %*% beta)
    # observed RESI from original model
    RESI.sq.ind = chisq2Ssq(chisq, sum(C), N)
  }

  # 3. EFFECTIVE SAMPLE SIZE (ESS)
  ## total num of observations
  tot_obs = nrow(data)
  w = RESI.sq.ori/RESI.sq.ind
  ess = tot_obs * w

  return(ess)
}

#' Calculating the effective sample size from a fitted GEE model.
#' @param obj a fitted `geeglm` object
#' @param constr matrix, linear constraint
#' @export

ess.geeglm <- function(obj, total = TRUE){
  # dataset
  data = obj$data
  # sample size
  N = length(unique(obj$id))
  # model form
  form = formula(obj)

  # 1. RESI (from original model)
  # WALD STATISTICS
  ## point estimates
  beta = as.matrix(obj$coefficients)
  ## sandwich covariance
  cov.beta = vcov(obj)
  ## wald statistics from the full model
  if (total){
    chisq = beta^2 / diag(cov.beta)
    # observed RESI from original model
    RESI.sq.ori = chisq2Ssq(chisq, 1, N)
  } else {
    C = diag(rep(1, length(beta)))
    chisq = t(C %*% beta) %*% solve(C %*% cov.beta %*% t(C)) %*% (C %*% beta)
    # observed RESI from original model
    RESI.sq.ori = chisq2Ssq(chisq, sum(C), N)
  }

  # 2. RESI (from independence model)
  mod.ind = glm(form, data = data)
    # geeglm(form, data = data, id = obj$id, constr = "independence")
  beta.ind = coef(mod.ind)
  cov.beta.ind = sandwich::vcovHC(mod.ind)
  if (total){
    chisq.ind = beta.ind^2 / diag(cov.beta.ind)
    # observed RESI from independence model
    RESI.sq.ind = chisq2Ssq(chisq.ind, 1, N)
  } else {
    C = diag(rep(1, length(beta.ind)))
    chisq.ind = t(C %*% beta.ind) %*% solve(C %*% cov.beta.ind %*% t(C)) %*% (C %*% beta.ind)
    # observed RESI from independence model
    RESI.sq.ind = chisq2Ssq(chisq.ind, sum(C), N)
    }

  # 3. EFFECTIVE SAMPLE SIZE (ESS)
  ## total num of observations
  tot_obs = nrow(data)
  w = RESI.sq.ori/RESI.sq.ind
  ess = tot_obs * w
  return(ess)
}


