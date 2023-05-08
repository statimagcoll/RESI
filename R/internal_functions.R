# Internal functions

#' Transformation from Wald test statistics to squared RESI
#' @param chisq The chi-square statistic for the parameter(s) of interest.
#' @param df Number of degrees of freedom of the chi-square statistic.
#' @param n Number of independent samples.
#' @return Returns a scalar or vector argument of the squared robust effect size index estimate.
#' @noRd
chisq2Ssq = function(chisq, df, n){
  S = (chisq - df)/n
}

#' Statistic for bootstrapping RESI estimates with boot package
#'
#' @param dat data.frame; original dataset
#' @param inds bootstrapped indices
#' @param mod.full full model passed to resi
#' @param mod.reduced reduced model passed to resi
#' @param boot.method non-parametric (default) or Bayesian bootstrap
#' @return Returns a vector of RESI estimates for a single bootstrap replicate
#' @noRd
resi_stat = function(dat, inds, mod.full, mod.reduced, boot.method = "nonparam",
                     cluster = FALSE, clvar = NULL, mod.dat = NULL, ...){
  # nonparametric bootstrap
  if (boot.method == "nonparam"){
    if(!cluster){
      boot.data = dat[inds,]
      mod.full = try(update(mod.full, data = boot.data), silent = T)
      if (!(is.null(mod.reduced))){
        mod.reduced = try(update(mod.reduced, data = boot.data), silent = T)}
    }
     else{
      # for clustered data, dat is unique ids and inds is which of the ids
      boot.data = mod.dat[unlist(lapply(inds, function(x) which(x == mod.dat[, clvar]))), ]
      boot.data$bootid = rep(1:length(unique(mod.dat[, clvar])),
                unlist(lapply(inds, function(x) length(which(x==mod.dat[,clvar])))))

      mod.full = try(update(mod.full, data = boot.data, id = bootid), silent = T)
      if (!(is.null(mod.reduced))){
        mod.reduced = try(update(mod.reduced, data = boot.data, id = bootid), silent = T)}
    }}

  else{
    # Bayesian bootstrap
    n = nrow(dat)
    repeat{
      # Generate the random numbers from unif(0, 1)
      u = runif(n-1)
      u.sort = sort(u)
      g = c(u.sort, 1) - c(0, u.sort)
      if (sum(g == 0) == 0) break
    } # end `repeat`
    boot.data = cbind(dat, g)
    mod.full = try(update(mod.full, data = boot.data, weights = g), silent = T)
    if (!(is.null(mod.reduced))){
      mod.reduced = try(update(mod.reduced, data = boot.data, weights = g), silent = T)}
  }

  if(!(inherits(mod.full, "try-error") | inherits(mod.reduced, "try-error"))){
    out = resi_pe(mod.full, model.reduced = mod.reduced, data = boot.data, ...)$estimates
  } else{
    out = NA
  }

  # for models with fail checking, count nboot-nrow(t0) for number of failures
  return(out)
}


#' Non-parametric bootstrap sampling
#'
#' @param data data frame; The data frame that need bootstrapping
#' @param id.var character; for clustered/longitudinal data, the name of id variable used as sampling unit.
#' @return Returns a data frame containing bootstrapped data.
#' @noRd
boot.samp = function(data, id.var = NULL) {
  params = as.list(match.call()[-1])
  if (is.matrix(data)) data = as.data.frame(data)
  if (is.null(id.var)) {
    boot.ind = sample(1:nrow(data), replace = TRUE)
    boot.data = data[boot.ind, ]
  } else {
    boot.ind = sample(unique(data[, id.var]), replace = TRUE)
    boot.data = data[unlist(lapply(boot.ind, function(x) which(x == data[, id.var]))), ]
    boot.data$bootid = rep(1:length(unique(data[, id.var])),
                           unlist(lapply(boot.ind, function(x) length(which(x==data[,id.var])))))
  }
  return(boot.data)
}

#' Bayesian bootstrap sampling (Rubin, 1981)
#' @param data data frame; The data frame that need bootstrapping
#' @return Returns a data frame containing weights generated via Bayesian bootstraps.
#' @noRd
bayes.samp = function(data) {
  if (is.matrix(data)) data = as.data.frame(data)
  n = nrow(data)
  repeat{
    # Generate the random numbers from unif(0, 1)
    u = runif(n-1)
    u.sort = sort(u)
    g = c(u.sort, 1) - c(0, u.sort)
    if (sum(g == 0) == 0) break
  } # end `repeat`
  boot.data = cbind(data, g)
  return(boot.data)
}

#' Copied regtools::nlshc (to be removed later, current reverse dependency issue)
#' @return Returns robust covariance for nls model
#' @noRd
vcovnls = function (nlsout, type = "HC") {
  b = coef(nlsout)
  m = nlsout$m
  resid = m$resid()
  hmat = m$gradient()
  xhm = hmat
  yresidhm = resid + hmat %*% b
  lmout = stats::lm(yresidhm ~ xhm - 1)
  sandwich::vcovHC(lmout, type)
}
