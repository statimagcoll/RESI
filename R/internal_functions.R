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
                     cluster = FALSE, clvar = NULL, mod.dat = NULL, overall = TRUE,
                     nest = NULL, ...){
  # nonparametric bootstrap
  if (boot.method == "nonparam"){
    if(!cluster){
      boot.data = as.data.frame(dat[inds,])
      colnames(boot.data) = colnames(dat)
      mod.full = try(update(mod.full, data = boot.data), silent = T)
      if (!(is.null(mod.reduced)) & overall){
        mod.reduced = try(update(mod.reduced, data = boot.data), silent = T)}
    }
     else{
      # for clustered data, dat is unique ids and inds is which of the ids
      boot.data = mod.dat[unlist(lapply(dat[inds], function(x) which(x == mod.dat[, clvar]))), ]
      bootid = rep(1:length(unique(mod.dat[, clvar])),
                    unlist(lapply(dat[inds], function(x) length(which(x==mod.dat[,clvar])))))
      boot.data$bootid = bootid

      mod.full = try(update(mod.full, data = boot.data, id = bootid), silent = T)
      if (!(is.null(mod.reduced)) & overall){
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
    if (!(is.null(mod.reduced)) & overall){
      mod.reduced = try(update(mod.reduced, data = boot.data, weights = g), silent = T)}
  }

  if(!(inherits(mod.full, "try-error") | (inherits(mod.reduced, "try-error") & overall))){
    out = try(resi_pe(mod.full, model.reduced = mod.reduced, data = boot.data, overall = overall, ...)$estimates, silent = T)
    if(inherits(out, "try-error")){
      out = NA
    } else{
      # for factored variables, it's possible the model will not contain all factor levels and not represent the original model
      if(length(out) != nest){
        out = NA
      }
    }
  } else{
    out = NA
  }

  # for models with fail checking, count nboot-nrow(t0) for number of failures
  return(out)
}

# function taken from regtools::nlshc (removing dependency that was causing error)
#' Statistic for bootstrapping RESI estimates with boot package
#'
#' @param nlsout nls model
#' @param type for sandwich variance
#' @return Returns robust covariance matrix for nls model
#' @noRd
r_nlshc <- function (nlsout, type = "HC")
{
  b <- coef(nlsout)
  m <- nlsout$m
  resid <- m$resid()
  hmat <- m$gradient()
  xhm <- hmat
  yresidhm <- resid + hmat %*% b
  lmout <- lm(yresidhm ~ xhm - 1)
  sandwich::vcovHC(lmout, type)
}
