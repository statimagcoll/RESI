#######
# Function to construct the CIs for
# different versions of RESI via Wild Bootstrap and non-parametric bootstrap.
######

#' function generating CIS via bootstrap
#' @export
#' @param r numeric, the number of boostrap replicates
#' @param x data frame containing covariates
#' @param y outcome variable
#' @param m1 the number of target parameter(s)
#' @param func the function used to compute the Chi-sq statistics. corresponds to different versions of RESI

## Questions:
# 1. funciton residual() is not jack-knife version, but hatvalues() is, problematic??
#



boot.ci <- function(r = 1000, x, y, m1, func, alpha = 0.05){

  fit0 <- lm(as.formula(paste0('y ~ -1 + ', paste0("V", 1:ncol(x), collapse = '+'))), data = x)
  S_hat_boot_rade = S_hat_boot_norm = S_hat_np.boot = NULL
  dat = cbind(y, x)

    for (b in 1:r){

      # multiplier
      z_rade = sample(c(1, -1), size = nrow(x), replace = TRUE) # rademacher
      z_norm = rnorm(nrow(x)) # standard normal

      # design matrix
      dsgn.mtx = x %>% as.matrix #cbind(1, x) %>% as.matrix

      ## Note: Jack-knife version of HC3 estimator was used to estimate the residuals

      # compute chi-sq statistics
      ## 1)rademacher
      y_boot_rade = dsgn.mtx%*%coefficients(fit0) + residuals(fit0)/(1-hatvalues(fit0)) *z_rade
      modelfull = lm( as.formula(paste0('y_boot_rade ~ -1 + ', paste0('V', 1:ncol(x), collapse='+') )), data=x )
      modelreduced = lm( as.formula(paste0('y_boot_rade ~ -1 + ', paste0('V', (m1+1):ncol(x), collapse='+') ) ), data=x )
      chi.test = lmtest::waldtest(modelreduced, modelfull, vcov = func(modelfull), test = "Chisq")
      chi.stat = chi.test[2, 'Chisq']
      resdf = chi.test[2, 'Res.Df']
      df = chi.test[2, 'Df']
      # compute the estimated RESI
      S_hat_boot_rade[b] = sqrt(max((chi.stat - df)/resdf, 0) )

      ## 2) Normal dist
      y_boot_norm = dsgn.mtx%*%coefficients(fit0) + residuals(fit0)/(1-hatvalues(fit0)) *z_norm
      modelfull = lm( as.formula(paste0('y_boot_norm ~ -1 + ', paste0('V', 1:ncol(x), collapse='+') )), data=x )
      modelreduced = lm( as.formula(paste0('y_boot_norm ~ -1 + ', paste0('V', (m1+1):ncol(x), collapse='+') ) ), data=x )
      chi.test = lmtest::waldtest(modelreduced, modelfull, vcov = func(modelfull), test = "Chisq")
      chi.stat = chi.test[2, 'Chisq']
      resdf = chi.test[2, 'Res.Df']
      df = chi.test[2, 'Df']
      # compute the estimated RESI
      S_hat_boot_norm[b] = sqrt(max((chi.stat - df)/resdf, 0) )

      ## 3) non-parametric bootstrap
      ind = sample(1:nrow(x), size = nrow(x), replace = TRUE)
      dat.boot = dat[ind,]
      modelfull = lm( as.formula(paste0('y ~ -1 + ', paste0('V', 1:ncol(x), collapse='+') )), data=dat.boot )
      modelreduced = lm( as.formula(paste0('y ~ -1 + ', paste0('V', (m1+1):ncol(x), collapse='+') ) ), data=dat.boot )
      chi.test = lmtest::waldtest(modelreduced, modelfull, vcov = func(modelfull), test = "Chisq")
      chi.stat = chi.test[2, 'Chisq']
      resdf = chi.test[2, 'Res.Df']
      df = chi.test[2, 'Df']
      # compute the estimated RESI
      S_hat_np.boot[b] = sqrt(max((chi.stat - df)/resdf, 0) )
    } # end of loop 'b'

    CIs = rbind(quantile(S_hat_boot_rade, probs = c(alpha/2, 1-alpha/2)),
                quantile(S_hat_boot_norm, probs = c(alpha/2, 1-alpha/2)),
                quantile(S_hat_np.boot, probs = c(alpha/2, 1-alpha/2)))
    rownames(CIs) = c("Rademacher", "Normal", "Non-parametric")
    # colnames(CIs) = c("LB", "UB")
    CIs
}









