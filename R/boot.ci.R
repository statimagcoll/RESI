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
# Note: for now, it only works for m1 = 1 and the target x are the first m1 columns;
# maybe we can ask the users to specify which columns are target covariates later
#' @param func the function used to compute the Chi-sq statistics. Corresponds to different versions of RESI. By default, HC3 (Long & Ervin, 2000) will be used.
#' @param multi the distribution from which the multipliers will be drawn: 'none' = the multipliers equal constant 1 (default); 'rad' = rademacher; 'normal' = Std Normal distribution
#' @param boot.type which type of bootstrap to use.
# 1: resampling covariates along with residuals;
# 2: fixing covariates and only bootstrapping residulas;
# 3: resampling covariates and residuals independently w/ replacements
# 4. no sampling, just multipliers
#' @param correct whether the residuals with bias correction will be used, by default, TRUE.
#'
# linear model w/o intercept?
# what about m = m1 = 1?



boot.ci <- function(r = 1000, x, y, m1, func = sandwich::vcovHC, multi = 'none', boot.type = 1, alpha = 0.05, correct = TRUE){
  # check input
  ## ncol(x) >= m1
  ## nrow(x) == nrow(y)

  x <- as.data.frame(x)
  # no intercept
  fit0 <- lm(as.formula(paste0('y ~ -1 + ', paste0("V", 1:ncol(x), collapse = '+'))), data = x)

  S_hat = NULL

  dat = cbind(y, x)

  # model output: covariates & (jack-knife version of HC3) estimates of residuals
  if (correct == TRUE) {
    mod.out = cbind(x, res = residuals(fit0)/(1-hatvalues(fit0)))
  } else {
    mod.out = cbind(x, res = residuals(fit0))
  }

    for (b in 1:r){

      # version 1: resample X along with residuals non-parametrically
      if (boot.type == 1) {
        # bootstrap indicator
        boot.ind <- sample(1:nrow(x), replace = TRUE)
        boot.x <- mod.out[boot.ind, m1]
        boot.res <- mod.out[boot.ind, ncol(mod.out)] # last column
      }

      # version 2: fix X and only resample residuals (w/ replacement)
      if (boot.type == 2) {
        # bootstrap indicator
        boot.ind <- sample(1:nrow(x), replace = TRUE)
        boot.x <- x
        boot.res <- mod.out[boot.ind, ncol(mod.out)] # last column
      }

      # version 3: resampling covariates and residuals independently w/ replacements
      if (boot.type == 3) {
        # bootstrap indicator
        boot.ind1 <- sample(1:nrow(x), replace = TRUE)
        boot.x <- mod.out[boot.ind1, m1]
        boot.ind2 <- sample(1:nrow(x), replace = TRUE)
        boot.res <- mod.out[boot.ind2, ncol(mod.out)] # last column
      }

      # version 4: no resampling, just multipliers
      if (boot.type == 4){
        boot.x <- x
        boot.res <- mod.out[, ncol(mod.out)]
      }

      boot.x <- data.frame(boot.x)
      colnames(boot.x) <- paste0("V", 1:ncol(boot.x))

      # Obtain multiplier
      n = nrow(x)
      if (multi == "none") w <- rep(1, n)
      if (multi == "rad") w <- sample(x = c(-1, 1), size = n, replace = TRUE)
      if (multi == "normal") w <- rnorm(n)

      # design matrix
      dsgn.mtx = as.matrix(boot.x) #cbind(1, x) %>% as.matrix

      # computing bootstrapped outcome values
      boot.y <- dsgn.mtx %*% coefficients(fit0) + boot.res * w

      # compute chi-sq statistics
      modelfull <- lm(as.formula(paste0('boot.y ~ -1 + ', paste0('V', 1:ncol(x), collapse='+') )), data = boot.x )

      if (ncol(x) > m1){ # when there is nuisance covariates
        modelreduced <- lm( as.formula(paste0('boot.y ~ -1 + ', paste0('V', (m1+1):ncol(x), collapse='+') ) ), data=boot.x )
        # chi-square test statistics
        chi.test <- lmtest::waldtest(modelreduced, modelfull, vcov = func(modelfull), test = "Chisq")
      }
      if (ncol(x) == m1) { # when there is no nuisance covariates
        # modelreduced <- lm(boot.y ~ -1, data=boot.x )
        # chi.test <- lmtest::waldtest(modelreduced, modelfull, vcov = func(modelfull), test = "Chisq")
        chi.stat <- coefficients(modelfull) %*% solve(func(modelfull)) %*% coefficients(modelfull)
      }
      # residual d.f.
      res.df = n - m1
      df = m1

      # computing the estimated RESI
      S_hat[b] = sqrt(max((chi.stat - df)/res.df, 0) )


      # ## 1) rademacher
      #
      # y_boot_rade = dsgn.mtx%*%coefficients(fit0) + boot.res * z_rade
      #
      # modelfull = lm( as.formula(paste0('y_boot_rade ~ -1 + ', paste0('V', 1:ncol(x), collapse='+') )), data = boot.x )
      #
      # if (m1 == ncol(x)) {
      #
      # }
      #
      # modelreduced = lm( as.formula(paste0('y_boot_rade ~ -1 + ', paste0('V', (m1+1):ncol(x), collapse='+') ) ), data=boot.x )
      # chi.test = lmtest::waldtest(modelreduced, modelfull, vcov = func(modelfull), test = "Chisq")
      #
      # chi.stat = chi.test[2, 'Chisq']
      # resdf = chi.test[2, 'Res.Df']
      # df = chi.test[2, 'Df']
      # # compute the estimated RESI
      # S_hat_boot_rade[b] = sqrt(max((chi.stat - df)/resdf, 0) )
      #
      #
      # ## 2) Normal dist
      # y_boot_norm = dsgn.mtx%*%coefficients(fit0) + residuals(fit0)/(1-hatvalues(fit0)) *z_norm
      # modelfull = lm( as.formula(paste0('y_boot_norm ~ -1 + ', paste0('V', 1:ncol(x), collapse='+') )), data=x )
      # modelreduced = lm( as.formula(paste0('y_boot_norm ~ -1 + ', paste0('V', (m1+1):ncol(x), collapse='+') ) ), data=x )
      # chi.test = lmtest::waldtest(modelreduced, modelfull, vcov = func(modelfull), test = "Chisq")
      # chi.stat = chi.test[2, 'Chisq']
      # resdf = chi.test[2, 'Res.Df']
      # df = chi.test[2, 'Df']
      # # compute the estimated RESI
      # S_hat_boot_norm[b] = sqrt(max((chi.stat - df)/resdf, 0) )
      #
      # ## 3) non-parametric bootstrap
      # ind = sample(1:nrow(x), size = nrow(x), replace = TRUE)
      # dat.boot = dat[ind,]
      # modelfull = lm( as.formula(paste0('y ~ -1 + ', paste0('V', 1:ncol(x), collapse='+') )), data=dat.boot )
      # modelreduced = lm( as.formula(paste0('y ~ -1 + ', paste0('V', (m1+1):ncol(x), collapse='+') ) ), data=dat.boot )
      # chi.test = lmtest::waldtest(modelreduced, modelfull, vcov = func(modelfull), test = "Chisq")
      # chi.stat = chi.test[2, 'Chisq']
      # resdf = chi.test[2, 'Res.Df']
      # df = chi.test[2, 'Df']
      # # compute the estimated RESI
      # S_hat_np.boot[b] = sqrt(max((chi.stat - df)/resdf, 0) )

    } # end of loop 'b'

    CI = quantile(S_hat, probs = c(alpha/2, 1-alpha/2))
    # colnames(CI) = c("LB", "UB")
    output = list(CI = CI, boot.type = boot.type, multiplier = multi, alpha = alpha, correct = correct)
    return(output)
}









