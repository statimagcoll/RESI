# Function to estimate the Robust effect size index (Simon Vandekar et.al, 2020) for linear regression models
# different version of estimated RESI will be outputted,
# along with the specified version of CIs

#' Robust Effect Size index (RESI) estimation for linear regression models
#' @param model model formula for
#' @param model0 the comparison model
#' @param data the data frame containing the variables speficiedd in the model formula(s)
#' @param .vcov the covariance estimator for coefficient tested, by default, HC3 will be used (according to the help document of sandwich::vcovHC)
#' @param method The method constructing CI. "wild.norml" = the wild bootstrapping with standard normal multiplier;
#' "wild.rademacher" = wild bootstrapping with multiplier following rademacher distribution; "f" = F distribution; "chisq" = Chi-square distribution
#' @param alpha The type I error rate for the constructed CI
#' @param nboot if a bootstrapping method were specified to construct the CI, nboot = the number of bootstraps that will be implemented.
#' @return sssss
#' @export
#'
# Note:
# 1) I don't know how to call the fucntion in another R file.
# 2) The design matrix used in bootstraps should accomodate to the model w/ or w/o intercept

anoes.lm <- function(model, model0 = NULL, data,
                     .vcov = sandwich::vcovHC,
                     method = c("f", "chisq", "wild.normal", "wild.rademacher", "nonparametric.boots"),
                     nboot = 1000, alpha = 0.05){

  # Computing estimated RESI

  ## variable names from the model formula
  var.names = all.vars(model)

  if (is.null(model0)) model0= as.formula(paste0(var.names[1], " ~ 1"))

  model.full = lm(model, data = data)
  model.reduced = lm(model0, data = data)

  wald.output = lmtest::waldtest(model.reduced, model.full, vcov=.vcov(model.full), test = 'Chisq')
  # Chi-sq statistic
  chistat = wald.output$Chisq[2]
  # residual d.f.
  res.df = wald.output$Res.Df[2]
  # model d.f.
  df = wald.output$Df[2]

  # estimated RESI
  # S.hat = sqrt(max((chistat - df)/res.df, 0) )
  S.hat = chisq2S(chisq = chistat, df = df, rdf = res.df)

  # CIs
  # via Chi-square distribution
  CI_chi = ncc.ints(sqrt(chistat), df, alpha=alpha)[1, ]/sqrt(res.df)

  # via F-distribution
  CI_f = sqrt(ncf.ci(y = chistat/df, df1 = df, df2 = res.df, alpha = alpha)/res.df)

  # via wild bootstraps
  S_hat_boot_rade = S_hat_boot_norm = S_hat_np.boot = NULL

  for (b in 1:nboot){
    # multiplier
    z_rade = sample(c(1, -1), size = nrow(data), replace = TRUE) # rademacher
    z_norm = rnorm(nrow(data)) # standard normal

    # design matrix
    dsgn.mtx = as.matrix(cbind(1, data[, var.names[-1]]) )

    ## Note: Jack-knife version of HC3 estimator was used to estimate the residuals

    # compute chi-sq statistics
    ## 1)rademacher
    y_boot_rade = dsgn.mtx %*% coefficients(model.full) + residuals(model.full)/(1-hatvalues(model.full)) *z_rade
    data0 = cbind(data, y_boot_rade)
    model.full.boot = lm(update(model, y_boot_rade ~ .), data = data0)
    model.reduced.boot = lm(update(model0, y_boot_rade ~ .), data = data0)
    wald.output.boot = lmtest::waldtest(model.reduced.boot, model.full.boot, vcov=.vcov(model.full.boot), test = 'Chisq')
    # Chi-sq statistic
    chi.stat.boot = wald.output.boot$Chisq[2]
    # compute the estimated RESI in that bootstrap
    S_hat_boot_rade[b] = sqrt(max((chi.stat.boot - df)/res.df, 0) )

    ## 2) Normal dist
    y_boot_norm = dsgn.mtx %*% coefficients(model.full) + residuals(model.full)/(1-hatvalues(model.full)) *z_norm
    data0 = cbind(data, y_boot_norm)
    model.full.boot = lm(update(model, y_boot_norm ~ .), data = data0)
    model.reduced.boot = lm(update(model0, y_boot_norm ~ .), data = data0)
    wald.output.boot = lmtest::waldtest(model.reduced.boot, model.full.boot, vcov=.vcov(model.full.boot), test = 'Chisq')
    # Chi-sq statistic
    chi.stat.boot = wald.output.boot$Chisq[2]
    # compute the estimated RESI in that bootstrap
    S_hat_boot_norm[b] = sqrt(max((chi.stat.boot - df)/res.df, 0) )

    ## 3) non-parametric bootstrap
    ind = sample(1:nrow(data), size = nrow(data), replace = TRUE)
    dat.boot = data[ind,]
    model.full.boot = lm(model, data=dat.boot )
    model.reduced.boot = lm(model0, data=dat.boot )
    wald.output.boot = lmtest::waldtest(model.reduced.boot, model.full.boot, vcov = .vcov(model.full.boot), test = "Chisq")
    chi.stat.boot = wald.output.boot$Chisq[2]
    # compute the estimated RESI
    S_hat_np.boot[b] = sqrt(max((chi.stat.boot - df)/res.df, 0) )

  } # end of loop 'b'

  CI_wild.rade = quantile(S_hat_boot_rade, probs = c(alpha/2, 1-alpha/2))
  CI_wild.norm = quantile(S_hat_boot_norm, probs = c(alpha/2, 1-alpha/2))
  CI_np = quantile(S_hat_np.boot, probs = c(alpha/2, 1-alpha/2))

  CIs = rbind(CI_f,
              CI_chi,
              CI_wild.rade,
              CI_wild.norm,
              CI_np
              )
  rownames(CIs) = c("F", "Chi-sq", "Wild.Rademacher", "Wild.Normal", "Non-parametric")
  # colnames(CIs) = c("LB", "UB")
  output = list(est = S_hat,
                CI = CIs)
  return(output)
}



## test
y = sample(20:100, 100, replace = TRUE)
x1 = rnorm(100, mean = 20)
x2 = runif(100)*5 - 20
data = data.frame(y, x1, x2)

anoes.lm(model = y ~ x1+x2, data = data)

fit  = lm(y ~ -1 + x1 + x2)

