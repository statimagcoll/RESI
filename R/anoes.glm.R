# Function to conduct Analysis of Effect Sizes based on the RESI (Simon Vandekar et.al, 2020)
# for (generalized) linear regression models
# based on glm()

#' Analysis of Effect Sizes based on the Robust Effect Size index (RESI) for (generalized) linear regression models
#' @param model.full model formula for full model
#' @param model.reduced the reduced model or comparison model. By dafult, NULL
#' @param data the data frame containing the variables speficiedd in the model formula(s)
#' @param .vcov the covariance estimator for tested coefficient(s), by default, sandwich::vcovHC will be used.
#' @param method The method constructing CI.
#' "wild.norml" = the wild bootstrapping with standard normal multiplier;
#' "wild.rademacher" = wild bootstrapping with multiplier following rademacher distribution;
#' "f" = non-central F distribution;
#' "chisq" = non-cental Chi-square distribution
#' By default, the non-parametric bootstrap will be used.
#' @param alpha The significance level of the constructed CI. By default, 5%
#' @param nboot if bootstrapped CI was specified in `method`, nboot = the number of bootstraps that will be implemented. By default, 1000 bootstraps will be applied.
#' @importFrom stats coefficients hatvalues pf quantile residuals update
#' @return
#' @export


anoes.glm <- function(model.full, model.reduced = NULL,
                     .vcov = sandwich::vcovHC,
                     method = "bootstrap",
                     nboot = 1000, alpha = 0.05){

  # data
  data = model.full$model
  # `model family` or `the linkage function` in glm()
  family = model.full$family$family

  # 1. Obtain the point estimate of the RESI
  ## variable names from the model formula
  var.names = all.vars(model.full$formula)
  ## build the null model
  if (is.null(model.reduced)) model.reduced = as.formula(paste0(var.names[1], " ~ 1"))
  model.reduced = glm(model.reduced, data = data, family = family)

  ## compute the Wald-type stat for the covariates of interest
  wald.test = lmtest::waldtest(model.reduced, model.full, vcov=.vcov(model.full), test = 'Chisq')
  ## Chi-sq statistic
  chi2stat = wald.test$Chisq[2]
  ## residual d.f.
  res.df = wald.test$Res.Df[2]
  ## model d.f.
  df = wald.test$Df[2]

  ## point estimate of the RESI
  S.hat = chisq2S(chisq = chistat, df = df, rdf = res.df)

  # CIs
  if (tolower(method) %in% c('wild.normal', 'wild.rad', 'f', 'chisq', 'bootstrap')){
    # non-parametric bootstrap
    if (tolower(method) == "bootstrap"){

    }
    # via Chi-square distribution
    if (tolower(method) == "chisq"){
      S.CI = ncc.ints(sqrt(chistat), df, alpha=alpha)[1, ]/sqrt(res.df)
    }
    # via F-distribution
    if (tolower(method) == "f"){
      S.CI = sqrt(ncf.ci(y = chistat/df, df1 = df, df2 = res.df, alpha = alpha)/res.df)
    }
    # via wild bootstraps


  } else{stop("Please correctly specify the method for CI:
              (1) 'bootstrap': non-parametric bootstrap (default);
              (2) 'Chisq': (non-central) Chi-squared CI;
              (3) 'F': (non-central) F CI;
              (4) 'wild.normal': wild bootstrap with std normal mutipliers;
              (5)  wild.rad': wild bootstrap with rademacher multipliers; ")
              }




  rownames(CIs) = c("F", "Chi-sq", "Wild.Rademacher", "Wild.Normal", "Non-parametric")
  # colnames(CIs) = c("LB", "UB")
  output = list(est = S.hat,
                CI = CIs)
  return(output)
}



