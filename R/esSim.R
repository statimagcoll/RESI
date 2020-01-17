#' Run one simulation analysis to assess bias and variance of effect size index.

#' Includes code to assess coverage and width of confidence intervals from
#' Kent, J. T., & Hainsworth, T. J. (1995). Confidence intervals for the noncentral
#' chi-squared distribution. Journal of Statistical Planning and Inference,
#' 46(2), 147–159. https://doi.org/10.1016/0378-3758(94)00104-4.
#' See image.statMap for additional arguments
#'
#' @export
#' @param simsetup list, output from simSetup function.
#' @param shape numeric, shape parameter for gamma distribution.
#' @param n integer, sample size for this simulation.
#' @param alpha numeric, probability for confidence interval.
#' @importFrom stats as.formula lm rgamma rnorm
#' @importFrom lmtest waldtest
#' @importFrom sandwich vcovHC
simFunc = function(simsetup, shape=0.5, n=500, alpha=0.05){
  S = simsetup[['S']]
  m1 = simsetup[['m1']]
  m = m1 + simsetup[['m0']]
  x = matrix(rnorm(n * m), nrow=n, ncol=m ) %*% simsetup[['Vsqrt']]
  beta = simsetup[['beta']]
  hetero = simsetup[['hetero']]
  if(hetero){
    # here, the average variance is 1, but the variance depends on the value of x
    y = x %*% beta + rgamma(n, shape = shape, rate = sqrt(shape/ x[,1]^2) ) - sqrt(shape * x[,1]^2) # mean of gamma will be sqrt(shape * x)
  } else {y = x %*% beta + rgamma(n, shape = shape, rate = sqrt(shape) ) - sqrt(shape)}
  x = as.data.frame(x)
  modelfull = lm( as.formula(paste0('y ~ -1 + ', paste0('V', 1:m, collapse='+') )), data=x )
  modelreduced = lm( as.formula(paste0('y ~ -1 + ', paste0('V', (m1+1):m, collapse='+') )), data=x )

  chistat = lmtest::waldtest(modelreduced, modelfull, vcov=sandwich::vcovHC(modelfull), test = 'Chisq')
  chi = sqrt(chistat[2,'Chisq'])
  resdf = chistat[2,'Res.Df']
  df = chistat[2,'Df']
  Shat = sqrt(max((chistat[2,'Chisq'] - df)/resdf, 0) )
  bias = (Shat - S) # bias
  CI = ncc.ints(chi, df, alpha=alpha)/sqrt(resdf)
  widths = diff(t(CI))
  c(bias=bias,
    coverage.central=(S>=CI[1,1] & S<CI[1,2]),
    coverage.sr=(S>=CI[2,1] & S<CI[2,2]),
    width.central=widths[1],
    width.sr=widths[2] ) # was 4, now 3
}

#' Generate design matrix for regression of Y on X with dependence among covariates.
#'
#' @export
#' @param S numeric, effect size index. See https://arxiv.org/abs/1902.07232.
#' @param m1 integer, number of target covariates.
#' @param m0 integer, number of nuisance covariates.
#' @param rhosq numeric in [0,1] controlling dependence between target and nuisance covariates.
#' @param hetero logical, should Y be simulated with heteroskedasticity?
#' @importFrom pracma sqrtm
simSetup = function(S=0.25, m1=3, m0=2, rhosq=0.6, hetero=TRUE){
  m = m0 + m1
  # First set of variables are covariates of interest. Next set are nuisance. No intercept
  Vsqrt = rbind(cbind(diag(m1), matrix(sqrt(rhosq/m1/m0), nrow=m1, ncol=m0)), cbind(matrix(sqrt(rhosq/m1/m0), nrow=m0, ncol=m1), diag(m0)))
  Vsqrt = pracma::sqrtm(Vsqrt)$B
  if(hetero){
    mat = diag(m1) - rhosq/m1 * matrix(1, nrow=m1, ncol=m1)
    # For checking my math
    #mat2 = Vsqrt %*% matrix(rnorm(m * 10000), nrow=m)
    # X0 = mat2[(m1+1):m,]
    # X1 = mat2[1:m1,]
    # A and B are classical notation from Boos and Stefanski of components of the sandwich estimator
    A = diag(m1) - matrix(1, nrow=m1, ncol=m1) * rhosq/m1
    X1TDeltaX0X0TX1 = matrix(c(3 * rhosq/m1, rep(rhosq/m1, m1-1)), nrow=m1, ncol=m1)
    B <- if(m1==1) matrix(3) else diag(c(3, rep(1, m1-1)) )
    B = B + rhosq/m0/m1 * (m0 + 2 * m0 * rhosq/m1) * matrix(1, nrow=m1, ncol=m1) - (X1TDeltaX0X0TX1  + t(X1TDeltaX0X0TX1 ) )
    Binv = solve(B)
    beta = rep(c(S/sqrt(sum(A %*% Binv %*% A)), 0), c(m1, m0) )
  } else {
    beta = rep(c(S/sqrt(m1 * (1-rhosq)), 0), c(m1, m0)  ) # assumes constant covariance between X1 and X0 variables and independence otherwise cov(X_1) = I, cov(X_0) = I, cov(X_1, X_0) = 1 1^T rho/sqrt(m0 * m1)
  }
  # return required objects and simulation parameters
  list(Vsqrt=Vsqrt, beta=beta, S=S, m1=m1, m0=m0, rhosq=rhosq, hetero=hetero)
}
