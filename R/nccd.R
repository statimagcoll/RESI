# routine to compute 4 prob intervals and 4 conf intervals for the
# non-central chi distribution.
# Code taken by SNV from https://www1.maths.leeds.ac.uk/~john/software/ncc/ncc.r.
# Taken from  Kent, J.T. and Hainsworth, T.J. (1995) Confidence intervals for the
# noncentral chi-squared distribution. J Stat Planning Inf 46, 147-159.

# Consider y^2 ~ \chi^2_p(\lambda^2)

# Given \lambda, we can find a (1-\alpha) prob interval for y, or
# given y, we can find a (1-\alpha) confidence interval for \lambda.

# Four methods for each are carried out by the function below.
# The results are summarized in an 8 x 3 matrix, where the first column
# gives the lower endpoint, the second column gives the upper endpoint,
# and the third column is a label for the type of interval.

#' Chi distribution function.
#'
#' Code taken by SNV from https://www1.maths.leeds.ac.uk/~john/software/ncc/ncc.r.
#'
#' @export
#' @param u numeric, value of chi random variable.
#' @param p integer, degrees of freedom.
#' @param lambda, noncentrality parameter of chi distribution.
#' @importFrom stats pchisq
F=function(u,p,lambda) pchisq(u^2,p,lambda^2) # cdf of ncc

#' chi quantile function.
#'
#' Code taken by SNV from https://www1.maths.leeds.ac.uk/~john/software/ncc/ncc.r.
#'
#' @export
#' @param prob numeric in [0,1], quantile of chi random variable.
#' @param p integer, degrees of freedom.
#' @param lambda, noncentrality parameter of chi distribution.
#' @importFrom stats qchisq
Finv=function(prob,p,lambda) sqrt(qchisq(prob,p,lambda^2)) # quantile
lambound=function(y,p,alpha) {
  lambda=1; obj=1
  while(obj>0) {
    lambda=2*lambda
    obj=F(y,p,lambda)-alpha
  }
  lambda
}

#' Finds root of distribution function minus alpha.
#'
#' Code taken by SNV from https://www1.maths.leeds.ac.uk/~john/software/ncc/ncc.r.
#'
#' @param y numeric, value of chi random variable.
#' @param p integer, degrees of freedom.
#' @param alpha probability.
#' @importFrom stats uniroot
lamfind=function(y,p,alpha) {
  if(alpha<1e-4 | F(y,p,0)<alpha) stop("bad lamfind")
  lbig=lambound(y,p,alpha)
  f=function(lambda) F(y,p,lambda) - alpha
  lambda=uniroot(f,c(0,lbig))$root
}


#' Central confidence interval for chi-square noncentrality parameter.
#'
#' Code taken by SNV from https://www1.maths.leeds.ac.uk/~john/software/ncc/ncc.r.
#'
#' @export
#' @param y numeric, value of chi random variable.
#' @param p integer, degrees of freedom.
#' @param alpha probability for confidence interval.
ncc.ci.central=function(y,p,alpha=0.05) {
  u1=Finv(alpha/2,p,0); u2=Finv(1-alpha/2,p,0)
  if(y<=u1) {ll=NaN; lu=NaN}
  else {
    lu=lamfind(y,p,alpha/2)
    #    f=function(lambda) F(y,p,lambda)-alpha/2
    #    lbig=lambound(y,p,alpha/2)
    #    lu=uniroot(f,c(0,lbig))$root
    if(y <=u2) ll=0
    else ll=lamfind(y,p,1-alpha/2)
    #    f=function(lambda) F(y,p,lambda)-(1-alpha/2)
    #    lbig=lambound(y,p,1-alpha/2)
    #    ll=uniroot(f,c(0,lbig))$root
  }
  int5=c(ll,lu)
  int5
}


#' Symmetric range confidence interval for chi-square noncentrality parameter.
#'
#' Code taken by SNV from https://www1.maths.leeds.ac.uk/~john/software/ncc/ncc.r.
#'
#' @export
#' @param y numeric, value of chi random variable.
#' @param p integer, degrees of freedom.
#' @param alpha probability for confidence interval.
ncc.ci.sr=function(y,p,alpha=0.05) {
  u0=Finv(1-alpha,p,0)
  if(y<=u0) ll=0
  else {
    f=function(b) F(y,p,y-b)-F(max(y-2*b,0),p,y-b)-(1-alpha)
    bb=uniroot(f,c(0,y))$root
    ll=y-bb
  }
  y = 10
  f = function(b) F(y+2*b,p,y+b)-F(y,p,y+b)-(1-alpha)
  #f2 = function(b) integrate(function(x) dchisq(x, df = p, ncp = (y+b)^2), lower=y^2, upper=(y+2*b)^2)$value - (1-alpha)
  big=1; obj=-1
  while(obj<0) {big=2*big; obj=f(big)}
  bb=uniroot(f,c(0,big))$root
  lu=y+bb
  int8=c(ll,lu)
  int8
}

#' Return two confidence intervals for noncentrality parameter of a chi-square distribution.
#'
#' Code taken by SNV from https://www1.maths.leeds.ac.uk/~john/software/ncc/ncc.r.
#'
#' @export
#' @param yl numeric, value of chi random variable.
#' @param p integer, degrees of freedom.
#' @param alpha, probability for confidence interval.
ncc.ints=function(yl,p,alpha=0.05) {
  int5=ncc.ci.central(yl,p,alpha)
  int8=ncc.ci.sr(yl,p,alpha)
  ans=rbind(int5,int8)
  rownames(ans)=c("Central conf int",
                  "Symmetric range conf int")
  ans
}
