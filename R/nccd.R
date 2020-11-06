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


#' Function for MD Bessel interval
#'
#' Code taken by SNV from https://www1.maths.leeds.ac.uk/~john/software/ncc/ncc.r.
#'
#' @param u something.
#' @param nu something else.
#' @param lambda parameter value.
gb=function(u,nu,lambda) {
  # noncentral chi pdf wrt  Bessel base
  case1=(u>0 & lambda>0)
  case2=(u>0 & lambda==0)
  case3=(u==0)
  if(case1) pdf=lambda^(-nu)* exp(-.5*lambda^2)*besselI(lambda*u,nu)/
    sqrt(besselI(u^2,nu))
  if(case2) pdf=(u/2)^nu/(gamma(nu+1)*sqrt(besselI(u^2,nu)))
  if(case3) pdf=exp(-.5*lambda^2)/sqrt(gamma(nu+1))/(2^(nu/2))
  pdf
}


#' Function for MD radial interval
#'
#' Code taken by SNV from https://www1.maths.leeds.ac.uk/~john/software/ncc/ncc.r.
#'
#' @param u something.
#' @param nu something else.
#' @param lambda parameter value.
gr=function(u,nu,lambda) {
  # noncentral chi pdf wrt radial base
  ul=lambda*u
  if(ul>0) pdf=ul^(-nu)*besselI(ul,nu)*exp(-.5*(u^2+lambda^2))
  else pdf= (.5^nu)*exp(-.5*(u^2+lambda^2))/gamma(nu+1)
  pdf
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
  if(y<=u1) {ll=0; lu=0}
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
  f = function(b) F(y+2*b,p,y+b)-F(y,p,y+b)-(1-alpha)
  #f2 = function(b) integrate(function(x) dchisq(x, df = p, ncp = (y+b)^2), lower=y^2, upper=(y+2*b)^2)$value - (1-alpha)
  big=1; obj=-1
  while(obj<0) {big=2*big; obj=f(big)}
  bb=uniroot(f,c(0,big))$root
  lu=y+bb
  int8=c(ll,lu)
  int8
}





#' MD Bessel confidence interval for chi-square noncentrality parameter.
#'
#' Code taken by SNV from https://www1.maths.leeds.ac.uk/~john/software/ncc/ncc.r.
#'
#' @export
#' @param y numeric, value of chi random variable.
#' @param p integer, degrees of freedom.
#' @param alpha probability for confidence interval.
#' @importFrom stats constrOptim
ncc.ci.mdb=function(y,p,alpha=0.05) {
  u0=Finv(1-alpha,p,0); nu=(p-2)/2; zz=-qnorm(alpha/2)
  if(y<=u0) ll=0
  else {
    ly=lamfind(y,p,1-alpha)
    #    f=function(lambda) F(y,p,lambda)-(1-alpha)
    #    lbig=lambound(y,p,1-alpha)
    #    ly=uniroot(f,c(0,lbig))$root
    g0=gb(0,nu,ly); gy=gb(y,nu,ly)
    if(g0>=gy) ll=ly
    else {
      f=function(cd) (log(gb(y,nu,cd[1])/gb(cd[2],nu,cd[1])))^2 +
        (F(y,p,cd[1])-F(cd[2],p,cd[1]) - (1-alpha))^2 # cd1=lam; cd2=cc
      start=c(y-zz,y-2*zz)
      if(start[1]<=0) start[1]=.1*(y+.1)
      if(start[2]<=0) start[2]=.05*(y+.1)
      ll=constrOptim(start,f,grad=NULL,ui=diag(2),ci=c(0,0))$par[1]
    }
  }
  f=function(cd) (log(gb(y,nu,cd[1])/gb(cd[2],nu,cd[1])))^2 +
    (F(cd[2],p,cd[1])-F(y,p,cd[1]) - (1-alpha))^2 # cd1=lam; cd2=d
  start=c(y+zz,y+2*zz)
  lu=constrOptim(start,f,grad=NULL,ui=diag(2),ci=c(0,0))$par[1]
  int6=c(ll,lu)
  int6
}

#' MD radial confidence interval for chi-square noncentrality parameter.
#'
#' Code taken by SNV from https://www1.maths.leeds.ac.uk/~john/software/ncc/ncc.r.
#'
#' @export
#' @param y numeric, value of chi random variable.
#' @param p integer, degrees of freedom.
#' @param alpha probability for confidence interval.
#' @importFrom stats constrOptim qnorm
ncc.ci.mdr=function(y,p,alpha=0.05) {
  u0=Finv(1-alpha,p,0); nu=(p-2)/2; zz=-qnorm(alpha/2)
  if(y<=u0) ll=0
  else {
    ly=lamfind(y,p,1-alpha)
    #    f=function(lambda) F(y,p,lambda)-(1-alpha)
    #    lbig=lambound(y,p,1-alpha)
    #    ly=uniroot(f,c(0,lbig))$root
    g0=gr(0,nu,ly); gy=gr(y,nu,ly)
    if(g0>=gy) ll=ly
    else {
      f=function(cd) (log(gr(y,nu,cd[1])/gr(cd[2],nu,cd[1])))^2 +
        (F(y,p,cd[1])-F(cd[2],p,cd[1]) - (1-alpha))^2 # cd1=lam; cd2=cc
      start=c(y-zz,y-2*zz)
      if(start[1]<=0) start[1]=.1*(y+.1)
      if(start[2]<=0) start[2]=.05*(y+.1)
      ll=constrOptim(start,f,grad=NULL,ui=diag(2),ci=c(0,0))$par[1]
    }
  }
  f=function(cd) (log(gr(y,nu,cd[1])/gr(cd[2],nu,cd[1])))^2 +
    (F(cd[2],p,cd[1])-F(y,p,cd[1]) - (1-alpha))^2 # cd1=lam; cd2=d
  start=c(y+zz,y+2*zz)
  lu=constrOptim(start,f,grad=NULL,ui=diag(2),ci=c(0,0))$par[1]
  int7=c(ll,lu)
  int7
}


#' Return four confidence intervals for noncentrality parameter of a chi-square distribution.
#'
#' Code taken by SNV from https://www1.maths.leeds.ac.uk/~john/software/ncc/ncc.r.
#'
#' @export
#' @param yl numeric, value of chi random variable.
#' @param p integer, degrees of freedom.
#' @param alpha, probability for confidence interval.
ncc.ints=function(yl,p,alpha=0.05) {
  int1=ncc.ci.central(yl,p,alpha)
  int2=ncc.ci.sr(yl,p,alpha)
  #int3 = ncc.ci.mdb(yl,p,alpha)
  #int4 = ncc.ci.mdr(yl, p, alpha)
  ans=rbind(int1,int2)#,int3,int4)
  rownames(ans)=c("Central",
                  "Symmetric range")
                  #"MD Bessel",
                  #"MD radial")
  colnames(ans) = c('L', 'U')
  ans
}
