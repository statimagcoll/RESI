
## ------------------------------------------
# Functions computing the CIs for non-centrality parameters of F distributions
# Modified based on Kent's code for CIs for noncentral chi-square distributions
## To compute the Confidence intervals for the non-centrality parameters from F distributions
## Original Code taken by SNV from https://www1.maths.leeds.ac.uk/~john/software/ncc/ncc.r.
## Taken from  Kent, J.T. and Hainsworth, T.J. (1995) Confidence intervals for the
## noncentral chi-squared distribution. J Stat Planning Inf 46, 147-159.

## ------------------------------------------

# Consider y^2 ~ \chi^2_p(\lambda^2)

# Given \lambda, we can find a (1-\alpha) prob interval for y, or
# given y, we can find a (1-\alpha) confidence interval for \lambda.

# Four methods for each are carried out by the function below.
# The results are summarized in an 8 x 3 matrix, where the first column
# gives the lower endpoint, the second column gives the upper endpoint,
# and the third column is a label for the type of interval.


#' Finding the bound of non-centrality parameter
#' @param y numeric, value of F random variable;
#' @param df1 integer, degrees of freedom 1.
#' @param df2 integer, degrees of freedom 2.
#' @param alpha significence level

lambound_f=function(y, df1, df2, alpha) {
  lambda=1; obj=1
  while(obj>0) {
    lambda=2*lambda
    obj=pf(y, df1, df2, lambda) - alpha
  }
  lambda
}


#' Finds root of distribution function minus alpha.
#'
#' Original code taken by SNV from https://www1.maths.leeds.ac.uk/~john/software/ncc/ncc.r.
#'
#' @param y numeric, value of F random variable.
#' @param df1 integer, degrees of freedom 1.
#' @param df2 integer, degrees of freedom 2.
#' @param alpha probability.
#' @importFrom stats uniroot

lamfind_f=function(y,df1, df2, alpha) {
  if(alpha<1e-4 | pf(y, df1, df2, 0) < alpha) stop("bad lamfind")
  lbig=lambound_f(y, df1, df2, alpha)
  f=function(lambda) pf(y, df1, df2, lambda) - alpha
  lambda=uniroot(f,c(0,lbig))$root
}




#' Confidence interval for F noncentraliy parameter
#' Written by Kaidi based on Kent's code
#' Original kent's Code taken by SNV from https://www1.maths.leeds.ac.uk/~john/software/ncc/ncc.r.
#'
#' Modified based on functions ncc.ci.sr() and ncc.ci.central()
#'
#' @export
#' @param y numeric, value of F random variable.
#' @param df1 integer, degrees of freedom 1.
#' @param df2 integer, degrees of freedom 2.
#' @param alpha two-tailed probability for confidence interval.
#' @return xxx
#' @export
ncf.ci=function(y, df1, df2, alpha=0.05) {

  p = pf(y, df1, df2, 0)

  if (p < 1-alpha/2) { # p<1-alpha/2 means ncp = 0 still fail to satisfy pf(y, ..., ncp = 0) >= 1 - alpha/2
    ll = 0
    if(p > alpha) { # when ll = 0, put all the probability on the upper bound side
      ## Then if pf(y, ..., ncp = 0) > alpha, upper bound > 0
      lu = lamfind_f(y, df1, df2, alpha)
    } else{ # ie, p<=alpha --> even putting all probability on the upper side,
            #still fail to get a value for ncp that is > 0
      lu = 0
    }
  } else{ # compute the central version of CI for the noncetrality parameter
    ll = lamfind_f(y, df1, df2, 1-alpha/2)
    lu = lamfind_f(y, df1, df2, alpha/2)
  }

  int5=c(ll,lu)
  int5
}




