# Internal functions

#' Transformation from Wald test statistics to squared RESI
#' @param chisq The chi-square statistic for the parameter(s) of interest.
#' @param df Number of degrees of freedom of the chi-square statistic.
#' @param rdf Model residual degrees of freedom or number of independent samples.
#' @return Returns a scalar or vector argument of the squared robust effect size index estimate.
#' @export
chisq2Ssq = function(chisq, df, rdf){
  S = (chisq - df)/rdf
}

#' Compute the robust effect size index estimator from chi-squared statistic.
#'
#' This function computes the robust effect size index from Vandekar, Rao, & Blume (2020).
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' For mixed effects models, RESI is conditional on the average correlation
#' structure within subjects.
#' @param chisq The chi-square statistic for the parameter of interest.
#' @param df Number of degrees of freedom of the chi-square statistic.
#' @param n Number of independent samples.
#' @keywords power
#' @return Returns a scalar or vector argument of the the robust effect size index estimate.
#' @export
chisq2S <- function(chisq, df, n){
  S = (chisq - df)/n
  sqrt(ifelse(S<0, 0, S))
}

#' Compute the robust effect size index estimator from t statistic
#'
#' This function computes the robust effect size index from Vandekar, Rao, & Blume (2020).
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' @param t The t statistic for the parameter of interest.
#' @param n Number of independent samples
#' @param rdf Model residual degrees of freedom or number of independent samples.
#' @keywords power
#' @return Returns a scalar or vector argument of the the robust effect size index estimate.
#' @export
t2S <- function(t, n, rdf){
  (t*sqrt(2))/(sqrt(n*rdf))*exp(lgamma(rdf/2) - lgamma((rdf-1)/2))
}

#' Compute the robust effect size index estimator from Z statistic
#'
#' This function computes the robust effect size index from Vandekar, Rao, & Blume (2020).
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' @param z The Z statistic for the parameter of interest.
#' @param n Number of independent samples.
#' @keywords power
#' @return Returns a scalar or vector argument of the the robust effect size index estimate.
#' @export
z2S <- function(z, n){
  z/sqrt(n)
}


#' Compute the robust effect size index estimator from F-statistic
#'
#' This function computes the robust effect size index from Vandekar, Rao, & Blume (2020).
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' @param f The chi-square statistic for the parameter of interest.
#' @param df Number of degrees of freedom of the chi-square statistic.
#' @param rdf Model residual degrees of freedom or number of independent samples.
#' @keywords power
#' @return Returns a scalar or vector argument of the the robust effect size index estimate.
#' @export
f2S <- function(f, df, rdf){
  S = (f*df*(rdf-2)/rdf - df)/rdf
  sqrt(ifelse(S<0, 0, S))
}

#' Convert S to Cohen's d
#'
#' Converts the robust effect size index (S) to Cohen's d using the formula from
#' Vandekar, Rao, & Blume (2020).
#' @param S Numeric, the robust effect size index.
#' @param pi Numeric, the sampling proportions.
#' @keywords power
#' @return Returns an estimate of Cohen's \emph{d} based on the RESI.
#' @details The pi parameter comes from the fact that Cohen's d doesn't account
#' for unequal sample proportions in the population, but S does.
#'
#' The default is set to a natural value 1/2, which corresponds to a case
#' control design, for example, where sampling proportions always are
#' controlled by the experimenter.
#' @export
S2d = function(S, pi=0.5){
  S * sqrt(1/pi + 1/(1-pi))
}

#' Covert Cohen's \emph{d} to S
#' Converts Cohen's \emph{d} robust effect size index (S) using the formula from
#' Vandekar, Rao, & Blume (2020).
#' @param d Numeric, value of Cohen's \emph{d}.
#' @param pi Numeric, the sampling proportions.
#' @keywords power
#' @return Returns an estimate the robust effect size index
#' @details The pi parameter comes from the fact that Cohen's d doesn't account
#' for unequal sample proportions in the population, but S does.
#'
#' The default is set to a natural value 1/2, which corresponds to a case
#' control design, for example, where sampling proportions always are
#' controlled by the experimenter.
#' @export
d2S <- function(d, pi = 0.5){
  abs(d)/sqrt(1/pi + 1/(1-pi))
}

#' Covert Cohen's \emph{f}^2 to S
#'
#' Converts Cohen's \emph{f}^2 to robust effect size index (S)
#' using the formula from Vandekar, Rao, & Blume (2020).
#' @param fsq Numeric, value of Cohen's \emph{f}^2.
#' @keywords power
#' @return Returns an estimate the robust effect size index
#' @export
fsq2S <- function(fsq){
  sqrt(fsq)
}

#' Covert S to Cohen's \emph{f}^2
#'
#' Converts robust effect size index (S) to Cohen's \emph{f}^2
#' (effect size for multiple regression) using the formula from Vandekar, Rao, & Blume (2020).
#' @param S Numeric,the robust effect size index.
#' @keywords power
#' @return Returns an estimate of Cohen's \emph{f}^2 based on the RESI
#' @export
S2fsq <- function(S){
  S^2
}

#' Covert S to R^2
#'
#' Converts robust effect size index (S) to R^2, the partial
#' coefficient of determination, using the formula from Vandekar, Rao, & Blume (2020).
#' @param S Numeric, the robust effect size index.
#' @keywords power
#' @return Returns an estimate of R^2 based on the RESI
#' @export
S2Rsq <- function(S){
  S^2/(1+S^2)
}

#' Covert R^2 to S
#'
#' Converts R^2, the partial coefficient of determination, to
#' robust effect size index (S) using the formula from Vandekar, Rao, & Blume (2020).
#' @param Rsq Numeric, R^2
#' @keywords power
#' @return Returns an estimate of R^2 based on the RESI
#' @export
Rsq2S <- function(Rsq){
  sqrt((-Rsq)/(Rsq-1))
}


#' Non-parametric bootstrap sampling
#' @export
#' @param data data frame; The data frame that need bootstrapping
#' @param id.var character; for clustered/longitudinal data, the name of id variable used as sampling unit.
#' @return Returns a data frame containing bootstrapped data.
boot.samp <- function(data, id.var = NULL) {
  params = as.list(match.call()[-1])
  if (is.matrix(data)) data = as.data.frame(data)
  if (is.null(id.var)) {
    boot.ind = sample(1:nrow(data), replace = TRUE)
    boot.data = data[boot.ind, ]
  } else {
    boot.ind = sample(unique(data[, id.var]), replace = TRUE)
    boot.data = data[unlist(lapply(boot.ind, function(x) which(x == data[, id.var]))), ]
  }
  return(boot.data)
}

#' Bayesian bootstrap sampling (Rubin, 1981)
#' @export
#' @param data data frame; The data frame that need bootstrapping
#' @return Returns a data frame containing weights generated via Bayesian bootstraps.
bayes.samp <- function(data) {
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

#' @export
print.resi <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nAnalysis of Effect sizes (ANOES) based on RESI:")
  cat("\nConfidence level = ", x$alpha)
  if (is.null(x$model.reduced)){
    cat("\nCall:  ", paste(deparse(x$model.full$call), sep = "\n", collapse = "\n"),  "\n",sep = "")
  }
  else{
    cat(ifelse(x$model.reduced$formula == as.formula(paste(format(formula(x$model.full)[[2]]), "~ 1")), "\nCall:  ", "\nFull Model:"),
        paste(deparse(x$model.full$call), sep = "\n", collapse = "\n"), "\n",sep = "")
    if (!(x$model.reduced$formula == as.formula(paste(format(formula(x$model.full)[[2]]), "~ 1")))) cat("Reduced Model:", paste(deparse(x$model.reduced$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  }

  # summary table
  if (!is.null(x$coefficients)){
    cat("\nCoefficient Table \n")
    print(round(x$coefficients, digits = digits))
  }

  # anova table
  if (!is.null(x$anova)){
    cat("\n\n")
    print(round(x$anova, digits = digits))
  }

  # overall
  if (!(is.null(x$model.reduced))){
    if ((x$model.reduced$formula == as.formula(paste(format(formula(x$model.full)[[2]]), "~ 1")))){
      cat("\nOverall RESI comparing model to intercept-only model:\n\n")
      overall = as.data.frame(x$overall)[2,]
      rownames(overall) = NULL
      print(round(overall, digits = digits))
    }
    else{
      cat("\nOverall RESI comparing full model to reduced model:\n\n")
      overall = as.data.frame(x$overall)[2,]
      rownames(overall) = NULL
      print(round(overall, digits = digits))
    }
  }

  cat("\nNotes:")
  if (x$naive.var) cat("\n1. The RESI was calculated using the naive covariance estimator.")
  else  cat("\n1. The RESI was calculated using a robust covariance estimator.")
  if (x$boot.method == "nonparam") cat("\n2. Confidence intervals (CIs) constructed using", x$nboot,"non-parametric bootstraps. \n")
  if (x$boot.method == "bayes") cat("\n2. Credible intervals constructed using", x$nboot,"Bayesian bootstraps. \n")
  # if(nzchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep = "")

  invisible(x)
}


#' @export
print.summary.resi <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nAnalysis of Effect sizes (ANOES) based on RESI:")
  cat("\nConfidence level = ", x$alpha)
  cat("\nCall:  ", paste(deparse(x$model.full$call), sep = "\n", collapse = "\n"),  "\n",sep = "")
  cat("\nCoefficient Table \n")
  print(round(x$coefficients, digits = digits))
  invisible(x)
}

#' @importFrom graphics abline axis lines
#' @export
plot.resi <- function(x, ycex.axis = 1, alpha = NULL, ...){
  if (is.null(alpha)){
    alpha = x$alpha[1]
  }
  else{
    if (length(alpha) > 1){
      warning('\nOnly first alpha will be plotted')
      alpha = alpha[1]
    }
    if (!(alpha %in% x$alpha)){
      stop('\nSpecified alpha not found in the resi object')
    }
  }
  ll = paste(alpha/2*100, "%", sep = "")
  ul = paste((1-alpha/2)*100, "%", sep = "")
  lev <- as.factor(rownames(x$coefficients))
  plot(x = x$coefficients[,"RESI"], y = 1:length(levels(lev)),
       xlim = c(min(0, min(x$coefficients[,ll])), max(x$coefficients[,ul])),
       xlab = "RESI Estimate", yaxt = "n", ylab = "",
       main = paste("RESI Estimates and ", (1-alpha)*100, "%", " CIs", sep=""),...)
  for (i in 1:nrow(x$coefficients)){
    lines(x = c(x$coefficients[i,ll], x$coefficients[i,ul]), y = c(i,i))
  }
  axis(2, 1:length(levels(lev)), levels(lev), las = 1, cex.axis = ycex.axis)
  abline(v = 0, lty = 2)
}

#' @export
plot.summary.resi <- plot.resi

#' @export
plot.anova_resi <- function(x, alpha = NULL, ycex.axis = 1,...){
  cols = grep('%', colnames(x))
  if (is.null(alpha)){
    alpha = gsub("%", "", colnames(x)[cols[1]])
    alpha = as.numeric(alpha)*2/100
  }
  else{
    if (length(alpha) > 1){
      warning('\nOnly first alpha will be plotted')
      alpha = alpha[1]
    }
  }
  ll = paste(alpha/2*100, "%", sep="")
  ul = paste((1-alpha/2)*100, "%", sep="")
  if (!(ll %in% colnames(x))){
    stop('\nSpecified alpha not found in the resi object')
  }

  lev <- as.factor(rownames(x))
  plot(x = x$RESI, y = 1:length(levels(lev)),
       xlim = c(0, max(x[,ul])), xlab = "RESI Estimate",
       yaxt = "n", ylab = "", main = paste("RESI Estimates and ", (1-alpha)*100, "%", " CIs", sep=""),...)
  for (i in 1:nrow(x)){
    lines(x = c(x[i,ll], x[i,ul]), y = c(i,i))
  }
  axis(2, 1:length(levels(lev)), levels(lev), las = 1, cex.axis = ycex.axis)
}
