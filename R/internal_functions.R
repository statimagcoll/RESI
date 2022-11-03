# Internal functions

#' Transferation from Wald test statistics to squared RESI
#' @export
#' @param chisq The chi-square statistic for the parameter(s) of interest.
#' @param df Number of degrees of freedom of the chi-square statistic.
#' @param rdf Model residual degrees of freedom or number of independent samples.
#' @return Returns a scalar or vector argument of the squared robust effect size index estimate.
chisq2Ssq = function(chisq, df, rdf){
  S = (chisq - df)/rdf
  ifelse(S<0, 0, S)
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
  cat("\n Analysis of Effect sizes (ANOES) based on RESI: ")
  cat("\n Significance level = ", x$alpha)
  cat( ifelse(is.null(x$input.object.reduced), "\n Call:  ", "\n Full Model:"),
      paste(deparse(x$input.object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (!is.null(x$input.object.reduced)) cat("\n Reduced Model:", paste(deparse(x$input.object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if(length(coef(x$input.object))) {
    print.default(format(round(x$resi, digits = digits)),
                  print.gap = 2, quote = FALSE)
  } else cat("No coefficients\n\n")

  cat("\n Notes:")
  if (x$robust.var) cat("\n 1. The RESI is calculated using the robust (sandwich) covariance estimator.")
  else cat("\n 1. The RESI is calculated using the naive covariance estimator.")
  if (x$boot.method == "nonparam") cat("\n 2. Confidence intervals (CIs) constructed using", x$nboot,"non-parametric bootstraps \n")
  if (x$boot.method == "bayes") cat("\n 2. Credible intervals constructed using", x$nboot,"Bayesian bootstraps \n")
  if(nzchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep = "")

  invisible(x)
}


#' @export
AR1_str = function(rho.e, num_visit){
  time = 0:(num_visit - 1)
  mtx <- matrix(NA, length(time), length(time))
  for (j in 1:length(time)){
    for (k in 1:length(time)){
      mtx[j, k] <- rho.e^(abs(time[j] - time[k]))
    }
  }
  return(mtx)
}

#' @export
comp_str = function(diag = 1, off_diag = 0.5, num_visit){
  mtx = matrix(NA, ncol = num_visit, nrow = num_visit)
  diag(mtx) = diag
  mtx[upper.tri(mtx)|lower.tri(mtx)] = off_diag
  return(mtx)
}

#' @export
AR1_to_comp = function(AR1){
  m = nrow(AR1)
  FUN = function(rho) sum(AR1[upper.tri(AR1)]) - m*(m-1)/2 * rho
  rho = uniroot(FUN, c(0, 10))$root
  comp = comp_str(1, rho, num_visit = m)
  return(comp)
}

#' @export
comp_to_AR1 = function(comp){
  m = nrow(comp)
  comp_rho = unique(comp[lower.tri(comp)])
  eq = paste0((m-1):1, "*rho", "^", 1:(m-1))
  FUN = function(rho) m*(m-1)/2 * comp_rho -  eval(parse(text = paste0(eq, collapse = " + ")))
  ar1_rho = uniroot(FUN, c(0, 10))$root
  ar1 = AR1_str(rho.e = ar1_rho, num_visit = m)
  return(ar1)
}
