# Internal functions

#' Transformation from Wald test statistics to squared RESI
#' @param chisq The chi-square statistic for the parameter(s) of interest.
#' @param df Number of degrees of freedom of the chi-square statistic.
#' @param n Number of independent samples.
#' @return Returns a scalar or vector argument of the squared robust effect size index estimate.
#' @noRd
chisq2Ssq = function(chisq, df, n){
  S = (chisq - df)/rdf
}

#' Non-parametric bootstrap sampling
#'
#' @param data data frame; The data frame that need bootstrapping
#' @param id.var character; for clustered/longitudinal data, the name of id variable used as sampling unit.
#' @return Returns a data frame containing bootstrapped data.
#' @noRd
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
#' @param data data frame; The data frame that need bootstrapping
#' @return Returns a data frame containing weights generated via Bayesian bootstraps.
#' @noRd
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
