% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rsq2S.R
\name{Rsq2S}
\alias{Rsq2S}
\title{Covert R^2 to S}
\usage{
Rsq2S(Rsq)
}
\arguments{
\item{Rsq}{Numeric, R^2}
}
\value{
Returns an estimate of R^2 based on the RESI
}
\description{
Converts R^2, the partial coefficient of determination, to
robust effect size index (S) using the formula from Vandekar, Tao, & Blume (2020).
}
\details{
The formula for the conversion is:

\eqn{S = \sqrt((-R^2)/(R^2 - 1))}
}
\examples{
# consider a moderate effect size of R^2 = 0.1
Rsq2S(0.1)
# this corresponds to a RESI of 0.333
}
