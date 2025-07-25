% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resi_pe.R
\name{resi_pe}
\alias{resi_pe}
\alias{resi_pe.default}
\alias{resi_pe.glm}
\alias{resi_pe.lm}
\alias{resi_pe.nls}
\alias{resi_pe.survreg}
\alias{resi_pe.coxph}
\alias{resi_pe.hurdle}
\alias{resi_pe.zeroinfl}
\alias{resi_pe.geeglm}
\alias{resi_pe.glmgee}
\alias{resi_pe.gee}
\alias{resi_pe.lme}
\alias{resi_pe.lmerMod}
\alias{resi_pe.glmmTMB}
\alias{resi_pe.emmGrid}
\title{Robust Effect Size Index (RESI) Point Estimation}
\usage{
resi_pe(...)

\method{resi_pe}{default}(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  overall = TRUE,
  vcovfunc = sandwich::vcovHC,
  Anova.args = list(),
  vcov.args = list(),
  unbiased = TRUE,
  waldtype = 0,
  ...
)

\method{resi_pe}{glm}(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  overall = TRUE,
  vcovfunc = sandwich::vcovHC,
  Anova.args = list(),
  vcov.args = list(),
  unbiased = TRUE,
  waldtype = 0,
  ...
)

\method{resi_pe}{lm}(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  vcovfunc = sandwich::vcovHC,
  Anova.args = list(),
  vcov.args = list(),
  unbiased = TRUE,
  overall = TRUE,
  ...
)

\method{resi_pe}{nls}(
  model.full,
  model.reduced = NULL,
  data,
  coefficients = TRUE,
  anova = FALSE,
  vcovfunc = r_nlshc,
  vcov.args = list(),
  unbiased = TRUE,
  overall = TRUE,
  ...
)

\method{resi_pe}{survreg}(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  vcovfunc = vcov,
  Anova.args = list(),
  unbiased = TRUE,
  overall = TRUE,
  ...
)

\method{resi_pe}{coxph}(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  vcovfunc = vcov,
  Anova.args = list(),
  unbiased = TRUE,
  overall = TRUE,
  ...
)

\method{resi_pe}{hurdle}(
  model.full,
  model.reduced = NULL,
  data,
  coefficients = TRUE,
  anova = TRUE,
  vcovfunc = sandwich::sandwich,
  vcov.args = list(),
  unbiased = TRUE,
  overall = TRUE,
  ...
)

\method{resi_pe}{zeroinfl}(
  model.full,
  model.reduced = NULL,
  data,
  coefficients = TRUE,
  anova = TRUE,
  vcovfunc = sandwich::sandwich,
  vcov.args = list(),
  unbiased = TRUE,
  overall = TRUE,
  ...
)

\method{resi_pe}{geeglm}(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  overall = TRUE,
  unbiased = TRUE,
  ...
)

\method{resi_pe}{glmgee}(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  overall = TRUE,
  unbiased = TRUE,
  ...
)

\method{resi_pe}{gee}(model.full, data, unbiased = TRUE, ...)

\method{resi_pe}{lme}(
  model.full,
  anova = TRUE,
  vcovfunc = clubSandwich::vcovCR,
  Anova.args = list(),
  vcov.args = list(),
  ...
)

\method{resi_pe}{lmerMod}(
  model.full,
  anova = TRUE,
  vcovfunc = clubSandwich::vcovCR,
  Anova.args = list(),
  vcov.args = list(),
  ...
)

\method{resi_pe}{glmmTMB}(
  model.full,
  anova = TRUE,
  vcovfunc = clubSandwich::vcovCR,
  Anova.args = list(),
  vcov.args = list(),
  ...
)

\method{resi_pe}{emmGrid}(object, model, N = NULL, unbiased = TRUE, ...)
}
\arguments{
\item{...}{Ignored.}

\item{model.full}{\code{lm, glm, nls, survreg, coxph, hurdle, zeroinfl, gee, geeglm} or \code{lme} model object.}

\item{model.reduced}{Fitted model object of same type as model.full. By default `NULL`; the same model as the full model but only having intercept.}

\item{data}{Data.frame or object coercible to data.frame of model.full data (required for some model types).}

\item{anova}{Logical, whether to produce an Anova table with the RESI columns added. By default = `TRUE`.}

\item{coefficients}{Logical, whether to produce a coefficients (summary) table with the RESI columns added. By default = `TRUE`.}

\item{overall}{Logical, whether to produce an overall Wald test comparing full to reduced model with RESI columns added. By default = `TRUE`.}

\item{vcovfunc}{The variance estimator function for constructing the Wald test statistic. By default, sandwich::vcovHC (the robust (sandwich) variance estimator).}

\item{Anova.args}{List, additional arguments to be passed to Anova function.}

\item{vcov.args}{List, additional arguments to be passed to vcovfunc.}

\item{unbiased}{Logical, whether to use the unbiased or alternative T/Z statistic to RESI conversion. By default, `TRUE`. See details.}

\item{waldtype}{Numeric, indicates which function to use for overall Wald test. 0 (default) = lmtest::waldtest Chi-square, 1 = lmtest::waldtest F, 2 = aod::wald.test}

\item{object}{emmGrid object}

\item{model}{model used to generate emmeans}

\item{N}{sample size}
}
\value{
Returns a list containing RESI point estimates
}
\description{
This function will estimate the robust effect size (RESI) from Vandekar, Tao, & Blume (2020).
The overall RESI is estimated via a Wald test. RESI is (optionally) estimated for each factor in coefficients-style table.
RESI is (optionally) estimated for each variable/interaction in an Anova-style table
for models with existing Anova methods. This function is the building block for the \code{\link{resi}} function.
}
\details{
The Robust Effect Size Index (RESI) is an effect size measure based on M-estimators.
This function is called by \code{\link{resi}} a specified number of times to
form bootstrapped confidence intervals. Called by itself, this function will
only calculate point estimates.

The RESI, denoted as S, is applicable across many model types. It is a unitless
index and can be easily be compared across models. The RESI can also be
converted to Cohen's \emph{d} (\code{\link{S2d}}) under model homoskedasticity.

The RESI is related to the non-centrality parameter
of the test statistic. The RESI estimate is consistent for all four
(Chi-square, F, T, and Z) types of statistics used. The Chi-square and F-based
calculations rely on asymptotic theory, so they may be biased in small samples.
When possible, the T and Z statistics are used. There are two formulas for both
the T and Z statistic conversion. The first (default, unbiased = TRUE)
are based on solving the expected value of the T or Z statistic for the RESI.
The alternative is based on squaring the T or Z statistic and using the
F or Chi-square statistic conversion. Both of these methods are consistent, but
the alternative exhibits a notable amount of finite sample bias. The alternative
may be appealing because its absolute value will be equal to the RESI based on
the F or Chi-square statistic. The RESI based on the Chi-Square and F statistics
is always greater than or equal to 0. The type of statistic
used is listed with the output. See \code{\link{f2S}}, \code{\link{chisq2S}},
\code{\link{t2S}}, and \code{\link{z2S}} for more details on the formulas.

For GEE (\code{geeglm}) models, a longitudinal RESI (L-RESI) and a cross-sectional,
per-measurement RESI (CS-RESI) is estimated. The longitudinal RESI takes the
specified clustering into account, while the cross-sectional RESI is estimated
using a model where each measurement is its own cluster.
}
\section{Methods (by class)}{
\itemize{
\item \code{resi_pe(default)}: RESI point estimation

\item \code{resi_pe(glm)}: RESI point estimation for generalized linear models

\item \code{resi_pe(lm)}: RESI point estimation for linear models

\item \code{resi_pe(nls)}: RESI point estimation for nonlinear least squares models

\item \code{resi_pe(survreg)}: RESI point estimation for survreg

\item \code{resi_pe(coxph)}: RESI point estimation for coxph models

\item \code{resi_pe(hurdle)}: RESI point estimation for hurdle models

\item \code{resi_pe(zeroinfl)}: RESI point estimation for zeroinfl models

\item \code{resi_pe(geeglm)}: RESI point estimation for geeglm object

\item \code{resi_pe(glmgee)}: RESI point estimation for glmgee object

\item \code{resi_pe(gee)}: RESI point estimation for gee object

\item \code{resi_pe(lme)}: RESI point estimation for lme object

\item \code{resi_pe(lmerMod)}: RESI point estimation for lmerMod object

\item \code{resi_pe(glmmTMB)}: RESI point estimation for glmmTMB object - Gaussian only

\item \code{resi_pe(emmGrid)}: RESI point estimation for emmeans object

}}
\examples{
# This function produces point estimates for the RESI. The resi function will
# provide the same point estimates but adds confidence intervals. See resi for
# more detailed examples.

## resi_pe for a linear model
# fit linear model
mod <- lm(charges ~ region * age + bmi + sex, data = RESI::insurance)
# run resi_pe on the model
resi_pe(mod)

# if you want to have RESI estimates in the coefficient table that are equal in absolute
# value to those in the Anova table (except for those with >1 df and/or included in other
# interaction terms), you can specify unbiased = FALSE to use the alternate conversion.
resi_pe(mod, unbiased = FALSE)
}
\references{
Vandekar S, Tao R, Blume J. A Robust Effect Size Index. \emph{Psychometrika}. 2020 Mar;85(1):232-246. doi: 10.1007/s11336-020-09698-2.
}
