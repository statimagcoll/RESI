# Robust Effect Size Index (RESI) Point Estimation

This function will estimate the robust effect size (RESI) from Vandekar,
Tao, & Blume (2020). The overall RESI is estimated via a Wald test. RESI
is (optionally) estimated for each factor in coefficients-style table.
RESI is (optionally) estimated for each variable/interaction in an
Anova-style table for models with existing Anova methods. This function
is the building block for the
[`resi`](https://statimagcoll.github.io/RESI/reference/resi.md)
function.

## Usage

``` r
resi_pe(...)

# Default S3 method
resi_pe(
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

# S3 method for class 'glm'
resi_pe(
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

# S3 method for class 'lm'
resi_pe(
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

# S3 method for class 'nls'
resi_pe(
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

# S3 method for class 'survreg'
resi_pe(
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

# S3 method for class 'coxph'
resi_pe(
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

# S3 method for class 'hurdle'
resi_pe(
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

# S3 method for class 'zeroinfl'
resi_pe(
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

# S3 method for class 'lmrob'
resi_pe(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  vcovfunc = stats::vcov,
  Anova.args = list(),
  vcov.args = list(),
  unbiased = TRUE,
  overall = TRUE,
  ...
)

# S3 method for class 'glmrob'
resi_pe(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  vcovfunc = stats::vcov,
  Anova.args = list(),
  vcov.args = list(),
  unbiased = TRUE,
  overall = TRUE,
  ...
)

# S3 method for class 'geeglm'
resi_pe(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  overall = TRUE,
  unbiased = TRUE,
  ...
)

# S3 method for class 'glmgee'
resi_pe(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  overall = TRUE,
  unbiased = TRUE,
  ...
)

# S3 method for class 'gee'
resi_pe(model.full, data, unbiased = TRUE, ...)

# S3 method for class 'lme'
resi_pe(
  model.full,
  anova = TRUE,
  vcovfunc = clubSandwich::vcovCR,
  Anova.args = list(),
  vcov.args = list(),
  ...
)

# S3 method for class 'lmerMod'
resi_pe(
  model.full,
  anova = TRUE,
  vcovfunc = clubSandwich::vcovCR,
  Anova.args = list(),
  vcov.args = list(),
  unbiased = TRUE,
  ...
)

# S3 method for class 'glmmTMB'
resi_pe(
  model.full,
  anova = TRUE,
  vcovfunc = clubSandwich::vcovCR,
  Anova.args = list(),
  vcov.args = list(),
  ...
)

# S3 method for class 'emmGrid'
resi_pe(object, model, N = NULL, unbiased = TRUE, ...)
```

## Arguments

- ...:

  Ignored.

- model.full:

  `lm, glm, nls, survreg, coxph, hurdle, zeroinfl, gee, geeglm` or `lme`
  model object.

- model.reduced:

  Fitted model object of same type as model.full. By default \`NULL\`;
  the same model as the full model but only having intercept.

- data:

  Data.frame or object coercible to data.frame of model.full data
  (required for some model types).

- anova:

  Logical, whether to produce an Anova table with the RESI columns
  added. By default = \`TRUE\`.

- coefficients:

  Logical, whether to produce a coefficients (summary) table with the
  RESI columns added. By default = \`TRUE\`.

- overall:

  Logical, whether to produce an overall Wald test comparing full to
  reduced model with RESI columns added. By default = \`TRUE\`.

- vcovfunc:

  The variance estimator function for constructing the Wald test
  statistic. By default, sandwich::vcovHC (the robust (sandwich)
  variance estimator).

- Anova.args:

  List, additional arguments to be passed to Anova function.

- vcov.args:

  List, additional arguments to be passed to vcovfunc.

- unbiased:

  Logical, whether to use the unbiased or alternative T/Z statistic to
  RESI conversion. By default, \`TRUE\`. See details.

- waldtype:

  Numeric, indicates which function to use for overall Wald test. 0
  (default) = lmtest::waldtest Chi-square, 1 = lmtest::waldtest F, 2 =
  aod::wald.test

- object:

  emmGrid object

- model:

  model used to generate emmeans

- N:

  sample size

## Value

Returns a list containing RESI point estimates

## Details

The Robust Effect Size Index (RESI) is an effect size measure based on
M-estimators. This function is called by
[`resi`](https://statimagcoll.github.io/RESI/reference/resi.md) a
specified number of times to form bootstrapped confidence intervals.
Called by itself, this function will only calculate point estimates.

The RESI, denoted as S, is applicable across many model types. It is a
unitless index and can be easily be compared across models. The RESI can
also be converted to Cohen's *d*
([`S2d`](https://statimagcoll.github.io/RESI/reference/S2d.md)) under
model homoskedasticity.

The RESI is related to the non-centrality parameter of the test
statistic. The RESI estimate is consistent for all four (Chi-square, F,
T, and Z) types of statistics used. The Chi-square and F-based
calculations rely on asymptotic theory, so they may be biased in small
samples. When possible, the T and Z statistics are used. There are two
formulas for both the T and Z statistic conversion. The first (default,
unbiased = TRUE) are based on solving the expected value of the T or Z
statistic for the RESI. The alternative is based on squaring the T or Z
statistic and using the F or Chi-square statistic conversion. Both of
these methods are consistent, but the alternative exhibits a notable
amount of finite sample bias. The alternative may be appealing because
its absolute value will be equal to the RESI based on the F or
Chi-square statistic. The RESI based on the Chi-Square and F statistics
is always greater than or equal to 0. The type of statistic used is
listed with the output. See
[`f2S`](https://statimagcoll.github.io/RESI/reference/f2S.md),
[`chisq2S`](https://statimagcoll.github.io/RESI/reference/chisq2S.md),
[`t2S`](https://statimagcoll.github.io/RESI/reference/t2S.md), and
[`z2S`](https://statimagcoll.github.io/RESI/reference/z2S.md) for more
details on the formulas.

For GEE (`geeglm`) models, a longitudinal RESI (L-RESI) and a
cross-sectional, per-measurement RESI (CS-RESI) is estimated. The
longitudinal RESI takes the specified clustering into account, while the
cross-sectional RESI is estimated using a model where each measurement
is its own cluster.

## Methods (by class)

- `resi_pe(default)`: RESI point estimation

- `resi_pe(glm)`: RESI point estimation for generalized linear models

- `resi_pe(lm)`: RESI point estimation for linear models

- `resi_pe(nls)`: RESI point estimation for nonlinear least squares
  models

- `resi_pe(survreg)`: RESI point estimation for survreg

- `resi_pe(coxph)`: RESI point estimation for coxph models

- `resi_pe(hurdle)`: RESI point estimation for hurdle models

- `resi_pe(zeroinfl)`: RESI point estimation for zeroinfl models

- `resi_pe(lmrob)`: RESI point estimation for lmrob objects (robustbase
  package)

- `resi_pe(glmrob)`: RESI point estimation for glmrob objects
  (robustbase package)

- `resi_pe(geeglm)`: RESI point estimation for geeglm object

- `resi_pe(glmgee)`: RESI point estimation for glmgee object

- `resi_pe(gee)`: RESI point estimation for gee object

- `resi_pe(lme)`: RESI point estimation for lme object

- `resi_pe(lmerMod)`: RESI point estimation for lmerMod object

- `resi_pe(glmmTMB)`: RESI point estimation for glmmTMB object -
  Gaussian only

- `resi_pe(emmGrid)`: RESI point estimation for emmeans object

## References

Vandekar S, Tao R, Blume J. A Robust Effect Size Index. *Psychometrika*.
2020 Mar;85(1):232-246. doi: 10.1007/s11336-020-09698-2.

## Examples

``` r
# This function produces point estimates for the RESI. The resi function will
# provide the same point estimates but adds confidence intervals. See resi for
# more detailed examples.

## resi_pe for a linear model
# fit linear model
mod <- lm(charges ~ region * age + bmi + sex, data = RESI::insurance)
# run resi_pe on the model
resi_pe(mod)
#> 
#> Analysis of effect sizes based on RESI:
#> Confidence level = 
#> Call:  lm(formula = charges ~ region * age + bmi + sex, data = RESI::insurance)
#> 
#> Coefficient Table 
#>                       Estimate Std. Error t value Pr(>|t|)    RESI
#> (Intercept)         -5359.4352  2175.9439 -2.4630   0.0139 -0.0673
#> regionnorthwest     -2339.4433  2395.1507 -0.9767   0.3289 -0.0267
#> regionsoutheast     -3230.8512  2643.1099 -1.2224   0.2218 -0.0334
#> regionsouthwest      -232.4839  2574.2823 -0.0903   0.9281 -0.0025
#> age                   220.3325    40.2091  5.4797   0.0000  0.1497
#> bmi                   323.7725    58.0849  5.5741   0.0000  0.1523
#> sexmale              1328.0215   621.7421  2.1360   0.0329  0.0584
#> regionnorthwest:age    34.9040    57.2364  0.6098   0.5421  0.0167
#> regionsoutheast:age    83.6359    63.3258  1.3207   0.1868  0.0361
#> regionsouthwest:age   -33.6290    61.4065 -0.5476   0.5840 -0.0150
#> 
#> 
#> Analysis of Deviance Table (Type II tests)
#> 
#> Response: charges
#>            Df        F Pr(>F)   RESI
#> region      3   1.5959 0.1886 0.0365
#> age         1 117.7046 0.0000 0.2951
#> bmi         1  31.0708 0.0000 0.1498
#> sex         1   4.5624 0.0329 0.0515
#> region:age  3   1.1167 0.3412 0.0161
#> 
#> Overall RESI comparing model to intercept-only model:
#> 
#>   Res.Df Df       F Pr(>F)   RESI
#> 1   1328  9 20.2486      0 0.3595
#> 
#> Notes:
#> 1. The RESI was calculated using a robust covariance estimator.

# if you want to have RESI estimates in the coefficient table that are equal in absolute
# value to those in the Anova table (except for those with >1 df and/or included in other
# interaction terms), you can specify unbiased = FALSE to use the alternate conversion.
resi_pe(mod, unbiased = FALSE)
#> 
#> Analysis of effect sizes based on RESI:
#> Confidence level = 
#> Call:  lm(formula = charges ~ region * age + bmi + sex, data = RESI::insurance)
#> 
#> Coefficient Table 
#>                       Estimate Std. Error t value Pr(>|t|)    RESI
#> (Intercept)         -5359.4352  2175.9439 -2.4630   0.0139 -0.0615
#> regionnorthwest     -2339.4433  2395.1507 -0.9767   0.3289  0.0000
#> regionsoutheast     -3230.8512  2643.1099 -1.2224   0.2218 -0.0192
#> regionsouthwest      -232.4839  2574.2823 -0.0903   0.9281  0.0000
#> age                   220.3325    40.2091  5.4797   0.0000  0.1472
#> bmi                   323.7725    58.0849  5.5741   0.0000  0.1498
#> sexmale              1328.0215   621.7421  2.1360   0.0329  0.0515
#> regionnorthwest:age    34.9040    57.2364  0.6098   0.5421  0.0000
#> regionsoutheast:age    83.6359    63.3258  1.3207   0.1868  0.0235
#> regionsouthwest:age   -33.6290    61.4065 -0.5476   0.5840  0.0000
#> 
#> 
#> Analysis of Deviance Table (Type II tests)
#> 
#> Response: charges
#>            Df        F Pr(>F)   RESI
#> region      3   1.5959 0.1886 0.0365
#> age         1 117.7046 0.0000 0.2951
#> bmi         1  31.0708 0.0000 0.1498
#> sex         1   4.5624 0.0329 0.0515
#> region:age  3   1.1167 0.3412 0.0161
#> 
#> Overall RESI comparing model to intercept-only model:
#> 
#>   Res.Df Df       F Pr(>F)   RESI
#> 1   1328  9 20.2486      0 0.3595
#> 
#> Notes:
#> 1. The RESI was calculated using a robust covariance estimator.
```
