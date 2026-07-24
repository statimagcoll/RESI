# Robust Effect Size Index (RESI) point and interval estimation for models

This function will estimate the robust effect size (RESI) from Vandekar,
Tao, & Blume (2020) and its confidence interval in various ways for a
fitted model object. The overall RESI is estimated via a Wald test. RESI
is (optionally) estimated for each factor in coefficients-style table.
RESI is (optionally) estimated for each variable/interaction in an
Anova-style table for models with existing Anova methods. CIs can be
calculated using either non-parametric or Bayesian bootstrapping.

## Usage

``` r
resi(model.full, ...)

# Default S3 method
resi(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  overall = TRUE,
  nboot = 1000,
  boot.method = "nonparam",
  vcovfunc = sandwich::vcovHC,
  alpha = 0.05,
  store.boot = FALSE,
  Anova.args = list(),
  vcov.args = list(),
  unbiased = TRUE,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("boot.ncpus", 1L),
  long = FALSE,
  clvar = NULL,
  ci.method = c("boot", "qf", "cf", "normal"),
  ...
)

# S3 method for class 'glm'
resi(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  overall = TRUE,
  nboot = 1000,
  vcovfunc = sandwich::vcovHC,
  alpha = 0.05,
  store.boot = FALSE,
  Anova.args = list(),
  vcov.args = list(),
  unbiased = TRUE,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("boot.ncpus", 1L),
  ci.method = "qf",
  ...
)

# S3 method for class 'lm'
resi(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  overall = TRUE,
  nboot = 1000,
  boot.method = "nonparam",
  vcovfunc = sandwich::vcovHC,
  alpha = 0.05,
  store.boot = FALSE,
  Anova.args = list(),
  vcov.args = list(),
  unbiased = TRUE,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("boot.ncpus", 1L),
  ci.method = "qf",
  ...
)

# S3 method for class 'lmrob'
resi(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  overall = TRUE,
  nboot = 1000,
  boot.method = "nonparam",
  vcovfunc = stats::vcov,
  alpha = 0.05,
  store.boot = FALSE,
  Anova.args = list(),
  vcov.args = list(),
  unbiased = TRUE,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("boot.ncpus", 1L),
  ci.method = "boot",
  ...
)

# S3 method for class 'glmrob'
resi(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  overall = TRUE,
  nboot = 1000,
  vcovfunc = stats::vcov,
  alpha = 0.05,
  store.boot = FALSE,
  Anova.args = list(),
  vcov.args = list(),
  unbiased = TRUE,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("boot.ncpus", 1L),
  ci.method = "boot",
  ...
)

# S3 method for class 'nls'
resi(
  model.full,
  model.reduced = NULL,
  data,
  coefficients = TRUE,
  overall = TRUE,
  nboot = 1000,
  boot.method = "nonparam",
  anova = FALSE,
  vcovfunc = r_nlshc,
  alpha = 0.05,
  store.boot = FALSE,
  vcov.args = list(),
  unbiased = TRUE,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("boot.ncpus", 1L),
  ...
)

# S3 method for class 'survreg'
resi(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  overall = TRUE,
  nboot = 1000,
  vcovfunc = vcov,
  alpha = 0.05,
  store.boot = FALSE,
  Anova.args = list(),
  unbiased = TRUE,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("boot.ncpus", 1L),
  ...
)

# S3 method for class 'coxph'
resi(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  overall = TRUE,
  nboot = 1000,
  vcovfunc = vcov,
  alpha = 0.05,
  store.boot = FALSE,
  Anova.args = list(),
  unbiased = TRUE,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("boot.ncpus", 1L),
  ...
)

# S3 method for class 'hurdle'
resi(
  model.full,
  model.reduced = NULL,
  data,
  coefficients = TRUE,
  overall = TRUE,
  nboot = 1000,
  vcovfunc = sandwich::sandwich,
  anova = FALSE,
  alpha = 0.05,
  store.boot = FALSE,
  vcov.args = list(),
  unbiased = TRUE,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("boot.ncpus", 1L),
  ...
)

# S3 method for class 'zeroinfl'
resi(
  model.full,
  model.reduced = NULL,
  data,
  coefficients = TRUE,
  overall = TRUE,
  nboot = 1000,
  vcovfunc = sandwich::sandwich,
  anova = FALSE,
  alpha = 0.05,
  store.boot = FALSE,
  vcov.args = list(),
  unbiased = TRUE,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("boot.ncpus", 1L),
  ...
)

# S3 method for class 'geeglm'
resi(
  model.full,
  model.reduced = NULL,
  data,
  anova = TRUE,
  coefficients = TRUE,
  overall = TRUE,
  nboot = 1000,
  alpha = 0.05,
  store.boot = FALSE,
  unbiased = TRUE,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("boot.ncpus", 1L),
  ...
)

# S3 method for class 'glmgee'
resi(
  model.full,
  model.reduced = NULL,
  data,
  anova = FALSE,
  coefficients = TRUE,
  overall = TRUE,
  nboot = 1000,
  alpha = 0.05,
  store.boot = FALSE,
  unbiased = TRUE,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("boot.ncpus", 1L),
  ...
)

# S3 method for class 'gee'
resi(
  model.full,
  data,
  nboot = 1000,
  alpha = 0.05,
  store.boot = FALSE,
  unbiased = TRUE,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("boot.ncpus", 1L),
  ...
)

# S3 method for class 'lme'
resi(
  model.full,
  alpha = 0.05,
  nboot = 1000,
  anova = TRUE,
  vcovfunc = clubSandwich::vcovCR,
  vcov.args = list(),
  ...
)

# S3 method for class 'lmerMod'
resi(
  model.full,
  alpha = 0.05,
  nboot = 1000,
  anova = TRUE,
  vcovfunc = clubSandwich::vcovCR,
  vcov.args = list(),
  ...
)

# S3 method for class 'glmmTMB'
resi(
  model.full,
  alpha = 0.05,
  nboot = 1000,
  anova = TRUE,
  vcovfunc = clubSandwich::vcovCR,
  vcov.args = list(),
  ...
)
```

## Arguments

- model.full:

  `lm, glm, nls, survreg, coxph, hurdle, zeroinfl, gee, geeglm` or `lme`
  model object.

- ...:

  Ignored.

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

- nboot:

  Numeric, the number of bootstrap replicates. Used only when
  \`ci.method = "boot"\`. By default, 1000.

- boot.method:

  String, which type of bootstrap to use: \`nonparam\` = non-parametric
  bootstrap (default); \`bayes\` = Bayesian bootstrap.

- vcovfunc:

  The variance estimator function for constructing the Wald test
  statistic. By default,
  [vcovHC](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html)
  (the robust (sandwich) variance estimator).

- alpha:

  Numeric, significance level of the constructed CIs. By default, 0.05.

- store.boot:

  Logical, whether to store all the bootstrapped estimates. By default,
  \`FALSE\`.

- Anova.args:

  List, additional arguments to be passed to
  [Anova](https://rdrr.io/pkg/car/man/Anova.html) function.

- vcov.args:

  List, additional arguments to be passed to vcovfunc.

- unbiased:

  Logical, whether to use the unbiased or alternative T/Z statistic to
  RESI conversion. By default, \`TRUE\`. See details.

- parallel:

  See documentation for [boot](https://rdrr.io/pkg/boot/man/boot.html).

- ncpus:

  See documentation for [boot](https://rdrr.io/pkg/boot/man/boot.html).

- long:

  Logical, whether the data is longitudinal/clustered. By default,
  \`FALSE\`.

- clvar:

  Character, the name of the cluster/id variable if data is clustered.
  By default, \`NULL\`.

- ci.method:

  Character, the method used to compute confidence intervals. One of
  \`"boot"\` (bootstrap), \`"qf"\` (quadratic-form inversion, default
  for lm and glm), \`"normal"\` (normal approximation), or \`"cf"\`
  (Cornish-Fisher inversion). See \`resi_pe_asymptotic\` for details on
  the asymptotic methods.

## Value

Returns a list that includes function arguments, RESI point estimates,
and confidence intervals in coefficients/anova-style tables

## Details

The RESI, denoted as S, is applicable across many model types. It is a
unitless index and can be easily be compared across models. The RESI can
also be converted to Cohen's *d*
([`S2d`](https://statimagcoll.github.io/RESI/reference/S2d.md)) under
model homoskedasticity.

This function computes the RESI point estimates and bootstrapped
confidence intervals based on Chi-square, F, T, or Z statistics. The
robust (sandwich) variance is used by default, allowing for consistency
under model-misspecification. The RESI is related to the non-centrality
parameter of the test statistic. The RESI estimate is consistent for all
four (Chi-square, F, T, and Z) types of statistics used. The Chi-square
and F-based calculations rely on asymptotic theory, so they may be
biased in small samples. When possible, the T and Z statistics are used.
There are two formulas for both the T and Z statistic conversion. The
first (default, unbiased = TRUE) are based on solving the expected value
of the T or Z statistic for the RESI. The alternative is based on
squaring the T or Z statistic and using the F or Chi-square statistic
conversion. Both of these methods are consistent, but the alternative
exhibits a notable amount of finite sample bias. The alternative may be
appealing because its absolute value will be equal to the RESI based on
the F or Chi-square statistic. The RESI based on the Chi-Square and F
statistics is always greater than or equal to 0. The type of statistic
used is listed with the output. See
[`f2S`](https://statimagcoll.github.io/RESI/reference/f2S.md),
[`chisq2S`](https://statimagcoll.github.io/RESI/reference/chisq2S.md),
[`t2S`](https://statimagcoll.github.io/RESI/reference/t2S.md), and
[`z2S`](https://statimagcoll.github.io/RESI/reference/z2S.md) for more
details on the formulas.

For GEE ([geeglm](https://rdrr.io/pkg/geepack/man/geeglm.html),
[glmgee](https://rdrr.io/pkg/glmtoolbox/man/glmgee.html)) models, a
longitudinal RESI (L-RESI) and a cross-sectional, per-measurement RESI
(CS-RESI) is estimated. The longitudinal RESI takes the specified
clustering into account, while the cross-sectional RESI is designed to
estimate the effect size if a random observation for each participant
were collected cross-sectionally.

For most `lm` and `nls` model types, there is a Bayesian bootstrap
option available as an alternative to the default, standard
non-parametric bootstrap. The interpretation of a Bayesian bootstrapped
interval is similar to that of a credible interval.

Certain model types require the data used for the model be entered as an
argument. These are: `nls, survreg,` and `coxph`. Additionally, if a
model includes certain functions (splines, factor, I), the data needs to
be provided.

If running into convergence issues with nls models, it is advised to
refit the nls model with starting values equal to the estimates provided
by the model and then try rerunning `resi`.

## Methods (by class)

- `resi(default)`: RESI point and interval estimation for models

- `resi(glm)`: RESI point and interval estimation for models

- `resi(lm)`: RESI point and interval estimation for lm models

- `resi(lmrob)`: RESI point and interval estimation for lmrob models
  (robustbase)

- `resi(glmrob)`: RESI point and interval estimation for glmrob models
  (robustbase)

- `resi(nls)`: RESI point and interval estimation for nls models

- `resi(survreg)`: RESI point and interval estimation for survreg models

- `resi(coxph)`: RESI point and interval estimation for coxph models

- `resi(hurdle)`: RESI point and interval estimation for hurdle models

- `resi(zeroinfl)`: RESI point and interval estimation for zeroinfl
  models

- `resi(geeglm)`: RESI point and interval estimation for GEE models

- `resi(glmgee)`: RESI point and interval estimation for GEE models in
  glmtoolbox

- `resi(gee)`: RESI point and interval estimation for GEE models

- `resi(lme)`: RESI point and interval estimation for LME (nlme) models

- `resi(lmerMod)`: RESI point and interval estimation for lmerMod models

- `resi(glmmTMB)`: RESI point and interval estimation for glmmTMB
  models - Gaussian only

## References

Vandekar S, Tao R, Blume J. A Robust Effect Size Index. *Psychometrika*.
2020 Mar;85(1):232-246. doi: 10.1007/s11336-020-09698-2.

Kang, K., Armstrong, K., Avery, S., McHugo, M., Heckers, S., & Vandekar,
S. (2021). Accurate confidence interval estimation for non-centrality
parameters and effect size indices. *arXiv preprint arXiv:2111.05966*.

Jones, M., Kang, K., Vandekar, S. (2025). *Journal of Statistical
Software*. RESI: An R Package for Robust Effect
Sizes.\<doi:10.18637/jss.v112.i03\>

## See also

[`resi_pe`](https://statimagcoll.github.io/RESI/reference/resi_pe.md),
[vcovHC](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html),
[`f2S`](https://statimagcoll.github.io/RESI/reference/f2S.md),
[`chisq2S`](https://statimagcoll.github.io/RESI/reference/chisq2S.md),
[`z2S`](https://statimagcoll.github.io/RESI/reference/z2S.md),
[`t2S`](https://statimagcoll.github.io/RESI/reference/t2S.md)

## Examples

``` r
## for timing purposes, a small number of bootstrap replicates is used in the
## examples. Run them with a higher or default `nboot` argument for better performance

## RESI on a linear model
# fit linear model
mod = lm(charges ~ region * age + bmi + sex, data = RESI::insurance)

# run resi on fitted model with desired number of bootstrap replicates
# store bootstrap results for calculating different CIs later
resi_obj = resi(mod, nboot = 50, store.boot = TRUE)
# print output
resi_obj
#> 
#> Analysis of effect sizes based on RESI:
#> Confidence level =  0.05
#> Call:  lm(formula = charges ~ region * age + bmi + sex, data = RESI::insurance)
#> 
#> Coefficient Table 
#>                       Estimate Std. Error t value Pr(>|t|)    RESI    2.5%
#> (Intercept)         -5359.4352  2175.9439 -2.4630   0.0139 -0.0673 -0.1203
#> regionnorthwest     -2339.4433  2395.1507 -0.9767   0.3289 -0.0267 -0.0804
#> regionsoutheast     -3230.8512  2643.1099 -1.2224   0.2218 -0.0334 -0.0870
#> regionsouthwest      -232.4839  2574.2823 -0.0903   0.9281 -0.0025 -0.0560
#> age                   220.3325    40.2091  5.4797   0.0000  0.1497  0.0924
#> bmi                   323.7725    58.0849  5.5741   0.0000  0.1523  0.1067
#> sexmale              1328.0215   621.7421  2.1360   0.0329  0.0584  0.0047
#> regionnorthwest:age    34.9040    57.2364  0.6098   0.5421  0.0167 -0.0369
#> regionsoutheast:age    83.6359    63.3258  1.3207   0.1868  0.0361 -0.0175
#> regionsouthwest:age   -33.6290    61.4065 -0.5476   0.5840 -0.0150 -0.0688
#>                       97.5%
#> (Intercept)         -0.0144
#> regionnorthwest      0.0270
#> regionsoutheast      0.0201
#> regionsouthwest      0.0511
#> age                  0.2072
#> bmi                  0.1980
#> sexmale              0.1121
#> regionnorthwest:age  0.0703
#> regionsoutheast:age  0.0898
#> regionsouthwest:age  0.0389
#> 
#> 
#> Analysis of Deviance Table (Type II tests)
#> 
#> Response: charges
#>            Df        F Pr(>F)   RESI   2.5%  97.5%
#> region      3   1.5959 0.1886 0.0365 0.0000 0.1048
#> age         1 117.7046 0.0000 0.2951 0.2384 0.3548
#> bmi         1  31.0708 0.0000 0.1498 0.1067 0.1980
#> sex         1   4.5624 0.0329 0.0515 0.0000 0.1121
#> region:age  3   1.1167 0.3412 0.0161 0.0000 0.0928
#> 
#> Overall RESI comparing model to intercept-only model:
#> 
#>   Res.Df Df       F Pr(>F)   RESI
#> 1   1328  9 20.2486      0 0.3595
#> 
#> Notes:
#> 1. The RESI was calculated using a robust covariance estimator.
#> 2. Confidence intervals (CIs) constructed using 50 non-parametric bootstraps. 
#> 

# fit a reduced model for comparison
mod_red = lm(charges ~ bmi, data = RESI::insurance)

# running resi and including the reduced model will provide almost the exact same
# output as not including a reduced model. The difference is that the "overall"
# portion of the output will compare the full model to the reduced model.
# The "summary" and "anova" RESI estimates will be the same. (The bootstrapped
# confidence intervals may differ.)
resi(model.full = mod, model.reduced = mod_red, nboot = 10)
#> 
#> Analysis of effect sizes based on RESI:
#> Confidence level =  0.05
#> Full Model:lm(formula = charges ~ region * age + bmi + sex, data = RESI::insurance)
#> Reduced Model:lm(formula = charges ~ bmi, data = RESI::insurance)
#> 
#> 
#> Coefficient Table 
#>                       Estimate Std. Error t value Pr(>|t|)    RESI    2.5%
#> (Intercept)         -5359.4352  2175.9439 -2.4630   0.0139 -0.0673 -0.1203
#> regionnorthwest     -2339.4433  2395.1507 -0.9767   0.3289 -0.0267 -0.0804
#> regionsoutheast     -3230.8512  2643.1099 -1.2224   0.2218 -0.0334 -0.0870
#> regionsouthwest      -232.4839  2574.2823 -0.0903   0.9281 -0.0025 -0.0560
#> age                   220.3325    40.2091  5.4797   0.0000  0.1497  0.0924
#> bmi                   323.7725    58.0849  5.5741   0.0000  0.1523  0.1067
#> sexmale              1328.0215   621.7421  2.1360   0.0329  0.0584  0.0047
#> regionnorthwest:age    34.9040    57.2364  0.6098   0.5421  0.0167 -0.0369
#> regionsoutheast:age    83.6359    63.3258  1.3207   0.1868  0.0361 -0.0175
#> regionsouthwest:age   -33.6290    61.4065 -0.5476   0.5840 -0.0150 -0.0688
#>                       97.5%
#> (Intercept)         -0.0144
#> regionnorthwest      0.0270
#> regionsoutheast      0.0201
#> regionsouthwest      0.0511
#> age                  0.2072
#> bmi                  0.1980
#> sexmale              0.1121
#> regionnorthwest:age  0.0703
#> regionsoutheast:age  0.0898
#> regionsouthwest:age  0.0389
#> 
#> 
#> Analysis of Deviance Table (Type II tests)
#> 
#> Response: charges
#>            Df        F Pr(>F)   RESI   2.5%  97.5%
#> region      3   1.5959 0.1886 0.0365 0.0000 0.1048
#> age         1 117.7046 0.0000 0.2951 0.2384 0.3548
#> bmi         1  31.0708 0.0000 0.1498 0.1067 0.1980
#> sex         1   4.5624 0.0329 0.0515 0.0000 0.1121
#> region:age  3   1.1167 0.3412 0.0161 0.0000 0.0928
#> 
#> Overall RESI comparing full model to reduced model:
#> 
#>   Res.Df Df       F Pr(>F)   RESI
#> 1   1328  8 15.9113      0 0.2983
#> 
#> Notes:
#> 1. The RESI was calculated using a robust covariance estimator.
#> 2. Confidence intervals (CIs) constructed using 10 non-parametric bootstraps. 
#> 

# used stored bootstrap results to get a different alpha-level confidence interval
summary(resi_obj, alpha = c(0.01, 0.1))
#> 
#> Analysis of effect sizes based on RESI:
#> Confidence level =  0.01 0.1
#> Call:  lm(formula = charges ~ region * age + bmi + sex, data = RESI::insurance)
#> 
#> Coefficient Table 
#>                       Estimate Std. Error t value Pr(>|t|)    RESI    0.5%
#> (Intercept)         -5359.4352  2175.9439 -2.4630   0.0139 -0.0673 -0.1369
#> regionnorthwest     -2339.4433  2395.1507 -0.9767   0.3289 -0.0267 -0.0973
#> regionsoutheast     -3230.8512  2643.1099 -1.2224   0.2218 -0.0334 -0.1038
#> regionsouthwest      -232.4839  2574.2823 -0.0903   0.9281 -0.0025 -0.0729
#> age                   220.3325    40.2091  5.4797   0.0000  0.1497  0.0743
#> bmi                   323.7725    58.0849  5.5741   0.0000  0.1523  0.0924
#> sexmale              1328.0215   621.7421  2.1360   0.0329  0.0584 -0.0122
#> regionnorthwest:age    34.9040    57.2364  0.6098   0.5421  0.0167 -0.0538
#> regionsoutheast:age    83.6359    63.3258  1.3207   0.1868  0.0361 -0.0344
#> regionsouthwest:age   -33.6290    61.4065 -0.5476   0.5840 -0.0150 -0.0858
#>                      99.5%      5%     95%
#> (Intercept)         0.0022 -0.1118 -0.0229
#> regionnorthwest     0.0439 -0.0718  0.0184
#> regionsoutheast     0.0370 -0.0784  0.0115
#> regionsouthwest     0.0679 -0.0474  0.0425
#> age                 0.2253  0.1016  0.1980
#> bmi                 0.2124  0.1141  0.1907
#> sexmale             0.1290  0.0133  0.1035
#> regionnorthwest:age 0.0871 -0.0283  0.0617
#> regionsoutheast:age 0.1066 -0.0089  0.0811
#> regionsouthwest:age 0.0558 -0.0602  0.0302
car::Anova(resi_obj, alpha = c(0.01, 0.1))
#>            Df        F  Pr(>F)     RESI     0.5%   99.5%       5%     95%
#> region      3   1.5959 0.18856 0.036480 0.000000 0.12289 0.000000 0.09552
#> age         1 117.7046 0.00000 0.295111 0.220069 0.37313 0.247731 0.34547
#> bmi         1  31.0708 0.00000 0.149798 0.092383 0.21239 0.114071 0.19071
#> sex         1   4.5624 0.03286 0.051549 0.000000 0.12902 0.011851 0.10349
#> region:age  3   1.1167 0.34115 0.016056 0.000000 0.11090 0.000000 0.08339

# the result of resi, as well as the summary or Anova of a `resi` object can be plotted
# if the resi object was created with the store.boot = `TRUE` option, any alpha
# can be specified
plot(resi_obj, alpha = 0.01)

# if the variable names on the y-axis are too long, you can reduce their size with
# the ycex.axis argument (or use regular common solutions like changing the margins)
plot(resi_obj, alpha = 0.01, ycex.axis = 0.5)


# for some model types and formula structures, data argument is required
if(requireNamespace("splines")){
  # fit logistic regression model with splines
  mod = glm(smoker ~ splines::ns(age, df = 3) + region, data = RESI::insurance,
    family = "binomial")

  # specify additional arguments to the variance-covariance function via vcov.args
  resi_obj = resi(mod, data = RESI::insurance, alpha = 0.01,
    vcov.args = list(type = "HC0"), nboot = 25)
  summary(resi_obj)
  car::Anova(resi_obj)}
#> Analysis of Deviance Table (Type II tests)
#> 
#> Response: smoker
#>                          Df  Chisq Pr(>Chisq)     RESI 0.5%    99.5%
#> splines::ns(age, df = 3)  3 1.4735    0.68841 0.000000    0 0.089389
#> region                    3 7.2960    0.06304 0.056663    0 0.135462


## RESI on a survival model with alternate Z2S
if(requireNamespace("survival")){
  # fit coxph model on example data from survival package
  # Note: for survival models, you need to specify robust variance in the model
  # creation. resi will ignore the vcovfunc argument for this reason.
  mod.coxph =  survival::coxph(survival::Surv(time, status) ~ age + sex + wt.loss,
   data=survival::lung, robust = TRUE)

  # run resi on the model
  # to use the alternative Z to RESI formula (which is equal in absolute value to the
  # chi-square to RESI (S) formula), specify unbiased = FALSE.
  resi(mod.coxph, data = survival::lung, unbiased = FALSE, nboot = 10)}
#> Loading required namespace: survival
#> 
#> Analysis of effect sizes based on RESI:
#> Confidence level =  0.05
#> Call:  survival::coxph(formula = survival::Surv(time, status) ~ age + 
#>     sex + wt.loss, data = survival::lung, robust = TRUE)
#> 
#> Coefficient Table 
#>         Estimate Std. Error z value Pr(>|z|)    RESI    2.5%   97.5%
#> age       0.0201     0.0101  1.9915   0.0464  0.1177  0.0082  0.3016
#> sex      -0.5210     0.1670 -3.1202   0.0018 -0.2020 -0.3004 -0.0520
#> wt.loss   0.0008     0.0060  0.1264   0.8994  0.0000 -0.0827  0.0835
#> 
#> 
#> Analysis of Deviance Table (Type II tests)
#> 
#> Response: survival::Surv(time, status)
#>         Df  Chisq Pr(>Chisq)   RESI   2.5%  97.5%
#> age      1 3.9662     0.0464 0.1177 0.0082 0.3016
#> sex      1 9.7357     0.0018 0.2020 0.0520 0.3004
#> wt.loss  1 0.0160     0.8994 0.0000 0.0000 0.1021
#> 
#> Overall RESI comparing model to intercept-only model:
#> 
#>      chi2 df      P   RESI   2.5%  97.5%
#> 1 13.7717  3 0.0032 0.2244 0.0994 0.3694
#> 
#> Notes:
#> 1. The RESI was calculated using a robust covariance estimator.
#> 2. Confidence intervals (CIs) constructed using 10 non-parametric bootstraps. 
#> 
```
