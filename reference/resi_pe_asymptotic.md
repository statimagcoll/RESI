# Robust Effect Size Index with Asymptotic Confidence Intervals

Computes RESI point estimates and asymptotic confidence intervals using
either a normal approximation (Zhang et al., 2025) or a quadratic-form
(Imhof/Davies) approach.

## Usage

``` r
resi_pe_asymptotic(
  model.full,
  data,
  vcovfunc = sandwich::vcovHC,
  coefficients = TRUE,
  anova = TRUE,
  alpha = 0.05,
  ci.method = c("normal", "qf", "cf"),
  type = "HC3",
  unbiased = TRUE,
  Anova.args = list(),
  vcov.args = list(),
  ...
)
```

## Arguments

- model.full:

  Fitted `lm` or `glm` model object.

- data:

  Data frame of model data.

- vcovfunc:

  Variance estimator for RESI point estimates. Default:
  [`sandwich::vcovHC`](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html).

- coefficients:

  Logical; include coefficient table. Default `TRUE`.

- anova:

  Logical; include anova table. Default `TRUE`.

- alpha:

  Numeric; significance level. Default `0.05`.

- ci.method:

  Character; `"normal"` (truncated normal), `"qf"` (quadratic-form
  Imhof), or `"cf"` (Cornish-Fisher test inversion). Default `"normal"`.

- type:

  Character; HC type for the M-estimator sandwich used in CI
  construction. Default `"HC3"`.

- unbiased:

  Logical; use bias-corrected RESI point estimate. Default `TRUE`.

- Anova.args:

  List; additional arguments passed to
  [`car::Anova`](https://rdrr.io/pkg/car/man/Anova.html).

- vcov.args:

  List; additional arguments passed to `vcovfunc`.

- ...:

  Ignored.

## Value

A list of class `"resi"` with `coefficients` and/or `anova` tables
containing RESI point estimates and CIs.
