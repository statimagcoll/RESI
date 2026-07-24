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
  ci.method = c("qf", "normal", "cf"),
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

  Character; `"qf"` (quadratic-form Imhof, default), `"normal"`
  (truncated normal), or `"cf"` (Cornish-Fisher test inversion).

- type:

  Character; HC type for the sandwich variance used in CI construction
  (controls tau_i weights in SigmaXw). Default `"HC3"`. Supported:
  `"HC0"`–`"HC5"`, `"const"`.

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
