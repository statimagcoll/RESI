# Lyu-Welsh limiting RESI for a three-level lmer fit

Evaluates the g, m -\> Inf limit (Lyu & Welsh, 2022) of the per-unit
covariance of the fixed effects and the corresponding RESI at the fitted
parameters. Between-school columns: n_g Var(beta_b) -\> tau2_school
C^-1, C = g^-1 sum_s x_s x_s'. Within-school columns: the school
intercept acts as a fixed effect, so I_w / N = N^-1 sum_s (A_s - b_s
b_s' / c_s). The between and within blocks are asymptotically
independent.

## Usage

``` r
lyuLimitResi(fit, school = "schoolid", child = "childid:schoolid")
```

## Arguments

- fit:

  A \`lmerMod\` fit with a school random intercept and child-level
  random effects (intercept and optionally slopes), e.g. \`y ~ X + (1 \|
  school/child)\`.

- school, child:

  Names of the school and child grouping factors in \`lme4::getME(fit,
  "flist")\`.

## Value

List with \`coefficients\` and \`anova\` RESI tables, the logical
\`between\` indicator per column, and the per-unit limiting covariance
\`V\`.
