# Fixed-m RESI under the expected hybrid plasmode design

Exact RESI and per-unit SE when each school profile is represented in
proportion to its sampling probability with average \`m\` students per
school (the number of schools g cancels). Converges to
\[lyuLimitResi()\] as m -\> Inf; equals \`resi_pe(level = 2)\` on the
observed design at m = mean school size.

## Usage

``` r
expectedDesignResi(fit, m, school = "schoolid", child = "childid:schoolid")
```

## Arguments

- fit:

  A \`lmerMod\` fit with a school random intercept and child-level
  random effects (intercept and optionally slopes), e.g. \`y ~ X + (1 \|
  school/child)\`.

- m:

  Average number of students per school.

- school, child:

  Names of the school and child grouping factors in \`lme4::getME(fit,
  "flist")\`.

## Value

As \[lyuLimitResi()\] plus \`se_unit\`, sqrt(n_k Var_kk).
