# Design-averaged per-unit SE for the hybrid plasmode

sqrt(n_k E_X\[V(X)\]\_kk) at g schools and average m students per
school, averaging the conditional covariance V(X) = I(X)^-1 over
\`n_rep\` resampled sets of g school profiles. Because beta-hat is
conditionally unbiased this is the finite-(g, m) reference for the Monte
Carlo SD of beta-hat; it exceeds \`expectedDesignResi()\$se_unit\` since
E\[I^-1\] \>= E\[I\]^-1.

## Usage

``` r
designAveragedSe(
  fit,
  g,
  m,
  n_rep = 2000L,
  school = "schoolid",
  child = "childid:schoolid"
)
```

## Arguments

- fit:

  A \`lmerMod\` fit with a school random intercept and child-level
  random effects (intercept and optionally slopes), e.g. \`y ~ X + (1 \|
  school/child)\`.

- g:

  Number of schools.

- m:

  Average number of students per school.

- n_rep:

  Number of resampled school sets.

- school, child:

  Names of the school and child grouping factors in \`lme4::getME(fit,
  "flist")\`.
