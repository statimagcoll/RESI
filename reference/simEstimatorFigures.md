# Estimator Comparison Figures: t2S/f2S vs z2S/chisq2S

Reads the per-replicate raw `.rds` files from one simulation directory
and produces Bias + MSE comparison figures contrasting the RESI
estimators actually used for `lm` models against the simpler z-statistic
/ chi-squared alternatives:

- **Coefficients table**:
  [`t2S`](https://statimagcoll.github.io/RESI/reference/t2S.md) (used)
  vs [`z2S`](https://statimagcoll.github.io/RESI/reference/z2S.md)
  (alternative)

- **Anova table**:
  [`f2S`](https://statimagcoll.github.io/RESI/reference/f2S.md) (used)
  vs
  [`chisq2S`](https://statimagcoll.github.io/RESI/reference/chisq2S.md)
  (alternative)

Only `lm` models are included; GLM models are skipped because z2S and
chisq2S are the natural estimators there.

## Usage

``` r
simEstimatorFigures(
  sim.dir = "resiBootSim",
  figures.dir = NULL,
  alpha = 0.05,
  fixed.knots = FALSE
)
```

## Arguments

- sim.dir:

  Character, directory containing simulation output with a `sim_raw/`
  sub-directory. Default `"resiBootSim"`.

- figures.dir:

  Character, output directory for PDFs. Default
  `file.path(sim.dir, "figures", "estimator_compare")`.

- alpha:

  Numeric, nominal level (used only in figure titles). Default 0.05.

- fixed.knots:

  Logical. Must match the value used in the original
  [`insurancePlasmodeSim`](https://statimagcoll.github.io/RESI/reference/insurancePlasmodeSim.md)
  call. Default `FALSE`.

## Value

Invisibly returns the combined estimator-comparison `data.frame`.

## See also

[`simCompareMethodsFigures`](https://statimagcoll.github.io/RESI/reference/simCompareMethodsFigures.md),
[`simRecomputeSummary`](https://statimagcoll.github.io/RESI/reference/simRecomputeSummary.md)
