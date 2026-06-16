# Package index

## Model-based functions

Use these functions on fitted models to obtain RESI estimates and
confidence intervals.

- [`resi()`](https://statimagcoll.github.io/RESI/reference/resi.md) :
  Robust Effect Size Index (RESI) point and interval estimation for
  models
- [`resi_pe()`](https://statimagcoll.github.io/RESI/reference/resi_pe.md)
  : Robust Effect Size Index (RESI) Point Estimation

## Conversion functions

Manual functions for converting test statistics to RESI estimates and
converting RESI estimates to and from other effect size measures.

- [`chisq2S()`](https://statimagcoll.github.io/RESI/reference/chisq2S.md)
  : Compute the robust effect size index estimate from chi-squared
  statistic.

- [`d2S()`](https://statimagcoll.github.io/RESI/reference/d2S.md) :

  Covert Cohen's *d* to \|S\|

- [`f2S()`](https://statimagcoll.github.io/RESI/reference/f2S.md) :
  Compute the robust effect size index estimate from F-statistic

- [`fsq2S()`](https://statimagcoll.github.io/RESI/reference/fsq2S.md) :

  Covert Cohen's *f*^2 to S

- [`Rsq2S()`](https://statimagcoll.github.io/RESI/reference/Rsq2S.md) :
  Covert R^2 to S

- [`S2d()`](https://statimagcoll.github.io/RESI/reference/S2d.md) :
  Convert S to Cohen's d

- [`S2fsq()`](https://statimagcoll.github.io/RESI/reference/S2fsq.md) :

  Covert S to Cohen's *f*^2

- [`S2Rsq()`](https://statimagcoll.github.io/RESI/reference/S2Rsq.md) :
  Covert S to R^2

- [`S2chisq()`](https://statimagcoll.github.io/RESI/reference/S2chisq.md)
  : Convert non-zero S to Chi-square statistic

- [`S2z()`](https://statimagcoll.github.io/RESI/reference/S2z.md) :
  Convert RESI (S) estimate to Z statistic

- [`t2S()`](https://statimagcoll.github.io/RESI/reference/t2S.md) :
  Compute the robust effect size index estimate from t statistic
  (default)

- [`z2S()`](https://statimagcoll.github.io/RESI/reference/z2S.md) :
  Compute the robust effect size index estimate from Z statistic

## Additional functions

Functions for summarizing and visualization.

- [`anova(`*`<resi>`*`)`](https://statimagcoll.github.io/RESI/reference/anova.resi.md)
  : Anova method for resi objects
- [`ggplot(`*`<resi>`*`)`](https://statimagcoll.github.io/RESI/reference/ggplot.resi.md)
  : Plotting RESI Estimates and CIs
- [`omnibus()`](https://statimagcoll.github.io/RESI/reference/omnibus.md)
  : Omnibus (Overall) Wald Test for resi objects
- [`plot(`*`<resi>`*`)`](https://statimagcoll.github.io/RESI/reference/plot.resi.md)
  : Plotting RESI Estimates and CIs
- [`summary(`*`<resi>`*`)`](https://statimagcoll.github.io/RESI/reference/summary.resi.md)
  : Summary method for resi objects

## Datasets

- [`depression`](https://statimagcoll.github.io/RESI/reference/depression.md)
  : Depression Treatment Data
- [`insurance`](https://statimagcoll.github.io/RESI/reference/insurance.md)
  : US Health Insurance Data
