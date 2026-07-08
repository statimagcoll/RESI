# Per-Term CI Method Comparison Figures

Reads simulation summary tables from multiple output directories (one
per CI method) and produces PDF figures comparing CI methods across
sample sizes. For each (model type x variance estimator x table)
combination two PDFs are saved: an estimator comparison (Bias, MSE) and
a CI comparison (SE calibration, coverage, width).

## Usage

``` r
simCompareMethodsFigures(
  output.dirs = c(Bootstrap = "resiBootSim", Normal = "resiAsympNormalSim", QF =
    "resiAsympQFSim", CF = "resiAsympCFSim"),
  figures.dir = NULL,
  alpha = 0.05,
  fixed.knots = FALSE
)
```

## Arguments

- output.dirs:

  Named character vector mapping CI method labels to their simulation
  output directories. Default:
  `c(boot = "resiBootSim", normal = "resiAsympNormalSim", qf = "resiAsympQFSim", cf = "resiAsympCFSim")`.
  Directories that do not exist are silently skipped.

- figures.dir:

  Character, directory where comparison figures are saved. Default:
  `file.path(output.dirs[1], "figures", "method_comparison")`.

- alpha:

  Numeric, nominal CI level used for coverage reference lines. Default
  0.05.

- fixed.knots:

  Logical. Must match the value used in the original
  [`insurancePlasmodeSim`](https://statimagcoll.github.io/RESI/reference/insurancePlasmodeSim.md)
  call. Default `FALSE`.

## Value

Invisibly returns the combined summary `data.frame`. Saves
`estimator_<model>_<vcov>_<table>.pdf`,
`compare_<model>_<vcov>_<table>.pdf`, and
`covquant_<model>_<vcov>_<table>.pdf` to `figures.dir`.

## See also

[`insurancePlasmodeSim`](https://statimagcoll.github.io/RESI/reference/insurancePlasmodeSim.md),
[`simFigures`](https://statimagcoll.github.io/RESI/reference/simFigures.md),
[`simEstimatorFigures`](https://statimagcoll.github.io/RESI/reference/simEstimatorFigures.md)
