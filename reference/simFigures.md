# Simulation Performance Figures for RESI Evaluation

Produces performance figures from the output of
[`insurancePlasmodeSim`](https://statimagcoll.github.io/RESI/reference/insurancePlasmodeSim.md).
Creates one PDF figure per (model type) \\\times\\ (variance estimator)
combination (4 figures total by default). Each figure is a 2 \\\times\\
4 panel grid:

- **Top row**: Anova-table RESI metrics (Bias, MSE, CI Coverage, CI
  Width)

- **Bottom row**: Coefficients-table RESI metrics (same metrics;
  intercept excluded)

- **Colors**: one colored line per model term (high-contrast matte
  palette)

- **x-axis**: sample size on a log scale

Dashed reference lines are drawn at zero for Bias and at \\1 - \alpha\\
for CI Coverage.

## Usage

``` r
simFigures(output.dir = "resiBootSim", alpha = 0.05, ci.label = NULL)
```

## Arguments

- output.dir:

  Character, directory containing simulation output from
  [`insurancePlasmodeSim`](https://statimagcoll.github.io/RESI/reference/insurancePlasmodeSim.md).
  Default `"resiBootSim"`.

- alpha:

  Numeric, nominal CI level used for the coverage reference line.
  Default 0.05.

- ci.label:

  Character, label for the CI method used in figure titles and file
  names. Default `NULL`, which auto-detects from `output.dir`: `"boot"`
  if the directory name contains "Boot"/"boot", `"normal"` if it
  contains "Normal"/"normal", `"qf"` if it contains "QF"/"qf", otherwise
  `basename(output.dir)`.

## Value

Invisibly returns the summary metrics `data.frame`. Saves PDF figures to
`file.path(output.dir, "figures")`, named
`sim_<model>_<vcov>_<ci.label>.pdf`.

## See also

[`insurancePlasmodeSim`](https://statimagcoll.github.io/RESI/reference/insurancePlasmodeSim.md)
