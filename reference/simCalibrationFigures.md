# Calibration Figures for RESI Variance Estimates

Reads output from
[`simCalibrationSim`](https://statimagcoll.github.io/RESI/reference/simCalibrationSim.md)
and produces one PDF per model setting showing convergence of point
estimates, covariance estimates, RESI estimates, and RESI variance
estimates to their population targets.

## Usage

``` r
simCalibrationFigures(
  output.dir = "resiCalibrationSim",
  figures.dir = NULL,
  alpha = 0.05
)
```

## Arguments

- output.dir:

  Character, directory containing output from
  [`simCalibrationSim`](https://statimagcoll.github.io/RESI/reference/simCalibrationSim.md).
  Default `"resiCalibrationSim"`.

- figures.dir:

  Character, output directory for PDFs. Default
  `file.path(output.dir, "figures")`.

- alpha:

  Numeric, nominal CI level. Default 0.05.

## Value

Invisibly returns the summary data frame. Saves PDF figures to
`figures.dir`.

## See also

[`simCalibrationSim`](https://statimagcoll.github.io/RESI/reference/simCalibrationSim.md),
[`insurancePlasmodeSim`](https://statimagcoll.github.io/RESI/reference/insurancePlasmodeSim.md)
