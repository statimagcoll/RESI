# Asymptotic Calibration Check for RESI Variance Estimates

Runs a plasmode simulation to verify that the asymptotic normal CI
machinery is correctly calibrated for all four model settings (lm/glm x
parametric/robust). For each (setting, sample size) cell the function
checks:

1.  **Bias(theta)**: \\\hat\theta \to \theta\_{\rm true}\\ – raw
    coefficients (and \\\phi = \hat\sigma^2\\ for `lm`) converge to the
    full-dataset values.

2.  **vcov check A** (estimator consistency): \\n \cdot \bar V\_{\rm
    analytic} / (n\_{\rm full} \cdot V\_{\rm true}) \to 1\\.

3.  **vcov check B** (calibration): \\n \cdot \widehat{\rm
    Var}(\hat\beta) / (n\_{\rm full} \cdot V\_{\rm true}) \to 1\\.

4.  **Bias(R)**: RESI point estimates converge to the full-dataset
    values.

5.  **sigma2S check A** (estimator consistency): \\\bar{\hat\sigma}^2_S
    / \sigma^2\_{S,\rm true} \to 1\\.

6.  **sigma2S check B** (CI calibration): \\n \cdot \widehat{\rm
    Var}(\hat R) / \sigma^2\_{S,\rm true} \to 1\\.

True values are taken from the full
[`insurance`](https://statimagcoll.github.io/RESI/reference/insurance.md)
dataset using the same definitions as
[`insurancePlasmodeSim`](https://statimagcoll.github.io/RESI/reference/insurancePlasmodeSim.md).
\\\sigma^2_S\\ is extracted from the asymptotic-normal CI half-width:
\\\hat\sigma^2_S = n \cdot ({\rm hw}/z\_{\alpha/2})^2\\.

## Usage

``` r
simCalibrationSim(
  nsim = 500L,
  n.vec = c(100L, 200L, 500L, 1000L, 2000L, 5000L),
  alpha = 0.05,
  output.dir = "resiCalibrationSim",
  fixed.knots = FALSE,
  mc.cores.reps = 1L
)
```

## Arguments

- nsim:

  Integer, replicates per (setting, \\n\\) cell. Default 500.

- n.vec:

  Integer vector of sample sizes. Default
  `c(100, 200, 500, 1000, 2000, 5000)`.

- alpha:

  Numeric, nominal CI level. Default 0.05.

- output.dir:

  Character, directory for raw per-cell RDS files. Default
  `"resiCalibrationSim"`.

- fixed.knots:

  Logical. Fix spline knots at full-dataset quantiles. Default `FALSE`.

- mc.cores.reps:

  Integer, cores for within-cell parallelism. Default 1.

## Value

Invisibly returns the combined metrics `data.frame`.

## See also

[`insurancePlasmodeSim`](https://statimagcoll.github.io/RESI/reference/insurancePlasmodeSim.md),
[`simCompareMethodsFigures`](https://statimagcoll.github.io/RESI/reference/simCompareMethodsFigures.md)
