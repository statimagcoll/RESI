# Recompute Summary Table from Raw Simulation Output

Reads the per-replicate raw `.rds` files saved by
[`insurancePlasmodeSim`](https://statimagcoll.github.io/RESI/reference/insurancePlasmodeSim.md)
and recomputes the summary metrics table (bias, MSE, coverage,
upper_coverage, lower_coverage, width). Use this whenever
`.simComputeMetrics` has been updated (e.g. new metrics added) without
needing to rerun the full simulation.

## Usage

``` r
simRecomputeSummary(
  output.dir = "resiBootSim",
  alpha = 0.05,
  fixed.knots = FALSE
)
```

## Arguments

- output.dir:

  Character, directory containing simulation output. Default
  `"resiBootSim"`.

- alpha:

  Numeric, CI significance level used in the original simulation.
  Default 0.05.

- fixed.knots:

  Logical. Must match the value used in the original
  [`insurancePlasmodeSim`](https://statimagcoll.github.io/RESI/reference/insurancePlasmodeSim.md)
  call so that the true RESI values are computed from the same model
  formula. Default `FALSE`.

## Value

Invisibly returns the updated summary `data.frame`. Overwrites
`output.dir/summary_table.rds`.

## Details

True RESI values are re-estimated from the full
[`insurance`](https://statimagcoll.github.io/RESI/reference/insurance.md)
dataset using the same formulae and variance functions as the original
simulation.

## See also

[`insurancePlasmodeSim`](https://statimagcoll.github.io/RESI/reference/insurancePlasmodeSim.md),
[`simFigures`](https://statimagcoll.github.io/RESI/reference/simFigures.md)
