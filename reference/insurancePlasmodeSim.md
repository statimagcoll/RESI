# Insurance Plasmode Simulation for RESI Evaluation

Runs a plasmode simulation study using the
[`insurance`](https://statimagcoll.github.io/RESI/reference/insurance.md)
dataset to evaluate RESI confidence interval performance. In each
replicate, `n` observations are resampled with replacement from the full
insurance dataset (*N* = 1338). The RESI point estimates from
[`resi_pe`](https://statimagcoll.github.io/RESI/reference/resi_pe.md)
applied to the full dataset are treated as the true parameter values for
computing bias, MSE, CI coverage, and CI width.

## Usage

``` r
insurancePlasmodeSim(
  nsim = 1000L,
  n.vec = c(50, 100, 200, 500, 1000, 2000, 5000),
  nboot = 500L,
  alpha = 0.05,
  ci.method = c("boot", "normal", "qf", "cf"),
  output.dir = NULL,
  fixed.knots = FALSE,
  mc.cores.settings = 1L,
  mc.cores.reps = 1L
)
```

## Arguments

- nsim:

  Integer, number of simulation replicates per (setting, `n`) cell.
  Default 1000. Use 10 for initial testing.

- n.vec:

  Integer vector of sample sizes. Default
  `c(50, 100, 200, 500, 1000, 2000, 5000)`.

- nboot:

  Integer, bootstrap replicates per internal
  [`resi`](https://statimagcoll.github.io/RESI/reference/resi.md) call.
  Default 500. Use 10 for initial testing. Ignored when
  `ci.method != "boot"`.

- alpha:

  Numeric, CI significance level. Default 0.05.

- ci.method:

  Character, CI method passed to
  [`resi`](https://statimagcoll.github.io/RESI/reference/resi.md). One
  of `"boot"` (bootstrap, default), `"normal"` (asymptotic
  truncated-normal), or `"qf"` (asymptotic quadratic-form / Imhof). When
  `ci.method != "boot"`, `nboot` is ignored.

- output.dir:

  Character, path to the directory where all results are saved. Created
  if it does not exist. Defaults to `"resiBootSim"`,
  `"resiAsympNormalSim"`, or `"resiAsympQFSim"` based on `ci.method`
  when `NULL`.

- fixed.knots:

  Logical. If `TRUE`, spline knots are fixed at the empirical tertiles
  of `age` in the full insurance dataset rather than re-selected by
  `df = 3` in each bootstrap sample. Default `FALSE`.

- mc.cores.settings:

  Integer, cores for the outer `mclapply` over (setting \\\times\\
  sample size) combinations. Default 1.

- mc.cores.reps:

  Integer, cores for the inner `mclapply` over simulation replicates
  within each (setting, `n`) cell. Default 1.

## Value

Invisibly returns the summary metrics `data.frame`. Side effects:

- `output.dir/sim_raw/<setting>_n<n>.rds`: list of per-replicate `anova`
  and `coefficients` tables.

- `output.dir/summary_table.rds`: combined metrics table with columns
  `model, vcov, n, n_success, table, term, bias, mse, coverage, width`.

## Details

Two models are evaluated:

- **lm**: `log10(charges) ~ ns(age, df=3) * sex + bmi + smoker + region`

- **glm**:
  `I(charges > 10000) ~ ns(age, df=3) * sex + bmi + smoker + region`
  with `family = binomial()`

Each model is evaluated under both parametric (`vcovfunc = stats::vcov`)
and robust (`vcovfunc = sandwich::vcovHC`) variance settings, yielding
four simulation conditions.

Parallelization is via
[`mclapply`](https://rdrr.io/r/parallel/mclapply.html), which uses
forking and is not supported on Windows (falls back to sequential
evaluation on Windows).

## See also

[`simFigures`](https://statimagcoll.github.io/RESI/reference/simFigures.md),
[`simCompareMethodsFigures`](https://statimagcoll.github.io/RESI/reference/simCompareMethodsFigures.md),
[`resi`](https://statimagcoll.github.io/RESI/reference/resi.md),
[`resi_pe`](https://statimagcoll.github.io/RESI/reference/resi_pe.md)
