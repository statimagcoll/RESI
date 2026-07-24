# Asymptotic CI Validation: Simulation Results

## Overview

This article presents simulation-based validation of the asymptotic
confidence interval (CI) methods available in the RESI package via
[`resi_pe_asymptotic()`](https://statimagcoll.github.io/RESI/reference/resi_pe_asymptotic.md)
and [`resi()`](https://statimagcoll.github.io/RESI/reference/resi.md).

All asymptotic methods share the same $`\hat\Sigma_R`$ estimator — the
**extended influence-function** approach that treats the sandwich bread
$`\hat\Sigma_{XA}`$ and meat $`\hat\Sigma_{XB}`$ as independent
estimators (argument `deriv_method = "extended"` internally) — and
differ only in how they invert the asymptotic distribution of the Wald
statistic $`T^2 = n\tilde{S}^2`$.

| `ci.method` | Description |
|----|----|
| `"qf"` **(default)** | Quadratic-form CI: inverts the exact weighted chi-square distribution via the Imhof/Davies algorithm |
| `"normal"` | Truncated-normal CI using a bias-corrected normal approximation |
| `"cf"` | Cornish-Fisher CI: inverts a Cornish-Fisher expansion of the weighted chi-square |

Results are shown for four estimation settings:

| Setting | Model | Variance estimator |
|----|----|----|
| lm parametric | [`lm()`](https://rdrr.io/r/stats/lm.html) | [`stats::vcov`](https://rdrr.io/r/stats/vcov.html) |
| lm robust | [`lm()`](https://rdrr.io/r/stats/lm.html) | [`sandwich::vcovHC`](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html) (HC3) |
| glm parametric | `glm(..., family=binomial)` | [`stats::vcov`](https://rdrr.io/r/stats/vcov.html) |
| glm robust | `glm(..., family=binomial)` | [`sandwich::vcovHC`](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html) (HC3) |

For background on the RESI, see the [RESI package
website](https://statimagcoll.github.io/RESI/) and [Jones et
al. (2025)](https://doi.org/10.18637/jss.v112.i03).

------------------------------------------------------------------------

## Simulation design

### Plasmode simulation

All evaluations use a **plasmode simulation design** ([Gadbury et
al. 2008](https://doi.org/10.1177/0962280207080302)) built around the
`insurance` dataset ($`N = 1338`$) shipped with the package. Bootstrap
resampling (with replacement) from the full dataset gives the empirical
sampling distribution. RESI values estimated on the full dataset via
[`resi_pe()`](https://statimagcoll.github.io/RESI/reference/resi_pe.md)
serve as the population targets.

Two models are evaluated:

- **Linear model** (`lm`):
  `log10(charges) ~ ns(age, df=3) * sex + bmi + smoker + region`
- **Logistic model** (`glm`, binomial):
  `I(charges > 15000) ~ ns(age, df=3) * sex + bmi + smoker + region`

Each simulation cell (model × variance estimator × sample size) uses **1
000 replicates**. Results are shown for sample sizes
$`n \in \{50, 100, 200, 500, 1000, 2000, 5000\}`$.

### Metrics

For each model coefficient and ANOVA term the following quantities are
reported:

- **Coverage**: proportion of replicates in which the 95% CI contains
  the population RESI value.
- **Width**: median CI width.
- **Bias** and **MSE**: mean signed and mean squared deviation of
  $`\hat{S}`$ from the population value.

Bootstrap CIs (1 000 replicates, non-parametric) provide the reference
standard.

### Variance calibration

The calibration study examines whether the analytic $`\hat\sigma^2_R`$
(the CI width estimator) is consistent for the true Monte Carlo variance
of $`\hat{R}`$ across replications. Each calibration figure has four
panels:

1.  **Bias of $`\hat\theta`$**: signed deviation of the parameter
    estimates (model coefficients and, for `lm`,
    $`\hat\phi = \hat\sigma^2`$) from their full-data targets.
2.  **Sandwich covariance normalized bias**:
    $`(n\,\widehat{\mathrm{Var}}(\hat\beta_k) -
    n\,\mathrm{Var}_{\mathrm{MC}}(\hat\beta_k)) / n\,\mathrm{Var}_{\mathrm{MC}}(\hat\beta_k)`$.
3.  **Bias of $`\hat{R}`$** (RESI point estimate).
4.  **$`\sigma^2_R`$ ratios**: solid line = $`\sigma^2_{R,\text{MC}} /
    \sigma^2_{R,\text{extended}}`$ (Monte Carlo variance over extended
    estimator); dashed = ratio against the initial delta-method
    derivation. A ratio of 1 (dashed reference line) indicates perfect
    calibration.

------------------------------------------------------------------------

## How to reproduce these figures

The chunks below are **not evaluated** when the article is rendered. Run
them from a working directory of your choice; then copy the resulting
PDFs into `vignettes/articles/figures/` inside the package source tree.

### Step 1 — run plasmode simulations

``` r

library(RESI)

# Bootstrap reference  (~hours; adjust mc.cores to your hardware)
insurancePlasmodeSim(nsim = 1000, ci.method = "boot",
                     output.dir = "resiBootSim",
                     mc.cores.settings = 4)

# Normal approximation
insurancePlasmodeSim(nsim = 1000, ci.method = "normal",
                     output.dir = "resiAsympNormalSim",
                     mc.cores.settings = 4)

# Quadratic-form (default)
insurancePlasmodeSim(nsim = 1000, ci.method = "qf",
                     output.dir = "resiAsympQFSim",
                     mc.cores.settings = 4)

# Cornish-Fisher
insurancePlasmodeSim(nsim = 1000, ci.method = "cf",
                     output.dir = "resiAsympCFSim",
                     mc.cores.settings = 4)
```

### Step 2 — run variance calibration simulation

``` r

# 500 replicates per (model × vcov × n) cell; extended estimator
simCalibrationSim(nsim = 500,
                  output.dir   = "resiCalibrationSim",
                  deriv_method = "extended",
                  mc.cores.reps = 4)
```

### Step 3 — generate figures

``` r

# QF detail figures (one PDF per model/vcov setting)
simFigures(output.dir = "resiAsympQFSim")

# Variance calibration figures
simCalibrationFigures(output.dir  = "resiCalibrationSim",
                      deriv_method = "extended")

# CI method comparison, coverage quantile, and estimator bias/MSE figures
simCompareMethodsFigures(
  output.dirs = c(Bootstrap = "resiBootSim",
                  Normal    = "resiAsympNormalSim",
                  QF        = "resiAsympQFSim",
                  CF        = "resiAsympCFSim"),
  figures.dir = "method_comparison"
)

# t2S/f2S vs z2S/chisq2S estimator comparison (lm only)
simEstimatorFigures(sim.dir     = "resiBootSim",
                    figures.dir = "method_comparison")
```

### Step 4 — copy figures into the package

``` r

# Adjust pkg_dir to the root of your RESI package source tree
pkg_dir  <- "path/to/RESI"
figs_dir <- file.path(pkg_dir, "vignettes", "articles", "figures")
dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)

file.copy(list.files("resiCalibrationSim/figures/extended",
                     pattern = "\\.pdf$", full.names = TRUE), figs_dir)
file.copy(list.files("resiAsympQFSim/figures",
                     pattern = "\\.pdf$", full.names = TRUE), figs_dir)
file.copy(list.files("method_comparison",
                     pattern = "\\.pdf$", full.names = TRUE), figs_dir)
```

------------------------------------------------------------------------

## Variance estimator calibration

These figures show whether $`\hat\sigma^2_R`$ — the width estimator
shared by all three asymptotic CIs — is well-calibrated as a function of
sample size.

### LM parametric

### LM robust

### GLM parametric

### GLM robust

------------------------------------------------------------------------

## CI method comparison

Coverage (top) and width (bottom) across all four CI methods —
Bootstrap, Normal, QF, and Cornish-Fisher — as a function of sample
size. The nominal coverage is 95% (dashed reference line). Each panel
shows results for one regression term; terms are labelled with their
population RESI value.

### LM parametric — coefficients

### LM parametric — ANOVA

### LM robust — coefficients

### LM robust — ANOVA

### GLM parametric — coefficients

### GLM parametric — ANOVA

### GLM robust — coefficients

### GLM robust — ANOVA

------------------------------------------------------------------------

## Coverage quantile comparison

These figures show the distribution of per-term coverage across the
sample-size grid, displayed as quantile bands. Unlike section 4, which
shows mean coverage per term, this view reveals how variable coverage is
across terms and highlights any terms where a method systematically
under- or over-covers.

### LM parametric — coefficients

### LM parametric — ANOVA

### LM robust — coefficients

### LM robust — ANOVA

### GLM parametric — coefficients

### GLM parametric — ANOVA

### GLM robust — coefficients

### GLM robust — ANOVA

------------------------------------------------------------------------

## Point estimate bias and MSE

These figures show how the RESI **point estimate** behaves across CI
methods (bias and MSE per term) as a function of sample size. Since all
methods use the same plug-in estimator $`\hat{S}`$, bias and MSE should
be identical across methods; any divergence reflects numerical
differences in small samples.

### LM parametric — coefficients

### LM parametric — ANOVA

### LM robust — coefficients

### LM robust — ANOVA

### GLM parametric — coefficients

### GLM parametric — ANOVA

### GLM robust — coefficients

### GLM robust — ANOVA

------------------------------------------------------------------------

## t/F vs z/$`\chi^2`$ estimator comparison (lm)

The `t2S`/`f2S` functions convert test statistics using a **t- or
F-distribution** quantile (accounting for residual degrees of freedom);
`z2S`/`chisq2S` use the **normal or chi-square** quantile. These figures
show how the two estimator families compare in bias and MSE for the `lm`
model, based on bootstrap simulation replicates.

### LM parametric — coefficients

### LM parametric — ANOVA

### LM robust — coefficients

### LM robust — ANOVA

------------------------------------------------------------------------

## QF method detail

Per-term coverage and width for the quadratic-form (QF) CI — the package
default — across the full sample-size grid.

### LM parametric

### LM robust

### GLM parametric

### GLM robust

------------------------------------------------------------------------

## Session info

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
    ##  [5] xfun_0.60         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
    ##  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
    ## [13] pkgdown_2.2.1     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
    ## [17] compiler_4.6.1    tools_4.6.1       ragg_1.5.2        bslib_0.11.0     
    ## [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
    ## [25] rlang_1.3.0       fs_2.1.0
