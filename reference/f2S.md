# Compute the robust effect size index estimate from F-statistic

This function computes the robust effect size index from Vandekar, Tao,
& Blume (2020). Vector arguments are accepted. If different length
arguments are passed they are dealt with in the usual way of R.

## Usage

``` r
f2S(f, df, rdf, n)
```

## Arguments

- f:

  The F statistic for the parameter of interest.

- df:

  Number of degrees of freedom of the F statistic.

- rdf:

  Model residual degrees of freedom.

- n:

  Number of independent samples.

## Value

Returns a scalar or vector argument of the robust effect size index
estimate.

## Details

The formula for converting an F statistic to S is:

\\ S = \sqrt(max(0, (f \* df \* (rdf - 2)/rdf - df)/n))\\

The estimator is derived by setting the statistic equal to the expected
value of the test statistic and solving for S.

## Examples

``` r

# to obtain example F values, first fit a lm
mod = lm(charges ~ region * age + bmi + sex, data = RESI::insurance)

# run Anova, using a robust variance-covariance function
# get the F values and Df values
fs = car::Anova(mod, vcov. = sandwich::vcovHC)[1:5, "F"]
#> Coefficient covariances computed by sandwich::vcovHC
dfs = car::Anova(mod, vcov. = sandwich::vcovHC)[1:5, "Df"]
#> Coefficient covariances computed by sandwich::vcovHC

# get RESI estimates
f2S(fs, df = dfs, rdf = mod$df.residual, n = nrow(RESI::insurance))
#> [1] 0.03647985 0.29511133 0.14979819 0.05154917 0.01605603
```
