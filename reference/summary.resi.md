# Summary method for resi objects

After running the
[`resi`](https://statimagcoll.github.io/RESI/reference/resi.md) function
on a fitted model, this function can be used to print the coefficients
table component. If the resi function was run with the \`store.boot =
TRUE\` option to store the full matrix of bootstrapped estimates, the
user can specify a different alpha level for this function's confidence
intervals.

## Usage

``` r
# S3 method for class 'resi'
summary(object, alpha = NULL, ...)
```

## Arguments

- object:

  an object resulting from resi function

- alpha:

  an optional new specification for the confidence level. Can be
  vector-valued

- ...:

  ignored

## Value

Returns a \`summary_resi\` object containing the computed coefficients
table

## Examples

``` r
# fit a model
mod = lm(charges ~ bmi + sex, data = RESI::insurance)

# run resi with the store.boot = TRUE option
resi_obj = resi(mod, nboot = 100, store.boot = TRUE, alpha = 0.01)

# run summary, specifying a different alpha level if desired
summary(resi_obj, alpha = 0.05)
#> 
#> Analysis of effect sizes based on RESI:
#> Confidence level =  0.05
#> Call:  lm(formula = charges ~ bmi + sex, data = RESI::insurance)
#> 
#> Coefficient Table 
#>              Estimate Std. Error t value Pr(>|t|)   RESI    2.5%  97.5%
#> (Intercept)  739.4306  1669.9119  0.4428   0.6580 0.0121 -0.0417 0.0659
#> bmi          389.4347    57.8612  6.7305   0.0000 0.1839  0.1340 0.2340
#> sexmale     1166.9940   647.9397  1.8011   0.0719 0.0492 -0.0043 0.1027
```
