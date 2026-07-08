# Anova method for resi objects

After running the
[`resi`](https://statimagcoll.github.io/RESI/reference/resi.md) function
on a fitted model, this function can be used to print the Anova-style
table component. If the resi function was run with the \`store.boot =
TRUE\` option to store the full matrix of bootstrapped estimates, the
user can specify a different alpha level for this function's confidence
intervals.

## Usage

``` r
# S3 method for class 'resi'
anova(object, alpha = NULL, ...)
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

Returns an \`anova\` object containing the computed Anova-style table

## Details

The resi function uses the car::Anova function to compute the Anova
table.

## Examples

``` r
# fit a model
mod = lm(charges ~ bmi + sex, data = RESI::insurance)

# run resi with the store.boot = TRUE option
resi.obj = resi(mod, nboot = 100, store.boot = TRUE, alpha = 0.01)

# run anova, specifying a different alpha level if desired
anova(resi.obj, alpha = 0.05)
#>     Df       F   Pr(>F)     RESI    2.5%   97.5%
#> bmi  1 45.2997 0.000000 0.181819 0.13404 0.23396
#> sex  1  3.2439 0.071915 0.040908 0.00000 0.10274
```
