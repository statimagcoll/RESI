# Omnibus (Overall) Wald Test for resi objects

After running the
[`resi`](https://statimagcoll.github.io/RESI/reference/resi.md) function
on a fitted model, this function can be used to print the overall Wald
test component. If the resi function was run with the \`store.boot =
TRUE\` option to store the full matrix of bootstrapped estimates, the
user can specify a different alpha level for this function's confidence
intervals.

## Usage

``` r
omnibus(object, alpha = NULL, ...)
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

Returns a \`omnibus_resi\` object containing the computed omnibus Wald
test

## Examples

``` r
# fit a model
mod = lm(charges ~ bmi + sex, data = RESI::insurance)

# run resi with the store.boot = TRUE option (ci.method must be "boot")
resi_obj = resi(mod, nboot = 100, store.boot = TRUE, alpha = 0.01, ci.method = "boot")

# run summary, specifying a different alpha level if desired
omnibus(resi_obj, alpha = 0.05)
#> 
#> Analysis of effect sizes based on RESI:
#> Confidence level =  0.05 
#>   Res.Df Df      F Pr(>F)   RESI   2.5%  97.5%
#> 1   1337                                      
#> 2   1335  2 24.338      0 0.1866 0.1452 0.2359
```
