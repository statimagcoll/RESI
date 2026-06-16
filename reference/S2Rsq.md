# Covert S to R^2

Converts robust effect size index (S) to R^2, the partial coefficient of
determination, using the formula from Vandekar, Tao, & Blume (2020).

## Usage

``` r
S2Rsq(S)
```

## Arguments

- S:

  Numeric, the robust effect size index.

## Value

Returns an estimate of R^2 based on the RESI

## Details

The formula for the conversion is:

\\ R^2 = S^2 / (1 + S^2)\\

## Examples

``` r
# fit a simple linear regression with a binary predictor
mod = lm(charges ~ sex, data = RESI::insurance)

# calculate t-value
t = summary(mod)$coefficients[2, "t value"]

# calculate RESI (S)
S = t2S(t, n = 1338, rdf = 1336)

# convert S to R^2
S2Rsq(S)
#> [1] 0.003273823
```
