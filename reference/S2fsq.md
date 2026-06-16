# Covert S to Cohen's *f*^2

Converts robust effect size index (S) to Cohen's *f*^2 (effect size for
multiple regression) using the formula from Vandekar, Tao, & Blume
(2020).

## Usage

``` r
S2fsq(S)
```

## Arguments

- S:

  Numeric,the robust effect size index.

## Value

Returns an estimate of Cohen's *f*^2 based on the RESI

## Details

The formula for the conversion is:

\\ f^2 = S^2\\

## Examples

``` r
# fit a linear regression model with continuous outcome and predictor
mod = lm(charges ~ age, data = RESI::insurance)

# obtain t value for calculating RESI
t = summary(mod)$coefficients[2, "t value"]

# calculate RESI
S = t2S(t, n = 1338, rdf = 1336)

# convert to f^2
S2fsq(S)
#> [1] 0.09792731
```
