# Covert R^2 to S

Converts R^2, the partial coefficient of determination, to robust effect
size index (S) using the formula from Vandekar, Tao, & Blume (2020).

## Usage

``` r
Rsq2S(Rsq)
```

## Arguments

- Rsq:

  Numeric, R^2

## Value

Returns an estimate of R^2 based on the RESI

## Details

The formula for the conversion is:

\\S = \sqrt((-R^2)/(R^2 - 1))\\

## Examples

``` r
# consider a moderate effect size of R^2 = 0.1
Rsq2S(0.1)
#> [1] 0.3333333
# this corresponds to a RESI of 0.333
```
