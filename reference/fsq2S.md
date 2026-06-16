# Covert Cohen's *f*^2 to S

Converts Cohen's *f*^2 to robust effect size index (S) using the formula
from Vandekar, Tao, & Blume (2020).

## Usage

``` r
fsq2S(fsq)
```

## Arguments

- fsq:

  Numeric, value of Cohen's *f*^2.

## Value

Returns an estimate the robust effect size index

## Details

The formula for the conversion is:

\\S = \sqrt(f^2)\\

## Examples

``` r
# consider a moderate effect size of f^2 = 0.3
fsq2S(0.3)
#> [1] 0.5477226
# This corresponds to a RESI of 0.5477226
```
