# Covert Cohen's *d* to \|S\|

Converts Cohen's *d* robust effect size index (S) using the formula from
Vandekar, Tao, & Blume (2020).

## Usage

``` r
d2S(d, pi = 0.5)
```

## Arguments

- d:

  Numeric, value of Cohen's *d*.

- pi:

  Numeric, the sampling proportions.

## Value

Returns an estimate the robust effect size index

## Details

The pi parameter comes from the fact that Cohen's d doesn't account for
unequal sample proportions in the population, but S does.

The default is set to a natural value 1/2, which corresponds to a case
control design, for example, where sampling proportions always are
controlled by the experimenter.

The formula to convert Cohen's *d* to S is:

\\S = d/\sqrt( 1/\pi + 1/ (1 - \pi))\\

## Examples

``` r
# Consider an experiment with equal sampling proportions and a medium effect size
# corresponding to a Cohen's d of 0.5.
# convert to RESI (S)
d2S(d = 0.5)
#> [1] 0.25

# This corresponds to a RESI of 0.25.
```
