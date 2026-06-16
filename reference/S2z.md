# Convert RESI (S) estimate to Z statistic

Converts the robust effect size index (S) to Z statistic. Vector
arguments are accepted. If different length arguments are passed they
are dealt with in the usual way of R.

## Usage

``` r
S2z(S, n, unbiased = TRUE)
```

## Arguments

- S:

  The value of the RESI estimate.

- n:

  Number of independent samples.

- unbiased:

  Logical, whether the unbiased or alternative estimator was used to
  compute RESI estimate. Default is TRUE.

## Value

Returns a scalar or vector argument of the Chi-square statistic.

## Details

The formula for converting a RESI estimate to a corresponding Z
statistic depends on which estimator was used to compute the RESI
estimate (unbiased vs. alternative, see
[`z2S`](https://statimagcoll.github.io/RESI/reference/z2S.md)). For the
unbiased estimator, the RESI can be positive or negative and there is a
1-1 transformation from S to Z. The formula for converting S (unbiased)
to the Z statistic is:

\\\sqrt(n)\*S\\

For the alternative formula, if the RESI estimate is 0, the Z statistic
is only known within an interval, \[-1, 1\]. For a non-zero S, the
formula is:

\\\sqrt{S^2}/S\sqrt(n\*abs(S) + 1)\\

## Examples

``` r
# convert S estimates with corresponding degrees of freedom to
# Z statistics estimates (using unbiased formula)
S_ests = c(-0.2, 0, 0.1)
S2z(S = S_ests, n = 300, unbiased = TRUE)
#> [1] -3.464102  0.000000  1.732051

# convert S estimates with corresponding degrees of freedom to
# Z statistics estimates (using alernative formula)
S_ests = c(-0.2, 0, 0.1)
S2z(S = S_ests, n = 300, unbiased = FALSE)
#> Warning: Function is not 1-1 for S = 0, Z statistic is between -1 and 1
#> [1] -7.810250        NA  5.567764
```
