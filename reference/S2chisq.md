# Convert non-zero S to Chi-square statistic

Converts the robust effect size index (S) to Chi-square statistic, given
that S is greater than 0. For an S value of 0, only an upper bound on
the Chi-square statistic can be computed. Vector arguments are accepted.
If different length arguments are passed they are dealt with in the
usual way of R.

## Usage

``` r
S2chisq(S, df, n)
```

## Arguments

- S:

  The value of the RESI estimate.

- df:

  Number of degrees of freedom of the chi-square statistic.

- n:

  Number of independent samples.

## Value

Returns a scalar or vector argument of the Chi-square statistic.

## Details

The formula for converting a RESI estimate above 0 to Chi-square
statistic is:

\\ chisq = n\*S^2 + df\\

If the RESI estimate is 0, all that is known is that the Chi-square
statistic is less than or equal to the degrees of freedom.

## Examples

``` r
# convert S estimates with corresponding degrees of freedom to Chi-square estimates
S_ests = c(0.2, 0.4, 0.6)
dfs = c(2, 1, 3)
S2chisq(S = S_ests, df = dfs, n = 300)
#> [1]  14  49 111
```
