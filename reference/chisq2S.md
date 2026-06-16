# Compute the robust effect size index estimate from chi-squared statistic.

This function computes the robust effect size index from Vandekar, Tao,
& Blume (2020). Vector arguments are accepted. If different length
arguments are passed they are dealt with in the usual way of R. For
mixed effects models, RESI is conditional on the average correlation
structure within subjects.

## Usage

``` r
chisq2S(chisq, df, n)
```

## Arguments

- chisq:

  The chi-square statistic for the parameter of interest.

- df:

  Number of degrees of freedom of the chi-square statistic.

- n:

  Number of independent samples.

## Value

Returns a scalar or vector argument of the robust effect size index
estimate.

## Details

The formula for converting a Chi-square statistic to RESI is:

\\ S = \sqrt(max( 0, (chisq - df)/n))\\

## Examples

``` r
# obtain Chi-sq value by fitting an lm and running a Wald test
mod = lm(charges ~ region * age + bmi + sex, data = RESI::insurance)

# run a Wald test with robust variance
wt = lmtest::waldtest(mod, vcov = sandwich::vcovHC, test = "Chisq")

# get Chi-sq value and degrees of freedom
chisq = wt$Chisq[2]
df = abs(wt$Df[2])

# run chisq2S to convert to RESI
chisq2S(chisq, df = df, n = nrow(mod$model))
#> [1] 0.3598259
```
