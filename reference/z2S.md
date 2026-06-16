# Compute the robust effect size index estimate from Z statistic

This function computes the robust effect size index from Vandekar, Tao,
& Blume (2020). Vector arguments are accepted. If different length
arguments are passed they are dealt with in the usual way of R.

## Usage

``` r
z2S(z, n, unbiased = TRUE)
```

## Arguments

- z:

  The Z statistic for the parameter of interest.

- n:

  Number of independent samples.

- unbiased:

  Logical, whether to use unbiased or alternative estimator. See
  details.

## Value

Returns a scalar or vector argument of the robust effect size index
estimate.

## Details

This function computes S, the RESI, from a Z statistic. The formula for
the unbiased estimator (default) is derived by solving the expected
value of the Z statistic for S. It is unbiased and consistent.

The formula for the unbiased conversion is:

\\S = Z/\sqrt(n)\\

The formula for the alternative estimator is derived by squaring the Z
statistic and using the
[`chisq2S`](https://statimagcoll.github.io/RESI/reference/chisq2S.md)
formula. This estimator may be appealing for its intuitive relationship
to the Chi-square statistic; the absolute value of RESI estimates using
this formula will be equal to a RESI estimate using a Chi-square
statistic for the same model. However, this estimator does have finite
sample bias, which is an important consideration for the coverage of the
bootstrapping that `resi` uses.

The formula for the alternative conversion is:

\\ \sqrt(max(0, (Z^2 - 1)/n)) \* sign(Z)\\

## Examples

``` r
# to obtain example z values, first fit a glm
mod = glm(charges ~ region * age + bmi + sex, data = RESI::insurance)
# run coeftest to get z values using a robust variance-covariance function
zs = lmtest::coeftest(mod, vcov. = sandwich::vcovHC)[,'z value']

# get RESI estimates using unbiased estimator
z2S(zs, n = nrow(RESI::insurance))
#>         (Intercept)     regionnorthwest     regionsoutheast     regionsouthwest 
#>         -0.06733537         -0.02670248         -0.03341748         -0.00246893 
#>                 age                 bmi             sexmale regionnorthwest:age 
#>          0.14980504          0.15238719          0.05839381          0.01667149 
#> regionsoutheast:age regionsouthwest:age 
#>          0.03610640         -0.01497171 

# get RESI estimates usng alternative estimator
z2S(zs, n = nrow(RESI::insurance), unbiased = FALSE)
#> Warning: NaNs produced
#>         (Intercept)     regionnorthwest     regionsoutheast     regionsouthwest 
#>         -0.06153591          0.00000000         -0.01921832          0.00000000 
#>                 age                 bmi             sexmale regionnorthwest:age 
#>          0.14728939          0.14991488          0.05159896          0.00000000 
#> regionsoutheast:age regionsouthwest:age 
#>          0.02358576          0.00000000 
```
