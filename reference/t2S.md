# Compute the robust effect size index estimate from t statistic (default)

This function computes the robust effect size index from Vandekar, Tao,
& Blume (2020). Vector arguments are accepted. If different length
arguments are passed they are dealt with in the usual way of R.

## Usage

``` r
t2S(t, rdf, n, unbiased = TRUE)
```

## Arguments

- t:

  The t statistic for the parameter of interest.

- rdf:

  Model residual degrees of freedom/degrees of freedom of the t
  statistic.

- n:

  Number of independent samples.

- unbiased:

  Logical, whether to use unbiased or alternative estimator. See
  details.

## Value

Returns a scalar or vector argument of the robust effect size index
estimate.

## Details

This function computes S, the RESI, from a t statistic. The formula for
the unbiased estimator (default) is derived by solving the expected
value of the t statistic for S. It is unbiased and consistent.

The formula for the unbiased conversion is:

\\S = (t \* \sqrt(2) \* \Gamma(rdf/2)) / (\sqrt(n \* rdf) \*
\Gamma((rdf - 1)/2))\\

The formula for the alternative estimator is derived by squaring the t
statistic and using the
[`f2S`](https://statimagcoll.github.io/RESI/reference/f2S.md) formula.
This estimator may be appealing for its intuitive relationship to the F
statistic; the absolute value of RESI estimates using this formula will
be equal to a RESI estimate using an F statistic for the same model.
However, this estimator does have finite sample bias, which is an
important consideration for the coverage of the bootstrapping that
`resi` uses.

The formula for the alternative conversion is:

\\ \sqrt(max(0, (t^2 \* (rdf - 2)/rdf - 1)/rdf))\\

## Examples

``` r
# to obtain t values, first fit a lm
mod = lm(charges ~ region * age + bmi + sex, data = RESI::insurance)
# run lmtest::coeftest to get t values, using a robust variance-covariance formula
ts = lmtest::coeftest(mod, vcov. = sandwich::vcovHC)[,'t value']

# get RESI estimates using unbiased estimator
t2S(ts, n = nrow(RESI::insurance), rdf = mod$df.residual)
#>         (Intercept)     regionnorthwest     regionsoutheast     regionsouthwest 
#>        -0.067297337        -0.026687397        -0.033398602        -0.002467535 
#>                 age                 bmi             sexmale regionnorthwest:age 
#>         0.149720414         0.152301107         0.058360820         0.016662071 
#> regionsoutheast:age regionsouthwest:age 
#>         0.036086004        -0.014963251 

# get RESI estimates using alternative estimator
t2S(ts, n = nrow(RESI::insurance), rdf = mod$df.residual, unbiased = FALSE)
#>         (Intercept)     regionnorthwest     regionsoutheast     regionsouthwest 
#>         -0.06148040          0.00000000         -0.01917451          0.00000000 
#>                 age                 bmi             sexmale regionnorthwest:age 
#>          0.14717461          0.14979819          0.05154917          0.00000000 
#> regionsoutheast:age regionsouthwest:age 
#>          0.02354410          0.00000000 
```
