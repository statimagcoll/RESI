# Plotting RESI Estimates and CIs

This function uses ggplot2 graphics to plot robust effect size (RESI)
estimates and confidence intervals from \`resi\`, \`summary_resi\`, and
\`anova_resi\` objects.

## Usage

``` r
# S3 method for class 'resi'
ggplot(data, mapping, alpha = NULL, error.bars = TRUE, ..., environment)
```

## Arguments

- data:

  Object of \`resi\`, \`summary_resi\`, or \`anova_resi\` class

- mapping:

  Ignored, included for consistency with \`ggplot\` generic

- alpha:

  Numeric, desired alpha level for confidence intervals

- error.bars:

  Logical, whether to include end caps on the confidence intervals.
  Default = \`TRUE\`

- ...:

  Ignored

- environment:

  Ignored, included for consistency with \`ggplot\` generic

## Value

Returns a ggplot of RESI point estimates

## Examples

``` r
# create a resi object
resi_obj <- resi(lm(charges ~ region * age + bmi + sex, data = RESI::insurance),
nboot = 10)

# plot ANOVA table
ggplot2::ggplot(anova(resi_obj))
```
