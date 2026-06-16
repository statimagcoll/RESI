# Plotting RESI Estimates and CIs

This function uses base graphics to plot robust effect size (RESI)
estimates and confidence intervals from \`resi\`, \`summary_resi\`, and
\`anova_resi\` objects.

## Usage

``` r
# S3 method for class 'resi'
plot(
  x,
  alpha = NULL,
  ycex.axis = NULL,
  yaxis.args = list(),
  automar = TRUE,
  ...
)
```

## Arguments

- x:

  Object of \`resi\`, \`summary_resi\`, or \`anova_resi\` class

- alpha:

  Numeric, desired alpha level for confidence intervals

- ycex.axis:

  Numeric, scale specifically for the variable name labels

- yaxis.args:

  List, other arguments to be passed to
  [`axis`](https://rdrr.io/r/graphics/axis.html) for the y-axis

- automar:

  Logical, whether to automatically adjust the plotting margins to
  accommodate variable names. Default = \`TRUE\`

- ...:

  Other graphical parameters passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) and
  [`lines`](https://rdrr.io/r/graphics/lines.html)

## Value

Returns a plot of RESI point estimates

## Details

This function creates a forest-like plot with RESI estimates for each
variable or factor. The size of the left margin will be automatically
adjusted (and returned to original after plotting) unless \`automar =
FALSE\`. Additional graphics parameters will be passed to the main plot
function, the confidence intervals. Arguments specifically for the
y-axis (variable names) can be specified using \`yaxis.args\`. To
manually adjust the size of the y-axis labels without affecting the
x-axis, the user can specify a value for \`ycex.axis\`.

## Examples

``` r
# create a resi object
resi_obj <- resi(lm(charges ~ region * age + bmi + sex, data = RESI::insurance),
nboot = 10)

# plot coefficients table, changing size of labels for both axes in the usual way
plot(resi_obj, cex.axis = 0.7)


# plot ANOVA table, changing the size of just the y-axis
plot(resi_obj, ycex.axis = 0.8)
```
