# RESI

RESI is an R package designed to implement the Robust Effect Size Index
(RESI, denoted as S) described in Vandekar, Tao, & Blume (2020). The
RESI is a versatile effect size measure that can be easily computed and
added to common reports (such as summary and ANOVA tables). This package
currently supports `lm`, `glm`, `nls`, `survreg`, `coxph`, `hurdle`,
`zeroinfl`, `gee`, `geeglm`, `lme`, and `lmerMod` models. Nonparametric
bootstrapping is used to compute confidence intervals, although the
interval performance has not yet been evaluated for the longitudinal
models. A Bayesian bootstrap is also available for `lm` and `nls`
models. In addition to the main `resi` function, the package also
includes a point-estimate-only function (`resi_pe`), conversions from S
to other common effect size measures and vice versa, print methods, plot
methods, summary methods, and Anova/anova methods. A more detailed
vignette is being written.

If you would like to contribute to the package, please branch off of our
[GitHub](https://github.com/statimagcoll/RESI) and submit a pull request
describing the contribution. Please use the [GitHub
Issues](https://github.com/statimagcoll/RESI/issues) page to report any
problems and the
[Discussions](https://github.com/statimagcoll/RESI/discussions/categories/q-a)
page to seek additional support.

## References

Jones M, Kang K, Vandekar S (2025). RESI: An R Package for Robust Effect
Sizes. *Journal of Statistical Software, 112(3), 1–27.
*<https://doi.org/10.18637/jss.v112.i03>.**

Kang, K., Jones, M. T., Armstrong, K., Avery, S., McHugo, M., Heckers,
S., & Vandekar, S. Accurate Confidence and Bayesian Interval Estimation
for Non-centrality Parameters and Effect Size Indices. *Psychometrika.*
2023. 10.1007/s11336-022-09899-x. Advance online publication.
<https://doi.org/10.1007/s11336-022-09899-x>.

Vandekar S, Tao R, Blume J. A Robust Effect Size Index. *Psychometrika.*
2020 Mar;85(1):232-246. doi: 10.1007/s11336-020-09698-2.
