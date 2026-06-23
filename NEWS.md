## RESI 1.3.3 (in development)

### Bug Fixes
* Fixed error when passing a `data` argument containing `NA`s to `resi()` for
  GEE models (`geeglm`). Rows with `NA` in any model variable are now silently
  stripped before bootstrapping, matching the complete cases used to fit the
  model and preventing cluster-size mismatches during re-fitting (#51).
* Fixed misleading error in `summary.resi()` when a different `alpha` level is
  requested but confidence intervals cannot be recomputed. The message now
  clearly explains that either (a) bootstrapping is not supported for the model
  type, or (b) `resi()` was not run with `store.boot = TRUE`, and instructs the
  user to re-run `resi()` with the desired `alpha` directly (#53).
* Fixed `resi()` failing to compute the overall Wald test when the model
  response is a computed expression (e.g. `log10(charges)` or
  `I(charges > 10000)`). Previously the intercept-only reduced model could not
  be fitted inside forked parallel workers because `update()` tried to
  re-evaluate the expression in an environment where the underlying variable was
  not in scope. The fix constructs the reduced-model formula using the
  already-evaluated column name from `model.frame()`, avoiding any re-evaluation.
* Added an informative error when `vcov.args = list(type = "const")` is passed
  together with `vcovfunc = sandwich::vcovHC`. The message explains that
  `type = "const"` is the OLS sandwich (not robust) and directs users to use
  `vcovfunc = stats::vcov` for parametric variance estimation instead (#50).

## RESI 1.3.2

* Added RESI estimation for `emmeans` objects
* Added support for `glmgee` models from `glmtoolbox`
* Added support for Gaussian models from `glmmTMB`
* Added pdf vignette (`vignette("RESI_paper")`)
* Fixed bug with `geeglm` objects not assigning weights properly with missing data
* Documentation fix for Linux systems

## RESI 1.3.0

* Added citation for Journal of Statistical Software paper <doi:10.18637/jss.v112.i03>
* Bug fix in cluster bootstrapping: Now correctly assigns the clustered IDs in the bootstrap
* Minor documentation updates

## RESI 1.2.4

* Minor bug and documentation fixes
* Website deployed

## RESI 1.2.0

* Implemented `boot` package for bootstrapping
  + `boot.results` element of `resi` object contains full `boot` object
  + Allows parallelization
* Updated `geeglm`/`gee` methods
* Revamped `plot.resi` function
* Added `ggplot` methods
* Added `omnibus` function to extract overall Wald test from `resi` object
* Combined `t2S` and `t2S_alt` into `t2S` and `z2S` and `z2S_alt` into `z2S`
* Substantially reduced code duplication between methods
* Expanded error checking and messaging
* Fixed error in `resi` when using a model with only one predictor
* Various typo fixes and cleaned up documentation

## RESI 1.1.1

* Minor argument consistency fix

## RESI 1.1.0

* Expanded longitudinal methods
  + `gee` and `geeglm` reporting valid confidence intervals
  + `gee` and `geeglm` report both a longitudinal RESI estimate and a cross-sectional RESI estimate
  + `lme` and `lmerMod` reports point estimates for longitudinal RESI only (confidence intervals in development)
  + `anova` option added for `geeglm`, `lme`, and `lmerMod`
  + Updated bootstrap sampling to ensure the correct number of clusters
* Error in `f2S` and `t2S_alt` formula corrected (denominator is now n instead of the residual degrees of freedom)
* Added bootstrap fail counter on `resi` for `nls` and `geeglm` model types
* Changed class of `summary` on a `resi` object to `summary_resi` for consistency
* Small typos and inconsistencies fixed

## RESI 1.0.5

* First version
