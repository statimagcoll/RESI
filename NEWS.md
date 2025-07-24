## RESI 1.3.1

* 

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
