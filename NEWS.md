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
