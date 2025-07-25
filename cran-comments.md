## Test Environments
* local R installation, version 4.4.1
* rhub::rhub_check() (linux, macos, windows, atlas)
* win-builder

## R CMD check results
There were no ERRORs or WARNINGs.
There was one NOTE, a false positive for resi_pe.geeglm: no visible binding for global variable 'w'. w is used for the weights argument of the glm function, which expects a column name (without quotes).

## CRAN Submission 1.3.1
This update fixes a few outstanding ISSUEs from our GitHub and adds support for a couple of model classes, namely `glmgee`, `glmmTMB`, and `emmGrid`. It also provides the published JSS manuscript as a pdf vignette.

Words flagged by spell_check are not misspelled.

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Previous CRAN Comments

## Test Environments
* local R installation, version 4.4.1
* rhub::check() (linux, macos, windows)
* win-builder

## R CMD check results
There were no ERRORs or WARNINGs.
There was one NOTE with win-builder related to the DOI. The DOI in the CITATION is for a new JSS publication that will be registered after publication on CRAN. 

## CRAN Submission 1.3.0
This is a requested update by the editor of the *Journal of Statistical Software* to include citations to the forthcoming publication regarding the package. The DOI in the CITATION is for a new JSS publication that will be registered after publication on CRAN. There is also a minor bug fix included for assigning the bootstrap ids in the cluster bootstrap.

Words flagged by spell_check are not misspelled.

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Test Environments
* local R installation, version 4.3.1
* rhub::check_for_cran()
* win-builder

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## CRAN Submission 1.2.4
This is the same as submission 1.2.1/1.2.2/1.2.3 which did not pass incoming checks due to installation CPU time surpassing elapsed time. I removed the `RESIsetup.R` file because it did not solve the issue. I copied the one function we use from the `regtools` package to `internal_functions.R` and moved `regtools` to the Suggests field rather than imports, as `regtools` relies on `data.table` and may use more than 2 threads.

Words flagged by spell_check are not misspelled.

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Test Environments
* local R installation, version 4.3.1
* rhub::check_for_cran()
* win-builder

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## CRAN Submission 1.2.3
This is the same as submission 1.2.1/1.2.2 which did not pass incoming checks due to CPU time surpassing elapsed time. I have attempted to fix this by changing the function  in the `RESIsetup.R` file to `Sys.setenv("OMP_THREAD_LIMIT" = 1)` and added this line to the `RESI.Rmd` vignette as well.

Words flagged by spell_check are not misspelled.

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Test Environments
* local R installation, version 4.3.1
* rhub::check_for_cran()
* win-builder

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## CRAN Submission 1.2.2
This is the same as submission 1.2.1 which did not pass incoming checks due to CPU time surpassing elapsed time. This has been addressed by adding a function in the `RESIsetup.R` file to set `Sys.setenv("OMP_THREAD_LIMIT" = 2)`.

Words flagged by spell_check are not misspelled.

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Test Environments
* local R installation, version 4.3.1
* rhub::check_for_cran()
* win-builder

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## CRAN Submission 1.2.1
This update makes minor bug fixes and documentation fixes in response to reviewer comments received after submitting a manuscript regarding the package to the *Journal of Statistical Software.* The documentation for `f2S` and `plot.resi` have been edited slightly. The README has been updated with details for how to contribute to the package or report errors. A pkgdown website has been deployed.

Words flagged by spell_check are not misspelled.

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Test Environments
* local R installation, version 4.2.2
* rhub::check_for_cran()
* win-builder

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## CRAN Submission 1.2.0
This update reworks the package to address initial comments received after submitting a manuscript regarding the package to the *Journal of Statistical Software.* The main functions, `resi` and `resi_pe` have been reworked to minimize code duplication. The bootstrapping is now handled via the `boot` package, with an added internal function `resi_stat` to compute the estimates. `t2S_alt` and `t2S` were collapsed into one function (similarly for `z2S` and `z2S_alt`). An `omnibus` function was added to extract the overall Wald test from a `resi` object.  The plotting methods have been improved and `ggplot` methods have been added. The methods for `gee` and `geeglm` were updated based on ongoing theoretical work. A bug that was causing an error in `resi` on models containing only one predictor was resolved.

Words flagged by spell_check are not misspelled.

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Test Environments
* local R installation, version 4.2.2
* rhub::check_for_cran()
* win-builder

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## CRAN Submission 1.1.1
This update changes "object" to "mod" in Anova.resi for S3 method consistency with car::Anova, based on warning that was was given by CRAN's package check and makes a minor change to resi.geeglm.

Words flagged by spell_check are not misspelled.

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Test Environments
* local R installation, version 4.2.2
* rhub::check_for_cran()
* win-builder


## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## CRAN Submission 1.1.0
This is an update that expands the RESI estimation for longitudinal models and corrects several minor errors and bugs. There was a slight formula change in the f2S and t2S_alt function. The regtools package has been moved from Imports to Suggests for the moment due to receiving a notice of a downstream dependency issue regarding the float package. The one function used from that package is currently copied into the internal_functions.R here with comments in the code.

Words flagged by spell_check are not misspelled.

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Test Environments
* local R installation, version 4.1.2
* rhub::check_for_cran()
* win-builder


## R CMD check results
There were no ERRORs or WARNINGs. There was one NOTE, that this is a new submission and
some words are possibly misspelled. The words flagged are names and spelled correctly.


## CRAN Submission 1.0.5
Based on feedback received, for submission 1.0.5 I have:
* removed redundant "Tools to" in Description of DESCRIPTION
* corrected the format of the reference in Description of DESCRIPTION

For submission 1.0.4, I have: 
* removed a few of the tests that require many iterations for model types that are most similar to other methods being tested. I have also reduced the number of iterations in other tests.

For submission 1.0.3, I have:
* fixed a typo in function documentation

For submission 1.0.2, after failing the incoming checks I have:
* removed two tests that check for the estimates being in between the confidence 
interval limits for geeglm and lme models. The interval methods for the longitudinal 
models in this package have not been evaluated for performance and there are already 
warnings when producing intervals for these models. Since there is already a disclaimer 
for this, these tests have been removed.

For submission 1.0.1, after failing the incoming checks I have:
* shortened the time it take the examples in resi.R to run

## Downstream dependencies
There are currently no downstream dependencies for this package.
