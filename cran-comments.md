## CRAN Submission 1.4.2

### New features

- New asymptotic CI methods (`ci.method = "normal"`, `"qf"`, `"cf"`) via `resi_pe_asymptotic()`
- Cross-sectional (CS-RESI) and longitudinal (L-RESI) estimates for `lmerMod` models
- New simulation and calibration functions: `insurancePlasmodeSim()`, `simFigures()`, etc.

### Bug fixes

- Fixed intercept-only reduced model for computed responses in `resi_pe.default`
- Fixed derivative method to produce accurate coverage.
- Added statistical evaluation to package website.

## Test Environments

- local R installation, R 4.x.x (macOS)
- rhub::rhub_check() (linux, macos, windows, atlas)
- win-builder

## R CMD check results

There were no ERRORs, WARNINGs. There was a NOTE regarding URLs in /inst/doc/RESI_paper.pdf. The URLs reported in the NOTE originate from citations embedded in this PDF. The PDF is not generated during package build and the package code, examples, and vignettes do not depend on these URLs. The NOTE is confined to this archival document.

## Downstream dependencies

There are currently no downstream dependencies for this package.
