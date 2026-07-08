## CRAN Submission 1.4.0

### New features
- New asymptotic CI methods (`ci.method = "normal"`, `"qf"`, `"cf"`) via `resi_pe_asymptotic()`
- Cross-sectional (CS-RESI) and longitudinal (L-RESI) estimates for `lmerMod` models
- New simulation and calibration functions: `insurancePlasmodeSim()`, `simFigures()`, etc.

### Bug fixes
- Fixed intercept-only reduced model for computed responses in `resi_pe.default`

## Test Environments
* local R installation, R 4.x.x (macOS)
* rhub::rhub_check() (linux, macos, windows, atlas)
* win-builder

## R CMD check results
There were no ERRORs, WARNINGs. There was a NOTE regarding timestamps, unable to verify current time.
I resolve a global variable binding issues using utils::globalVariables, as previously recommended by the CRAN team.

## Downstream dependencies
There are currently no downstream dependencies for this package.
