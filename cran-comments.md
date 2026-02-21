## R CMD check results

0 errors | 0 warnings | 0 note

## Test environments
* local: macOS (R 4.4.3)
* GitHub Actions (ubuntu-latest): R-release, R-devel
* win-builder: R-release, R-devel

## Submission notes

This is a new release (v1.1.0) of the sffdr package.

### Changes in this version
* Improved computational efficiency by removing unnecessary dependencies (dplyr, tidyr, tibble)
* Enhanced robustness with comprehensive NA handling and input validation
* Added reproducibility via `seed` argument in main functions
* Improved algorithm stability in `fpi0est`, `pi0_model`, and `kernelEstimator` functions
* Added progress bars for long-running operations
* Updated documentation and vignettes

### Notes
There are no downstream dependencies to check.

## Reverse dependencies

There are currently no reverse dependencies for this package.