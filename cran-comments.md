## R CMD check results

0 errors | 0 warnings | 0 note

## Test environments
* local: macOS (R 4.4.3)
* GitHub Actions (ubuntu-latest): R-release, R-devel
* win-builder: R-release, R-devel

## Submission notes

This is a new release (v1.1.2) of the sffdr package.

### Changes in this version

* Added monotonicity option (running min or PAVA)
* Better handling of weights in kernel density estimation
* function docorrelate_informative to handle overlapping samples
* Added a weights argument which allows users to input inverse LD scores and fit to all SNPs.

### Notes
There are no downstream dependencies to check.

## Reverse dependencies

There are currently no reverse dependencies for this package.