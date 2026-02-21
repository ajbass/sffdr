# sffdr (development version)

# sffdr 1.1.1

* Added a weights argument which allows users to input inverse LD scores and fit to all SNPs.

# sffdr 1.1.0

## Major Changes

*  Now uses `fastglm` exclusively for significantly faster GLM fitting.

* Simplified `pi0_model()` function:

  - Removed `knots` parameter. Knots are now selected adaptively based on FDR 
    thresholds and signal strength.
  - Added `n_knots` parameter (default: 5) for target knot count.
  - Added `min_snps_per_knot` parameter (default: 100) to ensure stable spline fitting.
  - Added `max_knots` parameter (default: 5) to prevent overfitting.

* Removed `indep_snps` functionality from `sffdr()`. Kernel density estimation now always uses all SNPs with appropriate `nn` smoothing to handle LD.

* Improved speed of locfit by subsampling the null region of the distribution (accounted for using the weights argument).

* Added 98 unit tests covering `kernelEstimator`, `pi0_model`, `fpi0est`,  `sffdr`, and `fpvalues`.

# sffdr 1.0.0

* Initial CRAN submission.
