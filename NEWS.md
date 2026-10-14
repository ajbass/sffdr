# sffdr (development version)

* `decorrelate_informative()` now removes the primary's noise from the
  informative studies exactly, (z_y - rho z_x) / sqrt(1 - rho^2), instead of
  scaling the correction by the primary's local FDR, which left part of the
  overlap in the covariates and inflated false discoveries at larger overlap.
  It takes several informative studies at once and, without `rho_o`, estimates
  the overlap with the new `estimate_overlap()`. `p_x`, `p_y` and `lfdr_x` are
  no longer used.
* New `estimate_overlap()`: the null correlation of z-scores between studies
  from their products near zero.
* `overlap_null_density()`, `run_empirical_overlap_correction()`,
  `estimate_rho_overlap()` and `discover_empirical_nulls()` are defunct.

# sffdr 1.1.2

* Added monotonicity option (running min or PAVA)
* Better handling of weights in kernel density estimation
* function docorrelate_informative to handle overlapping samples

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
