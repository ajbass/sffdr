# sffdr 1.1.4

* The density window in `sffdr()` now spans the documented ~5,000 neighbours.
  `nn` was taken as a fraction of all SNPs but applied to the rows the density
  is fitted on (the tail plus 100,000 downsampled null SNPs), which left about
  600 neighbours at one million SNPs. Results at up to ~100,000 SNPs are
  unchanged; above that, calls at fq <= 0.05 move by up to 2% (fewer at one
  million SNPs) at the same ranking and the same type I error with independent
  SNPs.
* Perfect-LD copies (SNPs with identical p-values) now contribute one point
  each to the density fit, with the density evaluated for every SNP. With many
  copies the narrow window inflated functional p-values below 1e-3 by up to
  2.4-fold at one million SNPs; it is now within the level of the same number
  of independent SNPs. Ties from rounded p-values, which sit almost only in
  the bulk of the distribution, are recognised by that and left alone.
* `fpvalues()` gives tests with identical `lfdr` and `p` the same values (the
  set's), so the result no longer depends on the input order; negative `lfdr`
  is an error and all-zero `lfdr` gives functional p-values of 0 instead of NaN.
* `kernelEstimator()` no longer crashes R on tied or near-tied values in the
  univariate case (tied values are collapsed, the fit uses at most 20
  iterations and at least 10 distinct neighbours). Its `maxk` default is 5,000
  instead of 500,000 (35 MB instead of 2.7 GB of evaluation space), and the fit
  is repeated with more space if locfit runs out; the fit is unchanged. `...`
  is passed to `locfit()`, as it always was, and the arguments set by the
  estimator cannot be passed through it. The error after a failed cascade
  reports locfit's last message.
* `pi0_model()` and `fpi0est()` accept `indep_snps` as a logical vector, row
  indices or SNP ids (row names). Row indices were summed as if logical in
  `pi0_model()` (far too many SNPs per knot) and were recycled in `fpi0est()`.
* `pi0_model()` stops on non-numeric input (an rsID column was ranked as text
  and modelled) and warns when a column is not a p-value column or `qvalue()`
  fails, instead of silently using the fallback knot.
* `fpi0est()` checks the lengths of `p`, `weights` and `indep_snps`, requires
  `lambda` in (0, 1), warns when the chosen model did not converge and when
  weights above 1 are clamped (LD scores passed in place of inverse LD scores
  gave an unweighted fit). Tests with a missing p-value or covariate get
  fpi0 = 1 and the floor of 0.01 are documented. `tableLambda` no longer holds
  the 19 fitted models (~1 GB per million SNPs); the chosen model is returned
  as `model`.
* `sffdr()` stops on negative p-values (a log10 p-value column became the
  strongest discovery) and clamps exact zeros only.
* `plot.sffdr()` showed the wrong SNPs' fpi0 in panel (a) when some p-values
  were missing; also removed a duplicated scale and a deprecated legend
  position.
* The C++ smoothers check their input lengths and ordering.
* Documentation: monotone smoothing has been the default since 1.1.2; the
  bivariate `kernelEstimator()` takes the covariate first and the p-values
  second (the example and tests had them swapped).
* `tests/testthat.R` added, so `R CMD check` runs the test suite (it did not
  before).

# sffdr 1.1.3

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
* `sffdr()` no longer uses `weights` (a warning is given). Its density is
  estimated over all tested SNPs, the set the functional p-values count errors
  in; reweighting it by inverse LD scores biased the local FDR, because
  associated SNPs have larger LD scores (conservative with inverse-LD weights,
  anti-conservative with weights larger at associated SNPs). LD scores belong
  in the functional pi0 model as a covariate; inverse-LD weights remain
  available in `pi0_model()` and `fpi0est()` in place of LD pruning. Results
  without `weights` are unchanged.

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
