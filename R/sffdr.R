#' Surrogate Functional False Discovery Rate Analysis
#'
#' @description
#' Estimate functional p-values, q-values, and local false discovery rates (lfdr)
#' for GWAS data leveraging summary statistics from related traits. Functional p-values map from the
#' functional q-value (FDR-based measure) to a p-value for type I error rate control,
#' accounting for pleiotropy that impacts the prior probability of
#' association.
#'
#' @details
#' The surrogate functional FDR (sfFDR) methodology extends the functional FDR
#' framework to leverage multiple informative variables (e.g., functional annotations,
#' GWAS summary statistics) for increased power while controlling
#' false discovery rates.
#'
#' Workflow:
#'
#' \enumerate{
#'   \item Estimate functional pi0 (proportion of nulls) using \code{\link{fpi0est}}
#'   \item Call \code{sffdr()} with p-values and estimated functional pi0
#'   \item Use returned functional p-values/q-values/local FDRs for significance testing
#' }
#'
#' Surrogate Variable:
#' If not specified, the estimated functional pi0 is used as the surrogate variable.
#'
#' Weights and LD:
#' The joint density is estimated over all tested SNPs, which is the set the
#' functional p-values count errors in, so \code{sffdr()} no longer reweights
#' it. Inverse-LD weights change the density to that of independent signals;
#' because associated SNPs have larger LD scores, this biases the local FDR
#' (conservative with inverse-LD weights, anti-conservative with weights that
#' are larger at associated SNPs). LD scores carry information about
#' association, so their place is the functional pi0 model, as a covariate
#' (see the vignette); as inverse weights they remain available in
#' \code{\link{pi0_model}} and \code{\link{fpi0est}} in place of an LD-pruned
#' training set.
#'
#' Interpretation note: Functional p-values reflect the global ranking of all
#' SNPs by the functional local FDR. A SNP with a small functional p-value but
#' large flfdr (> 0.5) has weak individual evidence but ranks favorably
#' relative to null SNPs. This happens when there is strong prior evidence (fpi0) of the SNP
#' being non-null. Thus, users should also consider flfdr (and fqvalues) alongside fpvalues
#' when assessing individual SNP evidence. SNPs with flfdr > 0.5 should
#' be interpreted cautiously regardless of their functional p-value.
#'
#' @param p.value Numeric vector of p-values to analyze.
#' @param fpi0 Numeric vector of functional pi0 estimates, obtained
#'   from \code{\link{fpi0est}}. Values must be in [0, 1].
#' @param surrogate Optional numeric vector (same length as \code{p.value}) used
#'   as a surrogate variable for compression. If NULL (default), uses
#'   \code{fpi0} as the surrogate.
#' @param weights Ignored, with a warning (kept for compatibility with sffdr
#'   1.1.2); see Details.
#' @param epsilon Numeric; lower bound for p-value clamping during density estimation.
#'   Default is \code{.Machine$double.xmin}.
#' @param nn Numeric; nearest-neighbor bandwidth for \code{\link{kernelEstimator}}.
#'   If NULL (default), automatically selected as ~5000 neighbors.
#' @param monotone Logical; if TRUE, enforces monotonicity of the estimated density with respect to p-values within bins of the surrogate variable. Default is TRUE.
#' @param monotone_method Character; method for enforcing monotonicity. Options are "min" (running minimum) or "pava" (Pool Adjacent Violators Algorithm). Default is "min".
#'   We do not recommend "pava" for GWAS data sets due to LD, but it may be better for expression data.
#' @param fp_ties Logical; whether to break ties in functional p-values using the
#'   original p-value ordering. Default is TRUE.
#' @param seed Integer; random seed for reproducibility of rank tie-breaking.
#'   Default is 2026.
#' @param verbose Logical; print progress messages. Default is TRUE.
#' @param ... Additional arguments passed to \code{\link{kernelEstimator}}.
#'
#' @return
#' An S3 object of class \code{"sffdr"} containing:
#' \item{call}{The function call.}
#' \item{pvalues}{Original p-values.}
#' \item{fpvalues}{Functional p-values.}
#' \item{fqvalues}{Functional q-values.}
#' \item{flfdr}{Functional local false discovery rates.}
#' \item{fpi0}{Functional pi0 estimates.}
#' \item{fx}{Joint density estimates at observed (p-value, surrogate) pairs.}
#'
#' @examples
#' \dontrun{
#' # Import data
#' data(bmi)
#'
#' # Separate main p-values and conditioning p-values
#' p <- sumstats$bmi
#' z <- as.matrix(sumstats[, -1])
#'
#' # Apply pi0_model to create model
#' # (note: use indep_snps argument to specify independent SNPs for training)
#' fmod <- pi0_model(z)
#'
#' # Estimate functional pi0
#' # (note: use indep_snps argument to specify independent SNPs for training)
#' fpi0_out <- fpi0est(p, fmod)
#' fpi0 <- fpi0_out$fpi0
#'
#' # Apply sffdr
#' sffdr_out <- sffdr(p, fpi0)
#'
#' # Plot significance results
#' plot(sffdr_out)
#'
#' # Extract functional quantities
#' fp <- sffdr_out$fpvalues
#' fq <- sffdr_out$fqvalues
#' flfdr <- sffdr_out$flfdr
#' }
#'
#' @author Andrew J. Bass
#' @seealso \code{\link{fpi0est}}, \code{\link{plot.sffdr}}, \code{\link{kernelEstimator}}
#' @keywords sffdr
#' @aliases sffdr
#' @export
sffdr <- function(
  p.value,
  fpi0,
  surrogate = NULL,
  weights = NULL,
  epsilon = .Machine$double.xmin,
  nn = NULL,
  monotone = TRUE,
  monotone_method = c("min", "pava"),
  fp_ties = TRUE,
  seed = 2026,
  verbose = TRUE,
  ...
) {
  if (verbose) {
    message("==================================================")
    message(sprintf(
      "Running sfFDR Analysis (v%s)",
      utils::packageVersion("sffdr")
    ))
    message("==================================================")
  }

  # Set random seed for reproducibility
  if (!is.null(seed)) {
    withr::local_seed(seed)
  }

  # Input validation
  n <- length(p.value)
  monotone_method <- match.arg(monotone_method, c("min", "pava"))

  if (!is.numeric(p.value)) {
    stop("'p.value' must be a numeric vector.")
  }

  if (!is.numeric(fpi0)) {
    stop("'fpi0' must be a numeric vector.")
  }

  if (length(fpi0) != n) {
    stop("Length of 'fpi0' must match 'p.value'.")
  }

  if (!is.null(surrogate) && length(surrogate) != n) {
    stop("Length of 'surrogate' must match 'p.value'.")
  }

  if (any(p.value <= 0, na.rm = TRUE)) {
    warning("Found p-values <= 0. Clamping to machine epsilon.")
    p.value <- pmax(p.value, .Machine$double.xmin, na.rm = FALSE)
  }

  if (any(p.value > 1, na.rm = TRUE)) {
    stop("Found p-values > 1.")
  }

  if (any(fpi0 < 0 | fpi0 > 1, na.rm = TRUE)) {
    stop("'fpi0' values must be in [0, 1].")
  }

  # Handle NAs: identify valid indices, compute on valid subset, map back
  valid <- !is.na(p.value) & !is.na(fpi0)
  if (!is.null(surrogate)) {
    valid <- valid & !is.na(surrogate)
  }

  if (!any(valid)) {
    stop("No non-NA observations.")
  }

  # Work on valid subset
  p_valid <- p.value[valid]
  fpi0_valid <- fpi0[valid]
  n_valid <- sum(valid)

  # Weights no longer enter the density: it is estimated over all tested SNPs
  if (!is.null(weights)) {
    warning(
      "sffdr: 'weights' are ignored; the density is estimated over all tested SNPs. ",
      "Use 'weights' in pi0_model() and fpi0est(), or include the LD score as a covariate."
    )
  }

  # Set surrogate variable
  if (verbose) {
    message("  Transforming surrogate variable...")
  }

  surrogate_var <- if (is.null(surrogate)) fpi0_valid else surrogate[valid]

  # Rank-transform surrogate to percentile scale
  z <- rank(surrogate_var, ties.method = "random") / n_valid

  # Fit joint density
  if (verbose) {
    message("  Fitting kernel density...")
  }

  kd <- kernelEstimator(
    x = cbind(z, p_valid),
    eval.points = cbind(z, p_valid),
    nn = nn,
    epsilon = epsilon,
    weights = NULL,
    verbose = verbose,
    ...
  )
  fx_valid <- pmax(kd$fx, .Machine$double.xmin)

  if (monotone) {
    if (verbose) {
      msg <- sprintf(
        "  Enforcing conditional monotonicity using the '%s' method...",
        monotone_method
      )
      message(msg)
    }

    # Target ~1,000 SNPs per bin
    n_bins <- as.integer(max(10, min(1000, floor(n_valid / 1000))))

    if (verbose) {
      message(sprintf(
        "    Dividing data into %d conditional surrogate bins...",
        n_bins
      ))
    }

    groups <- as.integer(cut(z, breaks = n_bins, labels = FALSE))
    ord <- order(groups, p_valid)

    # Route to the appropriate C++ engine based on user choice
    if (monotone_method == "min") {
      smoothed_fx <- monoSmooth_conditional(
        pvalue = p_valid[ord],
        density = fx_valid[ord],
        group = groups[ord]
      )
    } else {
      smoothed_fx <- monoSmooth_pava(
        pvalue = p_valid[ord],
        density = fx_valid[ord],
        group = groups[ord]
      )
    }

    # Invert the order to map back to valid subset
    fx_valid[ord] <- smoothed_fx
  }

  # Compute functional local FDR
  if (verbose) {
    message("  Computing functional local FDRs...")
  }

  lfdr_valid <- pmin(fpi0_valid / fx_valid, 1)

  # Compute functional p-values and q-values
  if (verbose) {
    message("  Computing functional p-values and q-values...")
  }

  fpq <- if (!fp_ties) {
    fpvalues(lfdr_valid)
  } else {
    fpvalues(lfdr_valid, p_valid)
  }

  if (verbose) {
    message("  Done.")
    message("==================================================")
    message("")
  }

  # Map back to full length, NAs for invalid positions
  fpvalues_out <- rep(NA_real_, n)
  fqvalues_out <- rep(NA_real_, n)
  flfdr_out <- rep(NA_real_, n)
  fx_out <- rep(NA_real_, n)

  fpvalues_out[valid] <- fpq$fp
  fqvalues_out[valid] <- fpq$fq
  flfdr_out[valid] <- lfdr_valid
  fx_out[valid] <- fx_valid

  # Return results
  structure(
    list(
      call = match.call(),
      pvalues = p.value,
      fpvalues = fpvalues_out,
      fqvalues = fqvalues_out,
      flfdr = flfdr_out,
      fpi0 = fpi0,
      fx = fx_out
    ),
    class = "sffdr"
  )
}

#' @importFrom qvalue qvalue pi0est
#' @importFrom locfit locfit lp rbox
#' @importFrom splines ns
#' @useDynLib sffdr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
