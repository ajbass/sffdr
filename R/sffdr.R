#' @title
#' Estimate the functional p-values, q-values, and local false discovery rates given a set of p-values and
#' informative variables
#'
#' @description
#' Estimate the functional p-values, q-values, and local false discovery rates given a set of p-values and informative
#' variables. The functional p-values is mapping from the q-value (FDR-based measure) to a p-value for family-wise error rate control.
#'
#' @details
#' The function \code{\link{fpi0est}} should be called externally to estimate the
#' functional proportion of null tests given the set of informative variables.
#' The surrogate functional FDR methodology builds from the functional FDR (reference below)
#' methodology and implements some of the functions from the package (denoted in documentation).
#' Note that if there is only a single informative variable then the function reduces to fFDR framework.
#'
#' @param p.value A vector of p-values.
#' @param fpi0 An estimate of the function proportion of null tests using the \code{\link{fpi0est}} function.
#' @param surrogate A surrogate variable that compresses more than one informative variables.
#' Default is NULL. If \code{fpi0} is specified and \code{surrogate} is NULL then \code{fpi0} is used as the surrogate variable.
#' @param indep_snps A boolean vector (same size as p) specifying the set of independent tests. Default is NULL and all tests are treated independently.
#' @param monotone.window Enforce monotonicity at specified step size. See fFDR paper below. Default is NULL.
#' @param epsilon A numerical value the truncation for the p-values during density estimation. Default is 1e-15. You may want to consider decreasing this value if there are a substantial number of small p-values.
#' @param nn A numerical value specifying the nearest neighbor parameter in \code{\link{kernelEstimator}}. Default is NULL.
#' @param fp_ties A boolean specifying whether ties should be broken using the ordering of the p-values when calculating the fp-values. Only impacts the tests when the local FDR is tied. Default is TRUE.
#' @param \ldots Additional arguments passed to \code{\link{kernelEstimator}}.
#'
#'
#' @return
#' A list of object type "sffdr" containing:
#' \item{pvalues}{A vector of the original p-values.}
#' \item{fpvalues}{A vector of the estimated functional p-values.}
#' \item{fqvalues}{A vector of the estimated functional q-values.}
#' \item{flfdr}{A vector of the estimated functional local FDR values.}
#' \item{pi0}{An vector of the original functional proportion of null tests.}
#' \item{density}{An object containing the kernel density estimates from \code{kernelEstimator}.}
#'
#' @examples
#' # import data
#' data(bmi)
#'
#' # separate main p-values and conditioning p-values
#' p <- sumstats$bmi
#' z <- as.matrix(sumstats[, -1])
#'
#' # apply pi0_model to create model
#' knots <- c(0.005, 0.01, 0.025, 0.05, 0.1)
#' fmod <- pi0_model(z, knots = knots)
#'
#' # estimate functional pi0
#' fpi0_out <- fpi0est(p, z = fmod$zt, pi0_model = fmod$fmod)
#' fpi0 <- fpi0_out$fpi0
#'
#' # apply sffdr
#' # Note all tests are independent see 'indep_snps' argument
#' # The data has very small p-values, set epsilon to min of p
#' sffdr_out <- sffdr(p, fpi0, epsilon = min(p))
#'
#' # Plot significance results
#' plot(sffdr_out, rng = c(0, 5e-4))
#'
#' # Functional P-values, Q-values, and local FDR
#' fp <- sffdr_out$fpvalues
#' fq <- sffdr_out$fqvalues
#' flfdr <- sffdr_out$flfdr
#'
#' @author Andrew J. Bass
#' @seealso \code{\link{fpi0est}}, \code{\link{plot.sffdr}}
#' @keywords sffdr
#' @aliases sffdr
#' @import Rcpp
#' @export
sffdr <- function(p.value,
                  fpi0,
                  surrogate = NULL,
                  indep_snps = NULL,
                  monotone.window = NULL,
                  epsilon = 1e-15,
                  nn = NULL,
                  fp_ties = TRUE, ...) {

  if (is.null(indep_snps)) {
    indep_snps <- rep(TRUE, length(p.value))
    indep.check <- NULL
  } else {
    if (length(indep_snps) != length(p.value)){
      stop("Length of independent SNPs different than p-values")
    } else if (sum(indep_snps) < 100) {
      warning("Less than 100 independent SNPs. Estimation will be noisy.")
    }
    indep.check <- TRUE
  }

  if (is.null(surrogate)) {
    surrogate <- fpi0
  }
  # Rank transformation surrogate
  z <- rank(surrogate, ties.method = "random") / length(surrogate)

  # NN for density estimator
  if (is.null(nn)) {
    pi1 <- 1 - mean(fpi0)
    nn <- min(pi1, 1 - pi1)
    if (nn < 0.02) nn <- 0.02
  }

  kd <- kernelEstimator(cbind(z[indep_snps], p.value[indep_snps]),
                        nn = nn,
                        eval.points = cbind(z, p.value), epsilon = epsilon, ...)

  if (!is.null(monotone.window)) {
    if (length(p.value) > 500000) {
      warning("Enforcing monotonicity may take a few minutes.")
    }
    sp <- order(p.value)
    rank.p <- rank(p.value, ties.method = "random")
    original.fx <- kd$fx
    fx = monoSmooth(p.value[sp],  z[sp], original.fx[sp], monotone.window)
    fx <- fx[rank.p]
  } else {
    fx <- kd$fx
  }

  if (!is.null(indep.check)) {
    # account for LD
    kd_surrogate <- kernelEstimator(z[indep_snps],
                                    nn = nn,
                                    eval.points = z,
                                    epsilon = epsilon, ...)
    sfx <- kd_surrogate$fx
  } else {
    sfx <-  1
  }
  # Local FDR: sfx marginal density of surrogate, fx joint density
  lfdr <- pmin(fpi0 * sfx / fx, 1)

  if (!fp_ties) {
    fpq <- fpvalues(lfdr)
    fp <- fpq$fp
    fq <- fpq$fq
  } else {
    fpq <- fpvalues(lfdr, p.value)
    fp <- fpq$fp
    fq <- fpq$fq
  }

  ret <- list(call = match.call(),
              pvalues = p.value,
              fpvalues = fp,
              fqvalues = fq,
              flfdr = lfdr,
              fpi0 = fpi0,
              density = kd)
  class(ret) <- "sffdr"
  ret
}

#' @import qvalue
#' @import locfit
#' @import splines
#' @import dplyr
#' @import gam
#' @import qvalue
NULL
