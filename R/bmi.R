#' @name bmi
#' @aliases bmi
#' @title Subset of p-values from the UK Biobank analysis
#'
#' @usage data(bmi)
#'
#' @description The summary level data is a subset of independent SNPs
#' from the UK Biobank where we performed a GWAS of body mass index (BMI), body fat percentage (BFP),
#'cholesterol, and triglycerides. Note that BFP, cholesterol and triglycerides
#' are conditioning traits and were calculated using a separate set of individuals
#' than BMI. See manuscript for details.
#'
#' @return A list called \code{sumstats} containing:
#' \item{bmi}{Vector of 10,000 p-values for BMI.}
#' \item{bfp}{Vector of 10,000 p-values for BFP.}
#' \item{cho}{Vector of 10,000 p-values for cholesterol.}
#' \item{tri}{Vector of 10,000 p-values for triglycerides.}
#'
#' @seealso \code{\link{sffdr}}
#' @aliases sumstats
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
#' # The very small p-values, set epsilon to min of p
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
#' @keywords dataset, bmi
NULL
