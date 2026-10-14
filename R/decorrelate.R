# Sample overlap: decorrelated informative studies ----------------------------

#' Decorrelate the informative studies from the primary study
#'
#' @description
#' When the primary study shares samples with the informative studies, their
#' null z-scores are correlated. Used as covariates in \code{\link{pi0_model}},
#' the informative p-values then predict the primary's null tests, which
#' inflates false discoveries. \code{decorrelate_informative()} removes the
#' primary's noise from each informative study,
#' \deqn{u_k = (z_k - r_k z_1) / \sqrt{1 - r_k^2},}
#' where \eqn{r_k} is the null correlation of study \eqn{k} with the primary.
#' Under the primary's null, \eqn{u} is independent of \eqn{z_1} whatever the
#' informative studies' effects, so the primary's p-values are analysed
#' unchanged with the decorrelated p-values as covariates.
#'
#' @details
#' The correction needs only the overlap correlation of each informative study
#' with the primary: cross-trait LD score regression intercepts are the
#' recommended source for GWAS; otherwise they are estimated by
#' \code{\link{estimate_overlap}}. Errors in them enter at second order: in
#' simulations, a correlation of 0.5 given as 0.35 or 0.65 left the FDR within
#' 0.005 of its value at the true correlation.
#'
#' @param z_x Numeric vector; the primary study's z-scores.
#' @param z_y Numeric vector or matrix; signed z-scores of the informative
#'   studies (one column per study), in the same order as \code{z_x}.
#' @param p_x,p_y Ignored; kept for compatibility with sffdr 1.1.2.
#' @param rho_o Optional overlap correlations of the primary with each
#'   informative study (one number per study). If \code{NULL} (default),
#'   estimated with \code{\link{estimate_overlap}}.
#' @param lfdr_x Ignored. sffdr 1.1.2 scaled the correction by the primary's
#'   local FDR, which left part of the overlap in the covariates.
#' @param indep_snps Optional logical or integer vector selecting approximately
#'   independent SNPs, used to estimate \code{rho_o}.
#'
#' @return A list with
#' \describe{
#'   \item{z2}{decorrelated z-scores of the informative studies (vector or matrix, as \code{z_y})}
#'   \item{p2}{their p-values, for use as covariates}
#'   \item{rho_o}{the overlap correlations used}
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 20000
#' e1 <- rnorm(n)
#' e2 <- 0.5 * e1 + sqrt(1 - 0.5^2) * rnorm(n)
#' z_x <- e1 + ifelse(runif(n) < 0.05, 3, 0)            # primary
#' z_y <- e2 + ifelse(runif(n) < 0.1, 3, 0)             # informative, same samples
#' dec <- decorrelate_informative(z_x, z_y, rho_o = 0.5)
#' \donttest{
#' p_x <- 2 * pnorm(-abs(z_x))
#' mpi0 <- pi0_model(as.matrix(dec$p2))
#' fpi0 <- fpi0est(p_x, pi0_model_obj = mpi0)
#' res <- sffdr(p_x, fpi0 = fpi0$fpi0)
#' }
#'
#' @seealso \code{\link{estimate_overlap}}
#' @export
decorrelate_informative <- function(
  z_x,
  z_y,
  p_x = NULL,
  p_y = NULL,
  rho_o = NULL,
  lfdr_x = NULL,
  indep_snps = NULL
) {
  zy <- as.matrix(z_y)
  if (!is.numeric(z_x) || !is.numeric(zy)) {
    stop("'z_x' and 'z_y' must be numeric.")
  }
  n <- length(z_x)
  if (nrow(zy) != n) {
    stop("'z_y' must have one row (or element) per entry of 'z_x'.")
  }
  K <- ncol(zy)
  if (!is.null(lfdr_x)) {
    warning("decorrelate_informative: 'lfdr_x' is no longer used; the correction is not scaled by the local FDR.")
  }
  ok <- is.finite(z_x) & rowSums(!is.finite(zy)) == 0
  both <- function(x) any(x > 0) && any(x < 0)
  if (sum(ok) >= 20 && (!both(z_x[ok]) || !all(apply(zy[ok, , drop = FALSE], 2, both)))) {
    stop("'z_x' and 'z_y' must be signed z-scores (a study has only one sign).")
  }

  ## overlap correlations with the primary
  if (is.null(rho_o)) {
    rho_o <- estimate_overlap(z_x, zy, indep_snps = indep_snps)[1, -1]
  }
  if (!is.numeric(rho_o) || length(rho_o) != K || anyNA(rho_o) || any(abs(rho_o) >= 1)) {
    stop("'rho_o' must give one correlation in (-1, 1) per informative study.")
  }
  rho_o <- unname(rho_o)

  ## the primary's noise removed: independent of z_x under the primary's null
  u <- sweep(zy, 2, rho_o, function(z, r) z - r * z_x)
  u <- sweep(u, 2, sqrt(1 - rho_o^2), "/")
  p <- pmin(pmax(exp(log(2) + pnorm(-abs(u), log.p = TRUE)), .Machine$double.xmin), 1)
  u[!ok, ] <- NA_real_
  p[!ok, ] <- NA_real_
  if (is.null(dim(z_y))) {
    u <- drop(u)
    p <- drop(p)
  }
  list(z2 = u, p2 = p, rho_o = rho_o)
}
