# Defunct functions ------------------------------------------------------------

#' Defunct functions in sffdr
#'
#' The earlier sample-overlap functions have been removed. Use
#' \code{\link{decorrelate_informative}} to correct for sample overlap with
#' the informative studies, and
#' \code{\link{estimate_overlap}} to estimate the overlap correlation.
#'
#' \itemize{
#'   \item \code{overlap_null_density}: use \code{\link{decorrelate_informative}}
#'   \item \code{run_empirical_overlap_correction}: use \code{\link{decorrelate_informative}}
#'   \item \code{estimate_rho_overlap}: use \code{\link{estimate_overlap}}
#'   \item \code{discover_empirical_nulls}: use \code{\link{estimate_overlap}}
#' }
#'
#' @param ... Ignored.
#' @name sffdr-defunct
#' @keywords internal
NULL

.overlap_defunct <- function(old, new) {
  .Defunct(new, package = "sffdr", msg = sprintf(
    "'%s' has been removed from sffdr. Use '%s' instead (see ?%s).", old, new, new
  ))
}

#' @rdname sffdr-defunct
#' @export
overlap_null_density <- function(...) .overlap_defunct("overlap_null_density", "decorrelate_informative")

#' @rdname sffdr-defunct
#' @export
run_empirical_overlap_correction <- function(...) .overlap_defunct("run_empirical_overlap_correction", "decorrelate_informative")

#' @rdname sffdr-defunct
#' @export
estimate_rho_overlap <- function(...) .overlap_defunct("estimate_rho_overlap", "estimate_overlap")

#' @rdname sffdr-defunct
#' @export
discover_empirical_nulls <- function(...) .overlap_defunct("discover_empirical_nulls", "estimate_overlap")
