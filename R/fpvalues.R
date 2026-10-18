#' @title
#' Functional p-values
#'
#' @description
#' Calculate functional p-values from functional local FDRs.
#'
#' @details
#' Tests are ranked by \code{lfdr}, with ties broken by \code{p}. The
#' functional q-value of a test is the mean \code{lfdr} of the tests ranked at
#' or above it, and its functional p-value is the \code{lfdr} mass ranked at
#' or above it as a fraction of the total mass. Tests with identical
#' \code{lfdr} and \code{p} are one set: each receives the set's values, so
#' the result does not depend on the order of the input. When \code{p} is
#' NULL, ties are broken at random.
#'
#' @param lfdr A vector of functional local FDRs in [0, 1] (values above 1 are
#'   clamped to 1).
#' @param p A vector of p-values for ranking purposes. Default is NULL.
#'
#' @return
#' A list containing:
#' \item{fp}{Functional p-values.}
#' \item{fq}{Functional q-values.}
#' @export
fpvalues <- function(lfdr, p = NULL) {
  if (!is.numeric(lfdr)) {
    stop("'lfdr' must be a numeric vector.")
  }
  if (any(lfdr < 0, na.rm = TRUE)) {
    stop("'lfdr' must be in [0, 1].")
  }
  if (is.null(p)) {
    # create random ordering if p-values are not provided
    p <- runif(length(lfdr))
  }

  if (length(lfdr) != length(p)) {
    stop("Length of 'lfdr' and 'p' must match.")
  }

  n <- length(p)
  p_out <- rep(NA_real_, n)
  fdr_out <- rep(NA_real_, n)

  rm_na <- !is.na(p) & !is.na(lfdr)
  if (sum(rm_na) == 0) {
    out <- list(fp = p_out, fq = fdr_out)
    return(out)
  }

  p_sub <- p[rm_na]
  lfdr_sub <- lfdr[rm_na]
  lfdr_sub[lfdr_sub > 1] <- 1

  o <- order(lfdr_sub, p_sub)
  n_sub <- length(o)
  lfdr_o <- lfdr_sub[o]
  p_o <- p_sub[o]

  # Calculate FDR on sorted data
  fdr_sorted <- cumsum(lfdr_o) / seq_len(n_sub)
  max_fdr <- max(fdr_sorted)

  p_sorted <- if (max_fdr > 0) {
    pmin((fdr_sorted * (seq_len(n_sub) / n_sub)) / max_fdr, 1)
  } else {
    rep(0, n_sub) # no null mass anywhere
  }

  # Tests with identical (lfdr, p) are one set and take the set's values
  if (n_sub > 1) {
    new_set <- c(TRUE, lfdr_o[-1] != lfdr_o[-n_sub] | p_o[-1] != p_o[-n_sub])
    if (!all(new_set)) {
      set_id <- cumsum(new_set)
      last <- cumsum(tabulate(set_id))
      fdr_sorted <- fdr_sorted[last][set_id]
      p_sorted <- p_sorted[last][set_id]
    }
  }

  idx <- which(rm_na)
  fdr_out[idx[o]] <- fdr_sorted
  p_out[idx[o]] <- p_sorted

  out <- list(fp = p_out, fq = fdr_out)
  return(out)
}
