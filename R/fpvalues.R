.fpvalues_raw <- function(lfdr) {
  n <- length(lfdr)
  p_out <- rep(NA_real_, n)
  fdr_out <- rep(NA_real_, n)

  rm_na <- !is.na(lfdr)
  if (!any(rm_na)) {
    return(list(fp = p_out, fq = fdr_out))
  }

  lfdr_clean <- lfdr[rm_na]
  lfdr_clean[lfdr_clean > 1] <- 1

  o <- order(lfdr_clean)
  lfdr_sorted <- lfdr_clean[o]

  # Calculate cumulative mean (FDR)
  fdr_sorted <- cumsum(lfdr_sorted) / seq_along(lfdr_sorted)

  # Rank proportion
  rank_prop <- seq_along(lfdr_sorted) / length(lfdr_sorted)

  # P-value calculation
  p_sorted <- (fdr_sorted * rank_prop) / max(fdr_sorted)
  p_sorted <- pmin(p_sorted, 1)

  p_final <- numeric(length(lfdr_clean))
  fdr_final <- numeric(length(lfdr_clean))

  p_final[o] <- p_sorted
  fdr_final[o] <- fdr_sorted

  p_out[rm_na] <- p_final
  fdr_out[rm_na] <- fdr_final

  return(list(fp = p_out, fq = fdr_out))
}

#' @title
#' Functional p-values
#'
#' @description
#' Calculate functional p-values from functional local FDRs.
#'
#' @param lfdr A vector of functional local FDRs.
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
  if (is.null(p)) {
    out <- .fpvalues_raw(lfdr)
    return(out)
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

  r <- integer(length(o))
  r[o] <- seq_along(o)

  # Calculate FDR on sorted data
  fdr_sorted <- cumsum(lfdr_sub[o]) / seq_along(o)
  fdr_out[rm_na] <- fdr_sorted[r]

  valid_fdr <- fdr_out[rm_na]
  max_fdr <- max(valid_fdr)
  n_sub <- length(r)

  p_calculated <- (valid_fdr * (r / n_sub)) / max_fdr
  p_out[rm_na] <- pmin(p_calculated, 1)

  out <- list(fp = p_out, fq = fdr_out)
  return(out)
}
