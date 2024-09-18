fpvalues_raw <- function(lfdr) {
  p_out <- rep(NA, length(lfdr))
  rm_na <- !is.na(lfdr)

  lfdr <- lfdr[rm_na]
  lfdr[lfdr > 1] <- 1
  r <- rank(lfdr)

  emp_cdf <- ecdf(r / length(r))

  # correct ordering and calculate q-values
  out <- sort(lfdr, index.return = TRUE)
  out$x[out$x > 1] <- 1
  fdr <- cumsum(out$x) / (1:length(lfdr))
  fdr <- fdr[order(out$ix)]
  p_out[rm_na] <- pmin((fdr * emp_cdf(r / length(r))) / max(fdr), 1)
  return(list(fp = p_out,
              fq = fdr))
}

fpvalues <- function(lfdr, p = NULL) {
  if (is.null(p)) return(fpvalues_raw(lfdr)) # ignore order of p
  p_out <- fdr_out <- rep(NA, length(p))
  rm_na <- !is.na(p)
  p <- p[rm_na]
  lfdr <- lfdr[rm_na]
  lfdr[lfdr > 1] <- 1

  o = order(lfdr, p)
  r <- order(o)
  emp_cdf <- ecdf(r / length(r))

  # correct ordering and calculate q-values
  fdr <- (cumsum(lfdr[o]) / (1:length(lfdr)))[r]
  fdr_out[rm_na] <- fdr
  p_out[rm_na] <-  pmin((fdr * emp_cdf(r / length(r))) / max(fdr), 1)
  return(list(fp = p_out,
         fq = fdr_out))
}
