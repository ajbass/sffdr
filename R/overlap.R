# Sample overlap: the overlap correlation between studies ---------------------

#' Estimate the overlap correlation between studies
#'
#' @description
#' Estimates the null correlation of z-scores between studies that share
#' samples, for every pair of \code{cbind(z, z_informative)}. For each pair it
#' uses the mean product of z-scores inside the box \eqn{|z_1|, |z_2| < t},
#' which is unaffected by signal present in only one of the two studies, and
#' inverts it through the exact truncated bivariate normal moment. Cross-trait
#' LD score regression intercepts are preferable for GWAS when available.
#'
#' @param z Numeric vector; the primary study's z-scores.
#' @param z_informative Numeric vector or matrix of the informative studies' z-scores.
#' @param indep_snps Optional logical or integer vector of approximately independent SNPs to use.
#' @param threshold Numeric; the box half-width \eqn{t} (default 1.5). Smaller boxes
#'   are less affected by moderate effects and noisier.
#'
#' @return The \eqn{(K+1) \times (K+1)} correlation matrix, primary first.
#' @seealso \code{\link{decorrelate_informative}}
#' @export
estimate_overlap <- function(z, z_informative, indep_snps = NULL, threshold = 1.5) {
  Z <- cbind(z, as.matrix(z_informative))
  if (!is.numeric(Z)) stop("'z' and 'z_informative' must be numeric.")
  trn <- .train_set(indep_snps, nrow(Z))
  Z <- Z[trn & rowSums(!is.finite(Z)) == 0, , drop = FALSE]
  K1 <- ncol(Z)
  R <- diag(K1)
  for (i in seq_len(K1 - 1)) {
    for (j in (i + 1):K1) {
      R[i, j] <- R[j, i] <- .rho_box(Z[, i], Z[, j], threshold)
    }
  }
  ev <- eigen(R, symmetric = TRUE)                    # keep it a correlation matrix
  if (min(ev$values) < 1e-6) {
    R <- ev$vectors %*% diag(pmax(ev$values, 1e-6)) %*% t(ev$vectors)
    R <- R / sqrt(outer(diag(R), diag(R)))
  }
  dimnames(R) <- list(c("primary", paste0("informative", seq_len(K1 - 1))),
                      c("primary", paste0("informative", seq_len(K1 - 1))))
  R
}

## ---- internal helpers ----------------------------------------------------

## training set as a logical vector
.train_set <- function(indep_snps, n) {
  if (is.null(indep_snps)) return(rep(TRUE, n))
  if (is.logical(indep_snps)) {
    if (length(indep_snps) != n) stop("logical 'indep_snps' must have one entry per SNP.")
    return(!is.na(indep_snps) & indep_snps)
  }
  if (!is.numeric(indep_snps) || any(is.na(indep_snps)) || any(indep_snps < 1 | indep_snps > n)) {
    stop("'indep_snps' must be a logical vector or row indices between 1 and the number of SNPs.")
  }
  seq_len(n) %in% indep_snps
}

## E[w1 w2 | |w1| < t, |w2| < t] for a standard bivariate normal with correlation rho
.box_moment <- function(rho, t = 2, n = 401) {
  x <- seq(-t, t, length.out = n)
  h <- x[2] - x[1]
  wt <- c(1, rep(c(4, 2), length.out = n - 2), 1) * h / 3      # Simpson
  s <- sqrt(max(1 - rho^2, 1e-12))
  lo <- (-t - rho * x) / s
  hi <- (t - rho * x) / s
  pr <- pnorm(hi) - pnorm(lo)
  num <- sum(wt * dnorm(x) * x * (rho * x * pr + s * (dnorm(lo) - dnorm(hi))))
  den <- sum(wt * dnorm(x) * pr)
  num / den
}

## one pair: the box moment inverted through .box_moment
.rho_box <- function(z1, z2, t = 1.5) {
  k <- abs(z1) < t & abs(z2) < t
  if (sum(k) < 200) {
    stop("estimate_overlap: fewer than 200 SNPs with both |z| < ", t, "; supply the overlap correlation.")
  }
  m <- mean(z1[k] * z2[k])
  lim <- .box_moment(0.99, t)
  if (m >= lim) return(0.99)
  if (m <= -lim) return(-0.99)
  stats::uniroot(function(r) .box_moment(r, t) - m, c(-0.99, 0.99), tol = 1e-6)$root
}
