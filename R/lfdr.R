# Implements standard Q-value using local FDR for GWAS data
# Internal use at the moment, eventually will be public.
gwasQvalue <- function(p,
                       indep_snps = NULL,
                       trunc = TRUE,
                       monotone = TRUE,
                       nn = NULL,
                       pi0 = NULL,
                       transf = c("probit", "logit"),
                       adj = 0.5,
                       eps = 1e-15,
                       trim = 1e-15, ...) {
  # This code extends the qvalue::lfdr function to handle LD
  q_out <- lfdr_out <- p
  rm_na <- !is.na(p)
  if (is.null(indep_snps)) indep_snps <- rep(TRUE, length(p))
  indep_snps <- indep_snps[rm_na]
  p <- p[rm_na]
  # rm_na <- !is.na(p)
  p <- p
  if (min(p) < 0 || max(p) > 1) {
    stop("P-values not in valid range [0,1].")
  } else if (is.null(pi0)) {
    pi0 <- pi0est(p[indep_snps], ...)$pi0
  }
  # NN for density estimator
  if (is.null(nn)) {
    pi1 <- 1 - pi0
    nn <- min(pi1, 1 - pi1)
    if (nn < 0.02) nn <- 0.02
  }

  n <- length(p)
  transf <- match.arg(transf)
  if (transf == "probit") {
    p <- pmax(p, eps)
    p <- pmin(p, 1 - eps)
    kd <- kernelEstimator(p[indep_snps],
                          nn = nn,
                          trim = trim,
                          epsilon = eps,
                          eval.points = p)

    lfdr <- pi0  / kd$fx
  } else {
    x <- log((p + eps)/(1 - p + eps))
    myd <- density(x[indep_snps], adjust = adj)
    mys <- smooth.spline(x = myd$x, y = myd$y)
    y <- predict(mys, x)$y
    dx <- exp(x)/(1 + exp(x))^2
    lfdr <- (pi0 * dx) / y
  }
  if (trunc) {
    lfdr[lfdr > 1] <- 1
  }
  if (monotone) {
    o <- order(p, decreasing = FALSE)
    ro <- order(o)
    lfdr <- cummax(lfdr[o])[ro]
  }
  m <- length(p)
  i <- m:1L
  o <- order(p, decreasing = FALSE)
  ro <- order(o)
  lfdr[lfdr > 1] <- 1
  fdr <- cumsum(lfdr[o]) / (1:length(lfdr))
  fdr <- fdr[ro]
  q_out[rm_na] <- fdr
  lfdr_out[rm_na] <- lfdr
  return(list(qvalues = q_out,
              lfdr = lfdr_out,
              pi0 = pi0))
}
