#' Discover Empirical Negative Controls
#'
#' Uses marginal local FDRs to dynamically identify genes that are highly
#' likely to be null in both traits, serving as empirical negative controls
#' for batch effect/overlap correlation estimation.
#'
#' @param p1 Numeric vector; primary trait p-values.
#' @param p2 Numeric vector; surrogate trait p-values.
#' @param threshold Numeric; the lfdr threshold for confidence (default 0.95).
#' @param min_genes Integer; the minimum number of nulls required for a stable
#'   covariance estimate.
#'
#' @return Integer vector of indices representing empirical null genes.
#' @export
discover_empirical_nulls <- function(
  p1,
  p2,
  threshold = 0.95,
  min_genes = 500
) {
  # 1. Compute marginal local FDRs safely
  lfdr1 <- tryCatch(
    qvalue::lfdr(
      p1,
      eps = 1e-25,
      monotone = TRUE,
      transf = "probit",
      trunc = TRUE
    ),
    error = function(e) rep(1, length(p1))
  )
  lfdr2 <- tryCatch(
    qvalue::lfdr(
      p2,
      eps = 1e-25,
      monotone = TRUE,
      transf = "probit",
      trunc = TRUE
    ),
    error = function(e) rep(1, length(p2))
  )

  # 2. Intersect the high-confidence nulls
  empirical_null_idx <- which(lfdr1 > threshold & lfdr2 > threshold)

  # 3. Dynamic Fallback
  if (length(empirical_null_idx) < min_genes) {
    warning(sprintf(
      "Only %d double-nulls found at strict threshold %.2f. Relaxing to top %d most 'null-like' genes.",
      length(empirical_null_idx),
      threshold,
      min_genes
    ))
    combined_lfdr <- lfdr1 * lfdr2
    empirical_null_idx <- order(combined_lfdr, decreasing = TRUE)[seq_len(
      min_genes
    )]
  }

  return(empirical_null_idx)
}

#' Run Overlap Correction using Data-Driven Nulls
#'
#' @param z1 Numeric vector; primary trait z-scores.
#' @param z2 Numeric vector; surrogate trait z-scores.
#' @param p1 Numeric vector; primary trait p-values.
#' @param p2 Numeric vector; surrogate trait p-values.
#'
#' @return A list containing the estimated rho_o, the conditional f0 density,
#'   and the indices of the empirical nulls used.
#' @export
run_empirical_overlap_correction <- function(z1, z2, p1, p2) {
  cat("Discovering empirical double-null genes...\n")
  empirical_nulls <- discover_empirical_nulls(p1, p2)
  cat(sprintf(
    "Isolated %d empirical negative controls.\n",
    length(empirical_nulls)
  ))

  rho_hat_obj <- estimate_rho_mom(
    z1 = z1,
    z2 = z2,
    p1 = p1,
    p2 = p2,
    neg_controls = empirical_nulls
  )
  safe_rho_o <- rho_hat_obj$rho_hat
  cat(sprintf("Estimated batch/overlap rho_o: %.3f\n", safe_rho_o))

  lfdry <- tryCatch(
    qvalue::lfdr(
      p2,
      eps = 1e-25,
      monotone = TRUE,
      transf = "probit",
      trunc = TRUE
    ),
    error = function(e) rep(1, length(p2))
  )

  cat("Computing conditional null density...\n")
  cond_f0 <- overlap_null_density(
    z1 = z1,
    z2 = z2,
    lfdry = lfdry,
    rho_o = safe_rho_o,
    empirical_null_var = TRUE
  )

  list(
    rho_o = safe_rho_o,
    f0 = cond_f0,
    empirical_null_indices = empirical_nulls
  )
}

#' Soft whitening of informative trait z-scores
#'
#' Removes overlap-induced correlation from informative trait z-scores (`z_y`)
#' without altering primary p-values. The correction is scaled by the
#' primary trait's local false discovery rate (lfdr) to preserve true
#' biological signal while removing shared noise.
#'
#' @param z_x      Numeric vector; primary trait z-scores (used for weighting).
#' @param z_y      Numeric vector; informative trait z-scores (transformed).
#' @param p_x      Numeric vector; primary trait p-values.
#' @param p_y      Numeric vector; informative trait p-values.
#' @param rho_o    Numeric scalar; sample overlap correlation. If NULL, estimated
#'                 via \code{estimate_rho_overlap}. Supplying LDSC intercepts is recommended.
#' @param lfdr_x   Numeric vector; pre-computed marginal lfdr for primary trait.
#'                 If NULL, computed internally via \code{qvalue::lfdr}.
#'
#' @return A list containing:
#' \describe{
#'   \item{z2}{Numeric vector; decorrelated informative trait z-scores.}
#'   \item{p2}{Numeric vector; decorrelated informative trait p-values.}
#'   \item{lfdr_x}{Numeric vector; marginal lfdr values used for weighting.}
#'   \item{rho_o}{Numeric scalar; overlap correlation used.}
#' }
#'
#' @examples
#' \dontrun{
#' dec    <- decorrelate_informative(z_x, z_y, p_x, p_y)
#' mpi0   <- pi0_model(as.matrix(dec$p2))
#' fpi0   <- fpi0est(p_x, pi0_model_obj = mpi0)
#' result <- sffdr(p_x, fpi0 = fpi0$fpi0, surrogate = dec$p2)
#' }
#'
#' @seealso \code{\link{estimate_rho_overlap}}, \code{\link{sffdr}}
#' @export
decorrelate_informative <- function(
  z_x,
  z_y,
  p_x,
  p_y,
  rho_o = NULL,
  lfdr_x = NULL
) {
  # --- Input validation ---
  n <- length(z_x)
  if (!is.numeric(z_x)) {
    stop("'z_x' must be a numeric vector.")
  }
  if (!is.numeric(z_y)) {
    stop("'z_y' must be a numeric vector.")
  }
  if (!is.numeric(p_x)) {
    stop("'p_x' must be a numeric vector.")
  }
  if (!is.numeric(p_y)) {
    stop("'p_y' must be a numeric vector.")
  }

  if (length(z_y) != n || length(p_x) != n || length(p_y) != n) {
    stop(
      "Input vectors 'z_x', 'z_y', 'p_x', and 'p_y' must be the same length."
    )
  }
  if (
    any(p_x < 0 | p_x > 1, na.rm = TRUE) || any(p_y < 0 | p_y > 1, na.rm = TRUE)
  ) {
    stop("P-values must be in the range [0, 1].")
  }

  if (!is.null(rho_o)) {
    if (!is.numeric(rho_o) || length(rho_o) != 1) {
      stop("'rho_o' must be a numeric scalar.")
    }
    if (rho_o < -1 || rho_o > 1) stop("'rho_o' must be in [-1, 1].")
  }

  if (!is.null(lfdr_x) && (!is.numeric(lfdr_x) || length(lfdr_x) != n)) {
    stop("'lfdr_x' must be a numeric vector of length equal to 'z_x'.")
  }

  # --- NA handling ---
  valid <- !is.na(z_x) & !is.na(z_y) & !is.na(p_x) & !is.na(p_y)
  n_na <- sum(!valid)
  if (n_na > 0) {
    warning(sprintf(
      "decorrelate_informative: %d observation(s) with NA excluded.",
      n_na
    ))
  }

  # --- Estimate rho_o if not supplied ---
  if (is.null(rho_o)) {
    rho_o <- estimate_rho_overlap(
      z_x[valid],
      z_y[valid],
      p_x[valid],
      p_y[valid],
      method = "double_null"
    )
  }

  if (abs(rho_o) > 0.3) {
    warning(sprintf(
      "|rho_o| = %.3f is large. Soft whitening may be inaccurate. Supply LDSC intercept if available.",
      abs(rho_o)
    ))
  }

  # --- Marginal lfdr for primary trait ---
  if (is.null(lfdr_x)) {
    lfdr_x <- qvalue::lfdr(
      p_x,
      eps = 1e-25,
      monotone = TRUE,
      transf = "probit",
      trunc = TRUE
    )
  }
  lfdr_x <- pmin(pmax(lfdr_x, 0), 1)

  # --- Soft whitening ---
  z_y_star <- z_y - rho_o * z_x * lfdr_x
  p_y_star <- 2 * pnorm(-abs(z_y_star))

  # Map NA inputs back to NA in output
  if (n_na > 0) {
    z_y_star[!valid] <- NA_real_
    p_y_star[!valid] <- NA_real_
    lfdr_x[!valid] <- NA_real_
  }

  list(z2 = z_y_star, p2 = p_y_star, lfdr_x = lfdr_x, rho_o = rho_o)
}

#' Estimate sample overlap correlation
#'
#' Computes the empirical null correlation between two traits due to sample
#' overlap. Whenever possible, supplying an external \code{rho_o} (e.g., from
#' LDSC intercepts) is recommended over empirical estimation.
#'
#' @param z1,z2  Numeric vectors; z-scores for traits 1 and 2.
#' @param p1,p2  Numeric vectors; p-values for traits 1 and 2.
#' @param lfdry  Numeric vector; marginal lfdr for trait 2 (required only
#'               for \code{method = "weighted_cov"}).
#' @param method Character; the estimation strategy to use:
#'               \itemize{
#'                 \item \code{"weighted_cov"}: Weighted least squares via lfdr.
#'                 \item \code{"double_null"}: Pearson correlation strictly on the double-null stratum.
#'                 \item \code{"mom_scaled"}: Method of moments (scaled lambdas).
#'                 \item \code{"mom_rank"}: Method of moments (rank-normalised).
#'               }
#' @param thresh Numeric scalar; lfdr threshold for defining the double-null
#'               stratum in the \code{"double_null"} method. Default is 0.5.
#'
#' @return Numeric scalar representing the estimated overlap correlation.
#' @export
estimate_rho_overlap <- function(
  z1,
  z2,
  p1,
  p2,
  lfdry = NULL,
  method = "double_null",
  thresh = 0.5
) {
  method <- match.arg(
    method,
    c("weighted_cov", "double_null", "mom_scaled", "mom_rank")
  )

  if (!is.numeric(z1)) {
    stop("'z1' must be a numeric vector.")
  }
  if (!is.numeric(z2)) {
    stop("'z2' must be a numeric vector.")
  }
  if (!is.numeric(p1)) {
    stop("'p1' must be a numeric vector.")
  }
  if (!is.numeric(p2)) {
    stop("'p2' must be a numeric vector.")
  }

  n <- length(z1)
  if (length(z2) != n || length(p1) != n || length(p2) != n) {
    stop("Input vectors 'z1', 'z2', 'p1', and 'p2' must be the same length.")
  }

  valid <- !is.na(z1) & !is.na(z2) & !is.na(p1) & !is.na(p2)
  if (!is.null(lfdry)) {
    valid <- valid & !is.na(lfdry)
  }

  n_na <- sum(!valid)
  if (n_na > 0) {
    warning(sprintf(
      "estimate_rho_overlap: %d NA observation(s) excluded.",
      n_na
    ))
  }

  z1 <- z1[valid]
  z2 <- z2[valid]
  p1 <- p1[valid]
  p2 <- p2[valid]
  if (!is.null(lfdry)) {
    lfdry <- lfdry[valid]
  }

  if (length(z1) < 2) {
    return(NA_real_)
  }

  # --- Estimation ---
  if (method == "weighted_cov") {
    if (is.null(lfdry)) {
      stop("'weighted_cov' method requires 'lfdry'.")
    }
    w <- lfdry
    mx <- sum(w * z2) / sum(w)
    my <- sum(w * z1) / sum(w)
    w_cov <- sum(w * (z2 - mx) * (z1 - my))
    w_var <- sum(w * (z2 - mx)^2)
    rho_o <- if (w_var > 1e-8) w_cov / w_var else 0
    rho_o <- max(min(rho_o, 0.99), -0.99)
  } else if (method == "double_null") {
    lfdr1 <- tryCatch(
      qvalue::lfdr(
        p1,
        eps = 1e-25,
        monotone = TRUE,
        transf = "probit",
        trunc = TRUE
      ),
      error = function(e) rep(1, length(p1))
    )
    lfdr2 <- tryCatch(
      qvalue::lfdr(
        p2,
        eps = 1e-25,
        monotone = TRUE,
        transf = "probit",
        trunc = TRUE
      ),
      error = function(e) rep(1, length(p2))
    )

    lfdr1 <- pmin(pmax(lfdr1, 0), 1)
    lfdr2 <- pmin(pmax(lfdr2, 0), 1)

    # The > threshold filter inherently removes the heavy tails
    null_idx <- which(lfdr1 > thresh & lfdr2 > thresh)

    if (length(null_idx) < 1000 && length(z1) >= 1000) {
      warning(
        "Fewer than 1,000 double-null genes found. Estimate may be imprecise."
      )
    }

    if (length(null_idx) < 2) {
      warning("Too few null genes for correlation. Returning 0.")
      rho_o <- 0
    } else {
      # Use Pearson directly on the well-behaved null core
      rho_o <- cor(
        z1[null_idx],
        z2[null_idx],
        method = "pearson",
        use = "complete.obs"
      )
      rho_o <- max(min(if (is.na(rho_o)) 0 else rho_o, 0.99), -0.99)
    }
  } else {
    pi0_method <- if (method == "mom_scaled") {
      "2d_storey_scaled"
    } else {
      "2d_storey_rank"
    }
    rho_o <- estimate_rho_mom(z1, z2, p1, p2, pi0_method = pi0_method)$rho_hat
  }

  return(rho_o)
}

#' Estimate overlap correlation via method of moments
#'
#' @param z1,z2        Numeric vectors; z-scores.
#' @param p1,p2        Numeric vectors; p-values.
#' @param neg_controls Integer vector; indices of known null genes.
#' @param pi0_method   Character; method for joint pi0 estimation.
#' @param ...          Additional arguments passed to \code{estimate_pi0_joint}.
#' @noRd
estimate_rho_mom <- function(
  z1,
  z2,
  p1,
  p2,
  neg_controls = NULL,
  pi0_method = "2d_storey_scaled",
  ...
) {
  if (!is.null(neg_controls)) {
    rho_hat <- cor(
      z1[neg_controls],
      z2[neg_controls],
      method = "pearson",
      use = "complete.obs"
    )
    return(list(
      rho_hat = max(min(rho_hat, 0.99), -0.99),
      method = "neg_controls",
      pi0_joint = NA
    ))
  }

  if (pi0_method == "2d_storey_scaled") {
    pi0_obj <- estimate_pi0_joint(p1, p2, method = "scaled", ...)
  } else if (pi0_method == "2d_storey_rank") {
    pi0_obj <- estimate_pi0_joint(p1, p2, method = "rank", ...)
  } else if (pi0_method == "product") {
    pi0_obj <- list(
      pi0_joint = qvalue::qvalue(p1)$pi0 * qvalue::qvalue(p2)$pi0,
      method = "product"
    )
  } else {
    stop(
      "pi0_method must be '2d_storey_scaled', '2d_storey_rank', or 'product'"
    )
  }

  rho_hat <- max(min(mean(z1 * z2) / pi0_obj$pi0_joint, 0.99), -0.99)

  list(
    rho_hat = rho_hat,
    pi0_joint = pi0_obj$pi0_joint,
    raw_moment = mean(z1 * z2),
    pi0_method = pi0_method,
    pi0_obj = pi0_obj
  )
}

#' Estimate joint null proportion (pi0)
#'
#' @noRd
estimate_pi0_joint <- function(
  p1,
  p2,
  lambdas = seq(0.05, 0.8, by = 0.05),
  method = "scaled",
  spline_df = 3,
  min_n = 50,
  plot = FALSE
) {
  M <- length(p1)

  if (method == "rank") {
    p1_use <- rank(p1) / M
    p2_use <- rank(p2) / M
    lam1_seq <- lambdas
    lam2_seq <- lambdas
  } else if (method == "scaled") {
    p1_use <- p1
    p2_use <- p2
    pi0_1 <- qvalue::qvalue(p1)$pi0
    pi0_2 <- qvalue::qvalue(p2)$pi0
    lam1_seq <- pmin(lambdas + (1 - lambdas) * (1 - pi0_1), 0.95)
    lam2_seq <- pmin(lambdas + (1 - lambdas) * (1 - pi0_2), 0.95)
  }

  pi0_seq <- mapply(
    function(lam1, lam2) {
      idx <- sum(p1_use > lam1 & p2_use > lam2)
      if (idx < min_n) {
        return(NA_real_)
      }
      idx / (M * (1 - lam1) * (1 - lam2))
    },
    lam1_seq,
    lam2_seq
  )

  keep <- !is.na(pi0_seq) & is.finite(pi0_seq)
  if (sum(keep) < 4) {
    stop("Too few valid lambda values - reduce min_n or widen lambda range")
  }

  lam_use <- lambdas[keep]
  pi0_use <- pi0_seq[keep]

  spl <- smooth.spline(lam_use, pi0_use, df = spline_df)

  upper <- lam_use > median(lam_use)
  if (sum(upper) == 0) {
    pi0_est <- min(predict(spl, lam_use)$y)
  } else {
    pi0_est <- min(predict(spl, lam_use[upper])$y)
  }
  pi0_est <- max(min(pi0_est, 1.0), 0.0)

  list(
    pi0_joint = pi0_est,
    pi0_seq = pi0_seq,
    pi0_fit = predict(spl, lam_use)$y,
    lambdas = lambdas,
    lam1_seq = lam1_seq,
    lam2_seq = lam2_seq,
    method = method
  )
}


#' Compute conditional null density under sample overlap
#'
#' Computes the conditional null density f(z1 | H1=0, z2) for a primary
#' trait (z1) given a surrogate trait (z2) when the two studies share samples.
#'
#' @param z1 Numeric vector; z-scores for the primary trait.
#' @param z2 Numeric vector; z-scores for the surrogate trait.
#' @param lfdry Numeric vector; marginal local FDR of the surrogate trait.
#' @param rho_o Numeric scalar; expected sample overlap correlation. If \code{NULL}
#'   (default), it is estimated internally via \code{estimate_rho_overlap}.
#' @param p1 Numeric vector; primary trait p-values. Required if \code{rho_o = NULL}.
#' @param p2 Numeric vector; surrogate trait p-values. Required if \code{rho_o = NULL}.
#' @param se2 Numeric vector; standard errors of the surrogate trait effects.
#'   If provided, enables per-variant Empirical Bayes shrinkage (recommended for GWAS).
#'   If NULL, uses global shrinkage (recommended for RNA-seq).
#' @param empirical_null_var Logical; if TRUE (default), estimates null variance
#'   empirically from the double-null center to absorb genomic inflation.
#' @param ... Additional arguments passed to \code{estimate_rho_overlap}.
#'
#' @return A numeric vector representing the conditional null density for each
#'   observation. Pass this to the \code{f0} argument of \code{sffdr}.
#'
#' @export
overlap_null_density <- function(
  z1,
  z2,
  lfdry,
  rho_o = NULL,
  p1 = NULL,
  p2 = NULL,
  se2 = NULL,
  empirical_null_var = TRUE,
  ...
) {
  if (is.null(rho_o)) {
    if (is.null(p1) || is.null(p2)) {
      stop(
        "If `rho_o` is NULL, `p1` and `p2` must be provided to estimate it dynamically."
      )
    }
    rho_o <- estimate_rho_overlap(
      z1,
      z2,
      p1,
      p2,
      lfdry = lfdry,
      method = "double_null",
      ...
    )
  }

  if (!is.null(se2)) {
    W_hat <- max(mean((z2^2 - 1) * se2^2, na.rm = TRUE), 1e-4)
    r <- W_hat / (W_hat + se2^2)
  } else {
    sigma_s2_z <- max(mean(z2^2, na.rm = TRUE) - 1, 0)
    r <- sigma_s2_z / (sigma_s2_z + 1)
  }

  w_sum <- sum(lfdry, na.rm = TRUE)
  mx <- sum(lfdry * z2, na.rm = TRUE) / w_sum
  my <- sum(lfdry * z1, na.rm = TRUE) / w_sum

  z2_centered <- z2 - mx
  mean_null <- rho_o * z2_centered + my
  mean_signal <- rho_o * (1 - r) * z2_centered + my

  var_theoretical <- 1 - rho_o^2

  if (empirical_null_var) {
    resid <- z1 - mean_null
    var_null <- if (length(resid) > 30) {
      (mad(resid, na.rm = TRUE))^2
    } else {
      var_theoretical
    }
    var_null <- max(var_null, var_theoretical, 0.1)
  } else {
    var_null <- var_theoretical
  }

  var_signal <- pmax(1 - rho_o^2 * (1 - r), 0.1)

  mean_mix <- lfdry * mean_null + (1 - lfdry) * mean_signal
  var_within <- lfdry * var_null + (1 - lfdry) * var_signal
  var_between <- lfdry * (1 - lfdry) * (mean_null - mean_signal)^2
  var_mix <- pmax(var_within + var_between, 1e-6)

  dnorm(z1, mean = mean_mix, sd = sqrt(var_mix))
}
