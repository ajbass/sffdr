#' Kernel Density Estimation for GWAS P-values
#'
#' Performs kernel density estimation on p-values (univariate) or joint
#' p-value/covariate pairs (bivariate) using local regression on the probit scale.
#' The estimator is optimized for GWAS data with linkage disequilibrium (LD) and
#' uses adaptive downsampling to prioritize signal-rich regions while maintaining
#' computational efficiency.
#'
#' @param x Numeric vector of p-values (for 1D density) or a 2-column matrix where
#'   the first column contains p-values and the second column contains an informative
#'   covariate/surrogate. All p-values must be in (0, 1) and the covariate/surrogate must
#'   be rank-transformed to be (0,1].
#'
#' @param eval.points Points at which to evaluate the density estimate. Defaults to
#'   \code{x}. For custom evaluation points, must match the dimensionality of \code{x}
#'   (vector or 2-column matrix).
#'
#' @param epsilon Lower bound for p-values to prevent numerical issues. P-values
#'   below this are clamped to \code{epsilon}. Default: \code{1e-15}.
#'
#' @param epsilon.max Upper bound for p-values. P-values above this are clamped to
#'   \code{epsilon.max}. Default: \code{1 - 1e-4}.
#'
#' @param maxk Maximum number of fitting points passed to \code{\link[locfit]{locfit}}.
#'   Increase for very large datasets. Default: \code{500000}.
#'
#' @param maxit Maximum number of iterations for local regression fitting. Default:
#'   \code{2000}.
#'
#' @param target_null Maximum number of null SNPs to include in the weighted fit
#'   (bivariate case only). SNPs in the signal-enriched tail (defined by
#'   \code{tail_threshold}) are always retained; null SNPs are downsampled to this
#'   target and upweighted accordingly. Default: \code{100000}.
#'
#' @param trim For 1D estimation, fixes the density estimate to constant values on
#'   the intervals \code{(0, trim)} and \code{(1 - trim, 1)} to reduce boundary
#'   variance. Default: \code{0} (no trimming).
#'
#' @param nn Nearest-neighbor bandwidth parameter for \code{\link[locfit]{locfit}},
#'   expressed as a fraction of the data. If \code{NULL} (default), automatically
#'   determined based on effective sample size to span approximately 5000 neighbors
#'   (which corresponds to multiple LD blocks in GWAS). Larger values increase smoothing.
#'
#' @param tail_threshold Z-score threshold on the probit-transformed covariate scale
#'   (bivariate case only). SNPs with z < \code{tail_threshold} are treated as
#'   signal-enriched and prioritized in the adaptive fit. Default: \code{-2}
#'   (approximately 2.3\% of standard normal distribution).
#'
#' @param weights Optional numeric vector of weights for density estimation. For GWAS data, this should be inverse LD scores.
#'
#' @param ... Additional arguments passed to \code{\link[locfit]{lp}} for controlling
#'   local polynomial fitting (e.g., \code{deg}, \code{kern}).
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{x}{Evaluation points (original scale).}
#'     \item{fx}{Estimated density at evaluation points (original scale).}
#'     \item{s}{Evaluation points on probit scale (\code{qnorm(x)}).}
#'     \item{fs}{Estimated density on probit scale.}
#'   }
#'   The returned object has an attribute \code{"lfit"} containing the fitted
#'   \code{locfit} object for diagnostics.
#'
#' @details
#' The function implements a multi-stage density estimation procedure:
#'
#' \enumerate{
#'   \item \strong{Probit transformation}: P-values are transformed to the normal
#'     scale via \code{qnorm} to stabilize variance and handle extreme values.
#'
#'   \item \strong{Adaptive downsampling (bivariate only)}: To handle large GWAS
#'     datasets efficiently, the null region (where the covariate suggests low signal)
#'     is downsampled to \code{target_null} SNPs, with inverse-probability weighting
#'     to preserve the density. Signal-enriched SNPs (tail) are always retained.
#'
#'   \item \strong{Cascade fitting}: Multiple fitting strategies are attempted in
#'     sequence, with decreasing resolution and increasing robustness, until a
#'     valid fit is obtained.
#'
#'   \item \strong{Jacobian correction}: Density estimates are transformed back to
#'     the original p-value scale using the Jacobian of the probit transformation.
#' }
#'
#' The nearest-neighbor bandwidth \code{nn} controls smoothing and LD robustness.
#' By targeting ~5000 neighbors (default), the estimator naturally averages over
#' multiple LD blocks (~30-50 blocks in European ancestry populations), reducing
#' spurious local structure while preserving true signal-covariate relationships.
#'
#' @seealso \code{\link[locfit]{locfit}}, \code{\link{sffdr}}
#'
#' @examples
#' \dontrun{
#' # 1D density estimation
#' p <- runif(10000, 0, 1)
#' dens <- kernelEstimator(p)
#' plot(dens$x, dens$fx, type = "l")
#'
#' # 2D density with informative covariate
#' p <- runif(10000, 0, 1)
#' z <- runif(10000)  # rank-norm transformed covariate
#' x_mat <- cbind(p, z)
#' dens <- kernelEstimator(x_mat)
#' }
#'
#' @export
kernelEstimator <- function(
  x,
  eval.points = x,
  epsilon = 1e-15,
  epsilon.max = 1 - 1e-4,
  maxk = 500000,
  maxit = 200,
  target_null = 100000,
  trim = 0,
  nn = NULL,
  tail_threshold = -2,
  weights = NULL,
  ...
) {
  # Input validation
  is_matrix <- is.matrix(x)
  if (is_matrix && ncol(x) < 2) {
    stop("Matrix input must have at least 2 columns.")
  }

  # Transform to probit scale with clamping
  clamp_and_transform <- function(vals) {
    qnorm(pmin(pmax(vals, epsilon), epsilon.max))
  }

  train_s <- clamp_and_transform(x)
  eval_s <- if (identical(x, eval.points)) {
    train_s
  } else {
    clamp_and_transform(eval.points)
  }

  # Add jitter to break ties and improve convergence
  jitter_mag <- 1e-6
  if (is_matrix) {
    train_s[, 2] <- train_s[, 2] + runif(nrow(train_s), -jitter_mag, jitter_mag)
  } else {
    train_s <- train_s + runif(length(train_s), -jitter_mag, jitter_mag)
  }

  # Fit density model
  if (is_matrix) {
    lfit <- fit_bivariate_density(
      train_s = train_s,
      tail_threshold = tail_threshold,
      target_null = target_null,
      nn = nn,
      maxit = maxit,
      maxk = maxk,
      weights = weights,
      ...
    )
  } else {
    lfit <- fit_univariate_density(
      train_s = train_s,
      nn = nn,
      maxit = maxit,
      maxk = maxk,
      weights = weights,
      ...
    )
  }

  if (is.null(lfit)) {
    stop("Kernel estimation failed (all strategies exhausted).")
  }

  # Predict density on probit scale
  newdata <- if (is.matrix(eval_s)) {
    data.frame(V1 = eval_s[, 1], V2 = eval_s[, 2])
  } else {
    data.frame(train_s = eval_s)
  }
  fs_hat <- predict(lfit, newdata = newdata, log = FALSE)

  # Jacobian correction to transform back to original scale
  corrector <- if (is.matrix(eval_s)) {
    dnorm(eval_s[, 1]) * dnorm(eval_s[, 2])
  } else {
    dnorm(eval_s)
  }
  fx_hat <- fs_hat / pmax(corrector, .Machine$double.xmin)

  # Build result data frame
  res <- data.frame(
    x = eval.points,
    fx = fx_hat,
    s = eval_s,
    fs = fs_hat
  )

  # Apply boundary trimming for 1D case
  if (trim > 0 && !is_matrix) {
    val_lower <- res$fx[which.min(abs(res$x - trim))]
    val_upper <- res$fx[which.min(abs(res$x - (1 - trim)))]
    res$fx[res$x < trim] <- val_lower
    res$fx[res$x > (1 - trim)] <- val_upper
  }

  attr(res, "lfit") <- lfit
  res
}


#' Fit Bivariate Density with Adaptive Downsampling
#'
#' Internal helper for 2D kernel density estimation with signal prioritization.
#'
#' @param train_s 2-column matrix on probit scale (p-values, covariate).
#' @param tail_threshold Z-score threshold for signal detection.
#' @param target_null Maximum null SNPs to retain after downsampling.
#' @param nn Nearest-neighbor bandwidth (NULL for automatic).
#' @param maxit Maximum iterations for locfit.
#' @param maxk Maximum fitting points for locfit.
#' @param ... Additional arguments for lp().
#'
#' @return Fitted locfit object or NULL.
#' @keywords internal
fit_bivariate_density <- function(
  train_s,
  tail_threshold,
  target_null,
  nn,
  maxit,
  maxk,
  weights = NULL,
  ...
) {
  has_custom_weights <- !is.null(weights) && !all(weights == 1.0)

  if (is.null(weights)) {
    weights <- rep(1.0, nrow(train_s))
  }

  z <- train_s[, 2]
  idx_tail <- which(z < tail_threshold)
  n_tail <- length(idx_tail)

  if (n_tail == 0) {
    stop("No signal points found. Consider adjusting tail_threshold.")
  }

  n_total <- nrow(train_s)
  n_null <- n_total - n_tail

  if (n_null > target_null) {
    idx_null <- which(z >= tail_threshold)
    idx_keep <- sample(idx_null, target_null)
    idx_final <- c(idx_tail, idx_keep)
    w_null <- n_null / target_null
    w_vec <- c(weights[idx_tail], weights[idx_keep] * w_null)
  } else {
    idx_final <- seq_len(n_total)
    w_vec <- weights
  }

  fit_data <- data.frame(
    V1 = train_s[idx_final, 1],
    V2 = train_s[idx_final, 2],
    w = w_vec
  )

  if (has_custom_weights) {
    fit_data$V1_round <- round(fit_data$V1, 3)
    fit_data$V2_round <- round(fit_data$V2, 3)
    fit_data <- aggregate(w ~ V1_round + V2_round, data = fit_data, FUN = sum)
    names(fit_data) <- c("V1", "V2", "w")
  }

  valid_data <- fit_data[fit_data$V2 < tail_threshold, , drop = FALSE]

  n_eff <- nrow(fit_data)
  if (!is.null(nn)) {
    nn_high <- nn_safe <- nn
  } else {
    nn_high <- max(0.005, min(0.1, 5000 / n_eff))
    nn_safe <- max(0.010, min(0.1, 5000 / n_eff))
  }

  fit_strategy_final(
    data = fit_data,
    valid_data = valid_data,
    nn_high = nn_high,
    nn_safe = nn_safe,
    maxit = maxit,
    maxk = maxk,
    ...
  )
}


#' Fit Univariate Density
#'
#' Internal helper for 1D kernel density estimation.
#'
#' @param train_s Numeric vector on probit scale.
#' @param nn Nearest-neighbor bandwidth (NULL for automatic).
#' @param maxit Maximum iterations for locfit.
#' @param maxk Maximum fitting points for locfit.
#' @param ... Additional arguments for lp().
#'
#' @return Fitted locfit object.
#' @keywords internal
fit_univariate_density <- function(
  train_s,
  nn,
  maxit,
  maxk,
  weights = NULL,
  ...
) {
  n_eff <- length(train_s)
  has_custom_weights <- !is.null(weights) && !all(weights == 1.0)

  if (is.null(weights)) {
    weights <- rep(1.0, n_eff)
  }

  fit_data <- data.frame(train_s = train_s, w = weights)

  if (has_custom_weights) {
    fit_data$train_s_round <- round(fit_data$train_s, 3)
    fit_data <- aggregate(w ~ train_s_round, data = fit_data, FUN = sum)
    names(fit_data) <- c("train_s", "w")
    n_eff <- nrow(fit_data)
  }

  nn_base <- if (!is.null(nn)) {
    nn
  } else {
    max(0.01, min(0.1, 10000 / n_eff))
  }

  locfit(
    ~ lp(train_s, nn = nn_base, h = 0.05, ...),
    data = fit_data,
    weights = fit_data$w,
    maxk = maxk,
    maxit = maxit
  )
}

#' Fit locfit with Multiple Fallback Strategies
#'
#' Internal helper function that attempts to fit a locfit model using multiple
#' strategies with different kernels, degrees, and smoothing parameters.
#'
#' @param data Data frame with columns V1, V2, and w (weights).
#' @param valid_data Subset of data used for validation (typically signal-rich region).
#' @param nn_high First (higher resolution) nn value to try.
#' @param nn_safe Second (safer, more conservative) nn value.
#' @param maxit Maximum iterations for locfit.
#' @param maxk Maximum number of fitting points for locfit.
#' @param verbose Logical; print progress messages.
#' @param ... Additional arguments passed to lp() in locfit.
#'
#' @return A fitted locfit object, or NULL if all strategies fail.
#' @keywords internal
fit_strategy_final <- function(
  data,
  valid_data,
  nn_high,
  nn_safe,
  maxit,
  maxk,
  verbose = TRUE,
  ...
) {
  # Subset for validation checks
  check_data <- if (nrow(valid_data) > 0) {
    valid_data
  } else {
    data[seq_len(min(100, nrow(data))), ]
  }

  # Validation function: check predictions are valid
  validate <- function(mod) {
    if (is.null(mod)) {
      return(FALSE)
    }
    p_log <- try(predict(mod, newdata = check_data, log = TRUE), silent = TRUE)
    if (
      inherits(p_log, "try-error") || !all(is.finite(p_log)) || max(p_log) > 709
    ) {
      return(FALSE)
    }
    p_lin <- exp(p_log)
    all(p_lin > 0)
  }
  strats <- list(
    list(
      name = "High-Res Tcub",
      deg = 2,
      kern = "tcub",
      ev = rbox(),
      nn = nn_high,
      h = 0.05
    ),
    list(
      name = "High-Res Gauss",
      deg = 2,
      kern = "gauss",
      ev = rbox(),
      nn = nn_high,
      h = 0.05
    ),
    list(
      name = "Safe Tcub",
      deg = 2,
      kern = "tcub",
      ev = rbox(),
      nn = nn_safe,
      h = 0.1
    ),
    list(
      name = "Linear Tree",
      deg = 1,
      kern = "tcub",
      ev = rbox(),
      nn = 2 * nn_safe,
      h = 0.2
    )
  )

  # Try each strategy
  for (s in strats) {
    if (verbose) {
      message(
        "    ",
        s$name,
        " (nn=",
        round(s$nn, 5),
        ", h=",
        s$h,
        ")...",
        appendLF = FALSE
      )
    }

    fatal <- FALSE
    fit <- tryCatch(
      {
        withCallingHandlers(
          {
            locfit(
              ~ lp(V1, V2, nn = s$nn, h = s$h, deg = s$deg), # <-- h goes inside lp()
              data = data,
              weights = data$w,
              kern = s$kern,
              ev = s$ev,
              maxit = maxit,
              maxk = maxk,
              ...
            )
          },
          warning = function(w) {
            msg <- conditionMessage(w)
            if (
              grepl(
                "procv|vertex|newsplit|maxk|max_nr|convergence|singular|zero variance",
                msg,
                ignore.case = TRUE
              )
            ) {
              fatal <<- TRUE
            }
            invokeRestart("muffleWarning")
          }
        )
      },
      error = function(e) NULL
    )

    if (!fatal && validate(fit)) {
      if (verbose) {
        message(" \u2713")
      }
      return(fit)
    } else {
      if (verbose) {
        message(" \u2717")
      }
    }
  }

  NULL
}
