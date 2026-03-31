#' Estimate Functional Proportion of Null Tests
#'
#' @description
#' Estimates the functional proportion of null tests (pi0) using a GLM approach
#' with a constrained binomial family. This function fits models across multiple
#' lambda thresholds and selects the optimal estimate via MISE minimization.
#'
#' @details
#' \strong{Algorithm:}
#'
#' \enumerate{
#'   \item For each lambda threshold, fit a binomial GLM: \eqn{P(p \ge \lambda | z)}
#'   \item Use constrained binomial family to ensure predictions in (0, 1)
#'   \item Select optimal lambda via MISE minimization
#' }
#'
#' \strong{Usage Patterns:}
#'
#' \strong{Pattern 1 (Recommended):} Use output from \code{\link{pi0_model}}
#' \preformatted{
#'   mpi0 <- pi0_model(z)
#'   # Clean syntax:
#'   fpi0_out <- fpi0est(p, mpi0)
#' }
#'
#' \strong{Pattern 2:} Manually specify formula and covariates
#' \preformatted{
#'   fpi0_out <- fpi0est(p, z = z_matrix, pi0_model = formula_obj)
#' }
#'
#' @param p Numeric vector of p-values.
#' @param pi0_model_obj Optional list object returned by \code{\link{pi0_model}}.
#'   Must contain \code{fmod} (formula) and \code{zt} (rank-transformed covariates).
#'   If provided, \code{z} and \code{pi0_model} are extracted automatically.
#'   Default is NULL.
#' @param z Optional data frame or matrix of rank-transformed covariates.
#'   Required if \code{pi0_model_obj} is not provided. Default is NULL.
#' @param pi0_model Optional formula (as character string or formula object).
#'   Required if \code{pi0_model_obj} is not provided. Default is NULL.
#' @param indep_snps Optional logical vector indicating independent SNPs for
#'   model fitting. Default is NULL (all SNPs used).
#' @param weights Optional numeric vector of weights for density estimation. For GWAS data, this should be inverse LD scores. If these are available
#'   then you don't have to use the indep_snps argument, as the weights will effectively prioritize independent SNPs in the fitting process.
#' @param lambda Numeric vector of lambda thresholds. Default is \code{seq(0.05, 0.95, 0.05)}.
#' @param constrained.p Logical; use constrained binomial family. Default is TRUE.
#' @param tol Numeric; convergence tolerance. Default is 1e-9.
#' @param maxit Integer; maximum iterations. Default is 200.
#' @param ncores Integer; number of cores for parallel lambda fitting via
#'   \code{parallel::mclapply}. Fork-based parallelism (Unix/macOS only) -
#'   on Windows this falls back to sequential execution. Default is 1L.
#' @param verbose Logical; print progress messages. Default is TRUE.
#' @param ... Additional arguments passed to \code{\link[fastglm]{fastglm}}.
#'
#' @examples
#' # Import data
#' data(bmi)
#'
#' # Separate main p-values and conditioning p-values
#' p <- sumstats$bmi
#' z <- as.matrix(sumstats[, -1])
#'
#' # Apply pi0_model to create model (uses adaptive knot selection)
#' fmod <- pi0_model(z)
#'
#' # Estimate functional pi0
#' fpi0_out <- fpi0est(p, fmod)
#' fpi0 <- fpi0_out$fpi0
#'
#' # Apply sffdr
#' sffdr_out <- sffdr(p, fpi0)
#'
#' @return An object of class \code{fpi0} (a list) containing:
#' \describe{
#'   \item{fpi0}{Numeric vector of functional pi0 estimates for each test.}
#'   \item{tableLambda}{A data frame summarizing results for each lambda value.}
#'   \item{MISE}{The Mean Integrated Squared Error (MISE) for the chosen model.}
#'   \item{lambda}{The selected optimal lambda value.}
#' }
#'
#' @importFrom stats mad pnorm cor family predict model.matrix qnorm dnorm runif
#'   quantile na.pass optimize binomial dbinom formula fitted fitted.values
#'   gaussian median model.frame model.offset model.weights napredict
#'   delete.response terms .checkMFClasses .getXlevels dlogis plogis aggregate
#'   smooth.spline
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
fpi0est <- function(
  p,
  pi0_model_obj = NULL,
  z = NULL,
  pi0_model = NULL,
  indep_snps = NULL,
  weights = NULL,
  lambda = seq(0.05, 0.95, 0.05),
  constrained.p = TRUE,
  tol = 1e-9,
  maxit = 200,
  ncores = 1L,
  verbose = TRUE,
  ...
) {
  if (
    !is.null(pi0_model_obj) &&
      (is.data.frame(pi0_model_obj) || is.matrix(pi0_model_obj))
  ) {
    if (is.null(z)) {
      z <- pi0_model_obj
      pi0_model_obj <- NULL
    } else {
      stop(
        "Ambiguous input: You provided a matrix for 'pi0_model_obj' AND a 'z' argument."
      )
    }
  }

  if (!is.null(pi0_model_obj)) {
    if (
      !is.list(pi0_model_obj) || !all(c("fmod", "zt") %in% names(pi0_model_obj))
    ) {
      stop(
        "'pi0_model_obj' must be the output list from pi0_model() containing 'fmod' and 'zt'."
      )
    }

    if (is.null(z)) {
      z <- pi0_model_obj$zt
    } else {
      warning("Both 'z' and 'pi0_model_obj' were provided. Using manual 'z'.")
    }

    # Extract Formula (pi0_model)
    if (is.null(pi0_model)) {
      pi0_model <- pi0_model_obj$fmod
    }
  }

  if (is.null(z)) {
    stop(
      "Input 'z' is missing. Please provide either 'z' directly or 'pi0_model_obj'."
    )
  }

  if (!is.matrix(z) && !is.data.frame(z)) {
    stop("'z' must be a matrix or data frame.")
  }

  if (is.null(pi0_model)) {
    stop(
      "Input 'pi0_model' is missing. Please provide either 'pi0_model' directly or 'pi0_model_obj'."
    )
  }

  if (min(p, na.rm = TRUE) < 0 || max(p, na.rm = TRUE) > 1) {
    stop("P-values must be in [0, 1].")
  }

  if (verbose) {
    message("==================================================")
    message("Estimating Functional Pi0")
    message("==================================================")
  }

  # Handle missing values
  na_p <- is.na(p)
  na_z <- apply(z, 1, anyNA)
  valid_idx <- !na_z & !na_p

  p_valid <- p[valid_idx]
  z_valid <- z[valid_idx, , drop = FALSE]

  n_total <- length(p)
  n_valid <- sum(valid_idx)

  if (verbose && n_total != n_valid) {
    message(sprintf(
      "  Removed %d observations with missing values (%d remaining)",
      n_total - n_valid,
      n_valid
    ))
  }

  fit_idx <- rep(TRUE, n_valid)

  if (!is.null(indep_snps)) {
    fit_idx <- fit_idx & indep_snps[valid_idx]
  }

  if (!is.null(weights)) {
    fit_idx <- fit_idx & !is.na(weights[valid_idx]) & weights[valid_idx] > 0
  }

  p_fit <- p_valid[fit_idx]
  z_fit <- z_valid[fit_idx, , drop = FALSE]
  w_fit <- if (!is.null(weights)) weights[valid_idx][fit_idx] else NULL

  if (!is.null(w_fit)) {
    w_fit <- pmin(pmax(w_fit, 1e-6), 1.0)
  }

  # Provide safe explicit weights for fastglm
  w_fit_glm <- if (is.null(w_fit)) rep(1.0, length(p_fit)) else w_fit

  # Build formula for binomial regression
  fm <- formula(paste("phi", paste(pi0_model, collapse = " ")))

  # --- Precompute design matrix (identical across all lambda values) ---
  z_fit_df <- as.data.frame(z_fit)
  z_fit_df$phi <- 0L # placeholder response
  mf <- model.frame(fm, data = z_fit_df, na.action = na.pass)
  X <- model.matrix(attr(mf, "terms"), data = mf)

  use_parallel <- ncores > 1L && .Platform$OS.type != "windows"

  if (verbose) {
    if (use_parallel) {
      message(sprintf(
        "  Fitting models across %d lambda values (ncores = %d)...",
        length(lambda),
        ncores
      ))
    } else {
      message(sprintf(
        "  Fitting models across %d lambda values...",
        length(lambda)
      ))
    }
  }

  # Fit pi0 model at a single lambda value using precomputed X
  fit_pi0_at_lambda <- function(lambda_val) {
    y <- as.numeric(p_fit >= lambda_val)
    fam <- if (constrained.p) {
      constrained.binomial(1 - lambda_val)
    } else {
      binomial()
    }

    rval <- suppressWarnings(
      fastglm::fastglm(
        x = X,
        y = y,
        family = fam,
        weights = w_fit_glm,
        method = 3L,
        tol = tol,
        maxit = maxit
      )
    )

    # Restore the attributes required by the custom predict.fastglm2 method
    rval$terms <- attr(mf, "terms")
    rval$xlevels <- .getXlevels(attr(mf, "terms"), mf)
    rval$call <- match.call()
    class(rval) <- "fastglm2"

    return(rval)
  }

  # Fit models across all lambda values (parallel if ncores > 1)
  if (use_parallel) {
    fits_raw <- parallel::mclapply(
      lambda,
      fit_pi0_at_lambda,
      mc.cores = ncores
    )
    # Check for errors from mclapply
    failed <- vapply(fits_raw, inherits, logical(1), "try-error")
    if (any(failed)) {
      warning(sprintf(
        "%d lambda fits failed in parallel; falling back to sequential.",
        sum(failed)
      ))
      for (i in which(failed)) {
        fits_raw[[i]] <- tryCatch(
          fit_pi0_at_lambda(lambda[i]),
          error = function(e) NULL
        )
      }
    }
    fpi0_models <- data.frame(lambda = lambda)
    fpi0_models$fpi0 <- fits_raw
  } else {
    fpi0_models <- data.frame(lambda = lambda)
    fpi0_models$fpi0 <- .lapply_pb(lambda, fit_pi0_at_lambda, verbose = verbose)
  }

  # Extract fitted pi0 values for all lambda
  fpi0_matrix <- sapply(seq_along(lambda), function(i) {
    pmin(fitted.values(fpi0_models$fpi0[[i]]) / (1 - lambda[i]), 1)
  })

  ref_pi0 <- pmin(fitted.values(fpi0_models$fpi0[[1]]) / (1 - lambda[1]), 1)

  # Select optimal lambda via closed-form MISE
  if (verbose) {
    message("  Selecting optimal lambda via MISE...")
  }

  # Compute global pi0
  if (is.null(w_fit)) {
    pi0_global <- qvalue::pi0est(p_fit, pi0.method = "bootstrap")$pi0
  } else {
    pi0_global <- pi0est_weighted(
      p_fit,
      weights = w_fit,
      pi0.method = "bootstrap"
    )$pi0
  }

  # Precompute shared terms for closed-form k
  one_minus_ref <- 1 - ref_pi0
  denom <- mean(one_minus_ref^2)

  mise_stats <- vapply(
    seq_along(lambda),
    function(i) {
      fpi0_i <- fpi0_matrix[, i]
      # Closed-form solution for k minimising mean((ref - k*(1-ref) - f)^2)
      k <- mean((ref_pi0 - fpi0_i) * one_minus_ref) / denom
      # Ensure k stays strictly within the [-1, 1] bounds
      k <- max(-1, min(1, k))
      phi_hat <- pmin(ref_pi0 - k * one_minus_ref, 1)
      omega <- mean((fpi0_i - phi_hat)^2)
      delta_sq <- (max(mean(fpi0_i) - pi0_global, 0))^2
      omega + delta_sq
    },
    numeric(1)
  )

  # Guard against all-NA MISE
  min_idx <- which.min(mise_stats)
  if (length(min_idx) == 0) {
    warning("All MISE values are NA; defaulting to largest lambda.")
    min_idx <- length(lambda)
  }
  lambda_hat <- lambda[min_idx]

  if (verbose) {
    message(sprintf("  Selected lambda: %.2f", lambda_hat))
  }

  fpi0_models$chosen <- fpi0_models$lambda == lambda_hat
  fpi0_selected <- fpi0_models[fpi0_models$chosen, ]

  # Predict for all valid observations (chunk-wise for large data memory safety)
  if (verbose) {
    message("  Generating predictions...")
  }

  chunk_size <- 500000
  n_valid <- nrow(z_valid)

  if (n_valid > chunk_size) {
    fpi0_pred <- numeric(n_valid)
    n_chunks <- ceiling(n_valid / chunk_size)

    for (i in seq_len(n_chunks)) {
      idx_start <- chunk_size * (i - 1) + 1
      idx_end <- min(chunk_size * i, n_valid)

      # Extract sub-matrix and match column names implicitly through data.frame
      z_chunk <- as.data.frame(z_valid[idx_start:idx_end, , drop = FALSE])

      fpi0_pred[idx_start:idx_end] <- pmin(
        predict(
          fpi0_selected$fpi0[[1]],
          newdata = z_chunk,
          type = "response"
        ) /
          (1 - lambda_hat),
        1
      )
    }
  } else {
    fpi0_pred <- pmin(
      predict(
        fpi0_selected$fpi0[[1]],
        newdata = as.data.frame(z_valid),
        type = "response"
      ) /
        (1 - lambda_hat),
      1
    )
  }

  # Apply floor to fpi0 predictions
  fpi0_pred <- pmax(fpi0_pred, 0.01)

  # Fill results including NAs
  fpi0_out <- rep(1, length(p))
  fpi0_out[valid_idx] <- fpi0_pred

  if (verbose) {
    message("  Done.")
    message("==================================================")
    message("")
  }
  # Return S3 object
  structure(
    list(
      fpi0 = fpi0_out,
      tableLambda = fpi0_models,
      MISE = mise_stats,
      lambda = lambda_hat
    ),
    class = "fpi0"
  )
}

constrained.binomial <- function(maximum) {
  link <- structure(
    list(
      name = paste0("constrained.logit (0, ", maximum, ")"),
      linkfun = function(mu) {
        epsilon <- .Machine$double.eps
        mu <- pmax(epsilon, pmin(maximum - epsilon, mu))
        prop <- mu / maximum
        log(prop / (1 - prop))
      },
      linkinv = function(eta) {
        maximum * plogis(eta)
      },
      mu.eta = function(eta) {
        val <- maximum * dlogis(eta)
        pmax(val, .Machine$double.eps)
      },
      valideta = function(eta) TRUE
    ),
    class = "link-glm"
  )

  fam <- binomial(link)

  # Variance for constrained scale
  fam$variance <- function(mu) {
    var <- mu * (1 - mu / maximum)
    pmax(var, 1e-10)
  }

  fam$validmu <- function(mu) all(mu >= 0) && all(mu <= maximum)
  fam$family <- paste0("constrained.binomial (0, ", maximum, ")")

  fam$d2link <- function(mu) {
    p <- mu / maximum
    (1 / (1 - p)^2 - 1 / p^2) / maximum^2
  }
  fam$d3link <- function(mu) {
    p <- mu / maximum
    (2 / (1 - p)^3 + 2 / p^3) / maximum^3
  }
  fam$d4link <- function(mu) {
    p <- mu / maximum
    (6 / (1 - p)^4 - 6 / p^4) / maximum^4
  }

  fam$dvar <- function(mu) 1 - 2 * mu / maximum
  fam$d2var <- function(mu) rep(-2 / maximum, length(mu))
  fam$d3var <- function(mu) rep(0, length(mu))

  epsilon <- 1e-4
  new_init <- bquote({
    mustart <- mustart * .(maximum)
    mustart <- pmax(mustart, .(epsilon))
    mustart <- pmin(mustart, .(maximum) - .(epsilon))
  })
  fam$initialize <- as.expression(c(fam$initialize, new_init))

  fam
}

updated_fastglm <- function(
  formula,
  data,
  method = 3,
  family = gaussian(),
  weights = NULL,
  offset = NULL,
  tol = 1e-08,
  maxit = 100L,
  ...
) {
  call <- match.call()
  M <- match.call(expand.dots = FALSE)
  m <- match(
    c("formula", "data", "subset", "weights", "na.action", "offset"),
    names(M),
    0L
  )
  M <- M[c(1L, m)]
  M$drop.unused.levels <- TRUE
  M[[1L]] <- quote(stats::model.frame)
  M <- eval(M, parent.frame())
  y <- M[[1]]
  tf <- attr(M, "terms")
  X <- model.matrix(tf, M)
  weights <- as.vector(model.weights(M))
  offset <- model.offset(M)
  if (is.null(offset)) {
    offset <- rep(0, length(y))
  }
  rval <- fastglm::fastglm(
    y = y,
    x = X,
    family = family,
    method = method,
    weights = weights,
    offset = offset,
    tol = tol,
    maxit = maxit,
    ...
  )
  if (ncol(M) > 1) {
    for (i in 2:ncol(M)) {
      if (is.factor(M[, i])) {
        rval$levels[[names(M)[i]]] <- levels(M[, i])
      }
    }
  }

  rval$terms <- tf
  rval$call <- call
  rval$xlevels <- .getXlevels(tf, M)
  rval$formula <- formula
  class(rval) <- "fastglm2"
  rval
}

#' @keywords internal
#' @exportS3Method family fastglm2
# Internal S3 method for fastglm2 family extraction
family.fastglm2 <- function(object, ...) {
  object$family
}

#' @keywords internal
#' @exportS3Method predict fastglm2
# Internal S3 method for fastglm2 predictions
predict.fastglm2 <- function(
  object,
  newdata,
  type = c("link", "response"),
  na.action = na.pass,
  ...
) {
  type <- match.arg(type)
  na.act <- object$na.action
  object$na.action <- NULL

  if (missing(newdata)) {
    pred <- if (type == "link") {
      object$linear.predictors
    } else {
      fitted(object)
    }
    if (!is.null(na.act)) pred <- napredict(na.act, pred)
  } else {
    pred <- get_predict(
      object,
      newdata,
      na.action = na.action
    )
    if (type == "response") {
      pred <- family(object)$linkinv(pred)
    }
  }
  pred
}

get_predict <- function(object, newdata, na.action = na.pass, ...) {
  tt <- terms(object)
  if (missing(newdata) || is.null(newdata)) {
    if (is.null(object$fitted.values)) {
      return(object$fitted.values)
    }
  } else {
    Terms <- delete.response(tt)
    m <- model.frame(
      Terms,
      newdata,
      na.action = na.action,
      xlev = object$xlevels
    )
    if (!is.null(cl <- attr(Terms, "dataClasses"))) {
      .checkMFClasses(cl, m)
    }
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    offset <- rep(0, nrow(X))
    if (!is.null(off.num <- attr(tt, "offset"))) {
      for (i in off.num) {
        offset <- offset + eval(attr(tt, "variables")[[i + 1]], newdata)
      }
    }
    if (!is.null(object$call$offset)) {
      offset <- offset + eval(object$call$offset, newdata)
    }
  }
  p <- object$rank
  ord <- colnames(X)
  if (p < ncol(X) && !(missing(newdata) || is.null(newdata))) {
    warning("prediction from a rank-deficient fit may be misleading")
  }
  beta <- object$coefficients
  beta[is.na(beta)] <- 0
  predictor <- drop(X[, ord, drop = FALSE] %*% beta[ord])
  if (!is.null(offset)) {
    predictor <- predictor + offset
  }
  if (missing(newdata) && !is.null(na.act <- object$na.action)) {
    predictor <- napredict(na.act, predictor)
  }
  predictor
}

#' @noRd
.lapply_pb <- function(X, FUN, verbose = TRUE, pb_width = 47, ...) {
  if (!verbose || length(X) == 0) {
    return(lapply(X, FUN, ...))
  }

  pb <- txtProgressBar(min = 0, max = length(X), style = 3, width = pb_width)

  # Safety: If code crashes, this runs to close the bar (printing \n) so the console doesn't break
  on.exit(close(pb))

  result <- vector("list", length(X))
  for (i in seq_along(X)) {
    result[[i]] <- FUN(X[[i]], ...)
    setTxtProgressBar(pb, i)
  }

  # Success!
  on.exit() # 1. Cancel the safety hook

  # 2. DO NOT call close(pb). It forces a newline.
  # 3. Manually overwrite the line with spaces
  cat("\r", strrep(" ", pb_width + 10), "\r", sep = "")

  result
}

pi0est_weighted <- function(
  p,
  weights = NULL,
  lambda = seq(0.05, 0.95, 0.05),
  pi0.method = c("smoother", "bootstrap"),
  smooth.df = 3,
  smooth.log.pi0 = FALSE,
  ...
) {
  # 1. NA Handling
  rm_na <- !is.na(p)
  if (!is.null(weights)) {
    rm_na <- rm_na & !is.na(weights)
  }

  p <- p[rm_na]

  if (!is.null(weights)) {
    weights <- weights[rm_na]
    weights <- pmin(pmax(weights, 1e-6), 1.0)
  } else {
    weights <- rep(1.0, length(p))
  }

  pi0.method <- match.arg(pi0.method)
  m_eff <- sum(weights)
  lambda <- sort(lambda)
  ll <- length(lambda)

  if (min(p) < 0 || max(p) > 1) {
    stop("ERROR: p-values not in valid range [0, 1].")
  } else if (ll > 1 && ll < 4) {
    stop(sprintf(paste(
      "ERROR:",
      paste("length(lambda)=", ll, ".", sep = ""),
      "If length of lambda greater than 1, you need at least 4 values."
    )))
  } else if (min(lambda) < 0 || max(lambda) >= 1) {
    stop("ERROR: Lambda must be within [0, 1).")
  }

  # Vectorized weighted pi0: sort-once + cumsum + findInterval
  ord <- order(p)
  p_sorted <- p[ord]
  w_sorted <- weights[ord]
  w_cumsum <- cumsum(w_sorted)
  total_w <- w_cumsum[length(w_cumsum)]

  # For each lambda, W = sum of weights where p >= lambda
  W <- vapply(
    lambda,
    function(l) {
      i <- findInterval(l, p_sorted, left.open = TRUE)
      if (i == 0L) total_w else total_w - w_cumsum[i]
    },
    numeric(1)
  )

  pi0 <- W / (m_eff * (1 - lambda))
  pi0.lambda <- pi0

  if (ll == 1) {
    pi0 <- min(pi0, 1)
    pi0Smooth <- NULL
  } else {
    if (pi0.method == "smoother") {
      if (smooth.log.pi0) {
        pi0_log <- log(pi0)
        spi0 <- smooth.spline(lambda, pi0_log, df = smooth.df)
        pi0Smooth <- exp(predict(spi0, x = lambda)$y)
        pi0 <- min(pi0Smooth[ll], 1)
      } else {
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- predict(spi0, x = lambda)$y
        pi0 <- min(pi0Smooth[ll], 1)
      }
    } else if (pi0.method == "bootstrap") {
      minpi0 <- quantile(pi0, prob = 0.1)

      mse <- (W / (m_eff^2 * (1 - lambda)^2)) *
        (1 - W / m_eff) +
        (pi0 - minpi0)^2

      pi0 <- min(pi0[mse == min(mse)], 1)
      pi0Smooth <- NULL
    } else {
      stop("ERROR: pi0.method must be one of \"smoother\" or \"bootstrap\".")
    }
  }

  if (pi0 <= 0) {
    stop(
      "ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use a different range of lambda."
    )
  }

  return(list(
    pi0 = pi0,
    pi0.lambda = pi0.lambda,
    lambda = lambda,
    pi0.smooth = pi0Smooth
  ))
}
