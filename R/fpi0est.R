#' Estimate Functional Proportion of Null Tests
#'
#' @description
#' Estimates the functional proportion of null tests (pi0) using a GLM approach
#' with a constrained binomial family. This function fits models across multiple
#' lambda thresholds and selects the optimal estimate via MISE minimization.
#'
#' @details
#' **Algorithm:**
#' 1. For each lambda threshold, fit a binomial GLM: P(p >= lambda | z)
#' 2. Use constrained binomial family to ensure predictions in (0, 1)
#' 3. Select optimal lambda via MISE
#'
#' **Model Fitting:**
#' Uses \code{\link[fastglm]{fastglm}} for efficient GLM fitting with a
#' constrained binomial family that bounds predictions away from 0 and 1.
#'
#' @param p Numeric vector of p-values.
#' @param z Data frame of covariates (output from \code{\link{pi0_model}$zt}).
#' @param lambda Numeric vector of lambda thresholds for estimating pi0.
#'   Default is \code{seq(0.05, 0.95, 0.05)}.
#' @param pi0_model Formula for the pi0 model (output from \code{\link{pi0_model}$fmod}).
#' @param indep_snps Optional logical vector indicating independent SNPs for
#'   model fitting. If NULL (default), all SNPs are used.
#' @param constrained.p Logical; use constrained binomial family ensuring
#'   predictions stay within (0, 1). Default is TRUE. Recommended for stability.
#' @param tol Numeric; convergence tolerance for \code{\link[fastglm]{fastglm}}.
#'   Default is 1e-9.
#' @param maxit Integer; maximum iterations for \code{\link[fastglm]{fastglm}}.
#'   Default is 200.
#' @param verbose Logical; print progress messages. Default is TRUE.
#' @param ... Additional arguments passed to \code{\link[fastglm]{fastglm}}.
#'
#' @return A list containing:
#' \describe{
#'   \item{fpi0}{Numeric vector of functional pi0 estimates for each test.}
#'   \item{tableLambda}{A data frame summarizing results for each lambda value.}
#'   \item{MISE}{The Mean Integrated Squared Error (MISE) for the chosen model.}
#'   \item{lambda}{The selected optimal lambda value.}
#' }
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
#' fpi0_out <- fpi0est(p, z = fmod$zt, pi0_model = fmod$fmod)
#' fpi0 <- fpi0_out$fpi0
#' # Apply sffdr
#' sffdr_out <- sffdr(p, fpi0)
#'
#' @importFrom stats family predict model.matrix qnorm dnorm runif
#'   quantile na.pass optimize binomial dbinom formula fitted fitted.values
#'   gaussian median model.frame model.offset model.weights napredict
#'   delete.response terms .checkMFClasses .getXlevels dlogis plogis
#' @export
fpi0est <- function(
  p,
  z,
  lambda = seq(0.05, 0.95, 0.05),
  pi0_model = NULL,
  indep_snps = NULL,
  constrained.p = TRUE,
  tol = 1e-9,
  maxit = 200,
  verbose = TRUE,
  ...
) {
  # Input checks
  if (min(p, na.rm = TRUE) < 0 || max(p, na.rm = TRUE) > 1) {
    stop("P-values must be in [0, 1].")
  }

  if (is.null(pi0_model)) {
    stop("'pi0_model' must be provided.")
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

  # Select fitting subset (independent SNPs if provided)
  if (!is.null(indep_snps)) {
    indep_valid <- indep_snps[valid_idx]
    p_fit <- p_valid[indep_valid]
    z_fit <- z_valid[indep_valid, , drop = FALSE]
  } else {
    p_fit <- p_valid
    z_fit <- z_valid
  }

  # Build formula for binomial regression
  fm <- formula(paste("phi", paste(pi0_model, collapse = " ")))

  if (verbose) {
    message(sprintf(
      "  Fitting models across %d lambda values...",
      length(lambda)
    ))
  }

  # Fit pi0 model at a single lambda value using fastglm
  fit_pi0_at_lambda <- function(lambda_val) {
    z_fit$phi <- as.numeric(p_fit >= lambda_val)

    fam <- if (constrained.p) {
      constrained.binomial(1 - lambda_val)
    } else {
      binomial()
    }

    suppressWarnings(
      updated_fastglm(
        formula = fm,
        data = z_fit,
        family = fam,
        tol = tol,
        maxit = maxit,
        ...
      )
    )
  }

  # Fit models across all lambda values
  fpi0_models <- data.frame(lambda = lambda)
  fpi0_models$fpi0 <- lapply(lambda, fit_pi0_at_lambda)

  # # Extract fitted pi0 values for all lambda
  fpi0_matrix <- sapply(seq_along(lambda), function(i) {
    pmin(fitted.values(fpi0_models$fpi0[[i]]) / (1 - lambda[i]), 1)
  })

  ref_pi0 <- pmin(fitted.values(fpi0_models$fpi0[[1]]) / (1 - lambda[1]), 1)

  # Select optimal lambda
  if (verbose) {
    message("  Selecting optimal lambda via MISE...")
  }
  # Compute global pi0
  pi0_global <- qvalue::pi0est(p_fit, pi0.method = "bootstrap")$pi0
  mise_stats <- vapply(
    seq_along(lambda),
    function(i) {
      fpi0_i <- fpi0_matrix[, i]
      k <- optimize(
        function(k) mean((ref_pi0 - k * (1 - ref_pi0) - fpi0_i)^2),
        interval = c(-1, 1)
      )$minimum
      phi_hat <- pmin(ref_pi0 - k * (1 - ref_pi0), 1)
      omega <- mean((fpi0_i - phi_hat)^2)
      delta_sq <- (max(mean(fpi0_i) - pi0_global, 0))^2
      omega + delta_sq
    },
    numeric(1)
  )

  lambda_hat <- lambda[which.min(mise_stats)]
  if (verbose) {
    message(sprintf("  Selected lambda: %.2f", lambda_hat))
  }
  fpi0_models$chosen <- fpi0_models$lambda == lambda_hat
  fpi0_selected <- fpi0_models[fpi0_models$chosen, ]

  # Predict for all valid observations (chunk-wise for large data)
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

      fpi0_pred[idx_start:idx_end] <- pmin(
        predict(
          fpi0_selected$fpi0[[1]],
          newdata = z_valid[idx_start:idx_end, , drop = FALSE],
          type = "response"
        ) /
          (1 - lambda_hat),
        1
      )
    }
  } else {
    fpi0_pred <- pmin(
      predict(fpi0_selected$fpi0[[1]], newdata = z_valid, type = "response") /
        (1 - lambda_hat),
      1
    )
  }

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
    pred <- switch(
      type,
      link = object$linear.predictors,
      response = fitted(object)
    )
    if (!is.null(na.act)) pred <- napredict(na.act, pred)
  } else {
    pred <- get_predict(
      object,
      newdata,
      type = "response",
      na.action = na.action
    )
    switch(
      type,
      response = {
        pred <- family(object)$linkinv(pred)
      },
      link =
    )
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
