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
#' @param indep_snps Optional independent SNPs for model fitting: a logical
#'   vector with one entry per p-value, row indices, or SNP ids matched against
#'   \code{rownames(z)}. Default is NULL (all SNPs used).
#' @param weights Optional numeric vector of weights in (0, 1], one per p-value,
#'   for the model fit. For GWAS data, this should be inverse LD scores (values
#'   above 1 are clamped to 1 with a warning). If these are available
#'   then you don't have to use the indep_snps argument, as the weights will effectively prioritize independent SNPs in the fitting process.
#' @param lambda Numeric vector of lambda thresholds in (0, 1). Default is \code{seq(0.05, 0.95, 0.05)}.
#' @param constrained.p Logical; use constrained binomial family. Default is TRUE.
#' @param tol Numeric; convergence tolerance. Default is 1e-9.
#' @param maxit Integer; maximum iterations of each fit, which also caps
#'   fastglm's step halvings. Default is 60.
#' @param ncores Integer; number of cores for parallel lambda fitting via
#'   \code{parallel::mclapply}. Fork-based parallelism (Unix/macOS only) -
#'   on Windows this falls back to sequential execution. Default is 1L.
#' @param verbose Logical; print progress messages. Default is TRUE.
#' @param folds Optional vector of fold labels (e.g., chromosomes), one per
#'   p-value. If given, lambda is chosen once and each fold's pi0 comes from the
#'   model refitted on the other folds, so a SNP's pi0 does not use its own
#'   p-value or its LD neighbours'. The fold models are kept with the estimates,
#'   and \code{sffdr()} with the same \code{folds} and \code{nfolds} scores each
#'   fold with its own model. Recommended with \code{sffdr(folds = )} when no
#'   LD-independent set (\code{indep_snps}) is available. Default is NULL.
#' @param nfolds Optional integer; with \code{folds}, the number of folds (the
#'   labels grouped as in \code{\link{sffdr}}). Default is NULL.
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
#'   \item{fpi0}{Numeric vector of functional pi0 estimates for each test,
#'     floored at 0.01. Tests with a missing p-value or covariate get 1. With
#'     \code{folds}, each test's estimate is from the model fitted without its
#'     fold, and the fold models are the attribute \code{"pi0_folds"}. With
#'     one covariate in the model, the covariate is the attribute
#'     \code{"surrogate"}, \code{sffdr()}'s default surrogate.}
#'   \item{tableLambda}{A data frame with one row per lambda value: the
#'     MISE, the mean fitted pi0 and whether it was chosen.}
#'   \item{MISE}{The Mean Integrated Squared Error (MISE) at each lambda.}
#'   \item{lambda}{The selected optimal lambda value.}
#'   \item{model}{The fitted model at the selected lambda (class
#'     \code{fastglm2}; \code{predict()} gives \eqn{P(p \ge \lambda | z)}).}
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
  maxit = 60,
  ncores = 1L,
  verbose = TRUE,
  folds = NULL,
  nfolds = NULL,
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

  if (!is.numeric(p)) {
    stop("'p' must be a numeric vector of p-values.")
  }
  if (length(p) != nrow(z)) {
    stop("'p' must have one entry per row of 'z'.")
  }
  if (min(p, na.rm = TRUE) < 0 || max(p, na.rm = TRUE) > 1) {
    stop("P-values must be in [0, 1].")
  }
  if (!is.null(weights)) {
    if (!is.numeric(weights) || length(weights) != length(p)) {
      stop("'weights' must be a numeric vector with one entry per p-value.")
    }
    if (any(weights > 1, na.rm = TRUE)) {
      warning(
        "'weights' above 1 are clamped to 1: weights are expected in (0, 1], e.g. inverse LD scores, not LD scores."
      )
    }
  }
  if (
    !is.numeric(lambda) || length(lambda) == 0 || anyNA(lambda) ||
      any(lambda <= 0 | lambda >= 1)
  ) {
    stop("'lambda' must be numeric values in (0, 1).")
  }
  lambda <- sort(unique(lambda))
  if (!is.null(folds)) {
    if (length(folds) != length(p)) stop("'folds' must have one label per p-value.")
    if (!is.null(nfolds)) folds <- .group_folds(folds, nfolds)
    if (length(unique(folds[!is.na(folds)])) < 2L) stop("'folds' must have at least two distinct labels.")
  }
  train <- .train_set(indep_snps, length(p), rownames(z))

  if (verbose) {
    message("==================================================")
    message("Estimating Functional Pi0")
    message("==================================================")
  }

  # Handle missing values
  na_p <- is.na(p)
  na_z <- rowSums(is.na(z)) > 0
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

  fit_idx <- train[valid_idx]

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

  # MISE and mean of one lambda's pi0, against the first lambda's (ref_pi0)
  ref_pi0 <- one_minus_ref <- denom <- NULL
  pi0_of <- function(fit, l) pmin(fitted.values(fit) / (1 - l), 1)
  stat_of <- function(fpi0_i) {
    # Closed-form solution for k minimising mean((ref - k*(1-ref) - f)^2)
    k <- mean((ref_pi0 - fpi0_i) * one_minus_ref) / denom
    # Ensure k stays strictly within the [-1, 1] bounds
    k <- max(-1, min(1, k))
    phi_hat <- pmin(ref_pi0 - k * one_minus_ref, 1)
    omega <- mean((fpi0_i - phi_hat)^2)
    delta_sq <- (max(mean(fpi0_i) - pi0_global, 0))^2
    c(mise = omega + delta_sq, mean = colMeans(matrix(fpi0_i)))
  }
  set_ref <- function(f) { ref_pi0 <<- f; one_minus_ref <<- 1 - f; denom <<- mean(one_minus_ref^2) }

  # Fit models across all lambda values (parallel if ncores > 1)
  fit_hat <- NULL
  if (use_parallel) {
    fits <- parallel::mclapply(
      lambda,
      fit_pi0_at_lambda,
      mc.cores = ncores
    )
    # Check for errors from mclapply
    failed <- vapply(fits, inherits, logical(1), "try-error")
    if (any(failed)) {
      warning(sprintf(
        "%d lambda fits failed in parallel; falling back to sequential.",
        sum(failed)
      ))
      for (i in which(failed)) {
        fits[[i]] <- tryCatch(
          fit_pi0_at_lambda(lambda[i]),
          error = function(e) NULL
        )
      }
    }
    set_ref(pi0_of(fits[[1]], lambda[1]))
    st <- lapply(seq_along(lambda), function(i) stat_of(pi0_of(fits[[i]], lambda[i])))
  } else {
    # One lambda at a time: a fit is kept only while its MISE is the lowest
    best <- Inf
    st <- .lapply_pb(seq_along(lambda), function(i) {
      fit <- fit_pi0_at_lambda(lambda[i]); f <- pi0_of(fit, lambda[i])
      if (i == 1L) set_ref(f)
      s <- stat_of(f)
      if (isTRUE(s[["mise"]] < best)) { best <<- s[["mise"]]; fit_hat <<- fit }
      s
    }, verbose = verbose)
  }
  mise_stats <- vapply(st, function(s) s[["mise"]], numeric(1))
  mean_fpi0 <- vapply(st, function(s) s[["mean"]], numeric(1))

  # Select optimal lambda via closed-form MISE
  if (verbose) {
    message("  Selecting optimal lambda via MISE...")
  }

  # Guard against all-NA MISE
  min_idx <- which.min(mise_stats)
  if (length(min_idx) == 0) {
    warning("All MISE values are NA; defaulting to largest lambda.")
    min_idx <- length(lambda)
    fit_hat <- NULL
  }
  lambda_hat <- lambda[min_idx]

  if (verbose) {
    message(sprintf("  Selected lambda: %.2f", lambda_hat))
  }

  if (use_parallel) { fit_hat <- fits[[min_idx]]; fits <- NULL }
  if (is.null(fit_hat)) fit_hat <- fit_pi0_at_lambda(lambda_hat)
  attr(fit_hat$terms, ".Environment") <- topenv()      # not this frame
  if (isFALSE(fit_hat$converged)) {
    warning(sprintf(
      "The pi0 model at lambda = %.2f did not converge in %d iterations; consider increasing 'maxit' or simplifying the model.",
      lambda_hat,
      maxit
    ))
  }

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
          fit_hat,
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
        fit_hat,
        newdata = as.data.frame(z_valid),
        type = "response"
      ) /
        (1 - lambda_hat),
      1
    )
  }

  # Apply floor to fpi0 predictions
  fpi0_pred <- pmax(fpi0_pred, 0.01)

  # Folds: refit at lambda_hat without each fold, keep the models for sffdr()
  pi0_folds <- NULL
  if (!is.null(folds)) {
    fv <- folds[valid_idx]; ff <- fv[fit_idx]; ids <- sort(unique(fv[!is.na(fv)]))
    y_hat <- as.numeric(p_fit >= lambda_hat)
    fam <- if (constrained.p) constrained.binomial(1 - lambda_hat) else binomial()
    if (verbose) message(sprintf("  Refitting in %d folds at lambda = %.2f...", length(ids), lambda_hat))
    refit <- function(f) {
      tr <- !is.na(ff) & ff != f
      suppressWarnings(fastglm::fastglm(x = X[tr, , drop = FALSE], y = y_hat[tr], family = fam, weights = w_fit_glm[tr],
                                        method = 3L, tol = tol, maxit = maxit))$coefficients
    }
    tt <- attr(mf, "terms"); attr(tt, ".Environment") <- topenv()          # not this frame (every lambda's fit)
    pi0_folds <- list(labels = folds, ids = ids, coef = do.call(cbind, .lapply_pb(ids, refit, verbose = verbose)), coef_all = fit_hat$coefficients, lambda = lambda_hat,
                      terms = tt, xlevels = .getXlevels(tt, mf), constrained = constrained.p, z = z)
    for (f in ids) {
      te <- which(!is.na(fv) & fv == f)
      fpi0_pred[te] <- .fpi0_of_fold(pi0_folds, f, rows = which(valid_idx)[te])
    }
  }

  # Fill results: tests with a missing p-value or covariate get fpi0 = 1
  fpi0_out <- rep(1, length(p))
  fpi0_out[valid_idx] <- fpi0_pred
  attr(fpi0_out, "pi0_folds") <- pi0_folds

  # One covariate: the default surrogate, missing values at the null end
  v <- intersect(all.vars(delete.response(attr(mf, "terms"))), colnames(z))
  if (length(v) == 1L && is.numeric(s <- if (is.data.frame(z)) z[[v]] else z[, v])) {
    if (anyNA(s)) {
      r <- suppressWarnings(cor(s, fpi0_out, method = "spearman", use = "complete.obs"))
      s <- if (isTRUE(r < 0)) -s else s
      s[is.na(s)] <- Inf
    }
    attr(fpi0_out, "surrogate") <- as.numeric(s)
  }

  if (verbose) {
    message("  Done.")
    message("==================================================")
    message("")
  }
  # Return S3 object; only the chosen model is kept (each fit holds several
  # vectors of length n)
  structure(
    list(
      fpi0 = fpi0_out,
      tableLambda = data.frame(
        lambda = lambda,
        MISE = mise_stats,
        mean_fpi0 = mean_fpi0,
        chosen = lambda == lambda_hat
      ),
      MISE = mise_stats,
      lambda = lambda_hat,
      model = fit_hat
    ),
    class = "fpi0"
  )
}

## pi0 under fold k's model (k NULL: the global model), 1 for a missing covariate
.fpi0_of_fold <- function(pf, k, rows = seq_len(nrow(pf$z)), chunk_size = 500000) {
  tt <- delete.response(pf$terms); b <- if (is.null(k)) pf$coef_all else pf$coef[, match(k, pf$ids)]; b[is.na(b)] <- 0
  out <- rep(1, length(rows)); zz <- pf$z[rows, , drop = FALSE]; ok <- which(stats::complete.cases(zz))
  for (i in split(ok, ceiling(seq_along(ok) / chunk_size))) {
    X <- model.matrix(tt, model.frame(tt, as.data.frame(zz[i, , drop = FALSE]), xlev = pf$xlevels))
    mu <- plogis(drop(X[, names(b), drop = FALSE] %*% b))
    out[i] <- pmax(pmin(if (pf$constrained) mu else mu / (1 - pf$lambda), 1), 0.01)
  }
  out
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

  # Progress bar with the time left
  n <- length(X)
  t0 <- proc.time()[["elapsed"]]
  show <- function(i) {
    k <- round(pb_width * i / n)
    left <- if (i > 0 && i < n) sprintf(", about %s left  ", .time_left((proc.time()[["elapsed"]] - t0) / i * (n - i))) else strrep(" ", 24)
    cat(sprintf("\r  |%s%s| %3d%%%s", strrep("=", k), strrep(" ", pb_width - k), round(100 * i / n), left))
  }

  # Safety: If code crashes, this ends the line so the console doesn't break
  on.exit(cat("\n"))
  show(0)

  result <- vector("list", n)
  for (i in seq_len(n)) {
    result[[i]] <- FUN(X[[i]], ...)
    show(i)
  }

  # Success: cancel the safety hook and overwrite the line with spaces
  on.exit()
  cat("\r", strrep(" ", pb_width + 40), "\r", sep = "")

  result
}

.time_left <- function(s) {
  if (s < 90) sprintf("%.0f s", s) else if (s < 5400) sprintf("%.0f min", s / 60) else sprintf("%.1f h", s / 3600)
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
