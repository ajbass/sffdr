#' Build the fpi0 Model for an Expected Chi-Squared Covariate
#'
#' @description
#' Builds the natural-spline model for the functional proportion of null tests
#' (fpi0) when the informative variable is an expected chi-squared statistic
#' per SNP, such as a prediction of \eqn{E[\chi^2]} from functional annotations
#' (the LD-score regression fit of Bass et al., 2026), rather than a p-value.
#' It is the companion of \code{\link{pi0_model}} for this kind of covariate:
#' the output goes to \code{\link{fpi0est}} in the same way, and the covariate
#' itself serves as the \code{surrogate} of \code{\link{sffdr}}.
#'
#' @param x Numeric vector, matrix or data frame of expected chi-squared
#'   statistics (rows = tests, one column per prediction). Further columns can
#'   be annotations (e.g., the LD score, MAF, annotation LD scores): each
#'   column gets its own spline. Any scale is accepted; only the ranks are
#'   modelled, and \code{null_value} is read on the scale of \code{x}.
#' @param indep_snps Independent SNPs (the training subset): a logical vector
#'   with one entry per row of \code{x}, row indices, or SNP ids matched
#'   against \code{rownames(x)}. If NULL, every test counts. Default: NULL.
#' @param weights Optional numeric vector of weights in (0, 1] (e.g., inverse
#'   LD scores); a bin's support is then its total weight. Default: NULL.
#' @param null_value The value of the covariate under no association, where
#'   the anchor knot is placed. Default: 1, the expectation of a chi-squared
#'   statistic with one degree of freedom under the null.
#' @param knots Candidate knots, as quantiles of the rank of \code{x}.
#'   Default: \code{c(0.5, 0.75, 0.9, 0.95, 0.99, 0.999)}.
#' @param min_snps_per_knot Minimum number of independent SNPs (or total
#'   weight) between consecutive knots. Default: 100.
#' @param verbose Print the knots kept for each column. Default: TRUE.
#' @param seed Seed for the random tie-breaking of the ranks. Default: 2026.
#'
#' @details
#' A p-value covariate puts its signal at rank 0, and \code{\link{pi0_model}}
#' places its knots at the quantiles of that column's discoveries. An expected
#' chi-squared puts its signal at the top and has a known null value: with no
#' association \eqn{E[\chi^2] = 1}, so SNPs predicted at or below
#' \code{null_value} carry no information and the proportion of nulls is flat
#' there. The model is a natural spline in the rank of the covariate (so its
#' scale is immaterial) with one knot anchored at the rank of
#' \code{null_value}, i.e. the share of SNPs predicted at or below it, and the
#' remaining knots at fixed upper quantiles of the rank, where the proportion
#' of nulls changes. Candidate knots at or below the anchor are dropped.
#'
#' Knots are then removed, sparsest bin first, until every bin between
#' consecutive knots holds at least \code{min_snps_per_knot} independent SNPs
#' (or that much weight), so no spline coefficient rests on a handful of SNPs;
#' the anchor is kept while a neighbouring knot can go instead. The anchor is
#' used only when \code{null_value} sits inside the covariate's range (rank
#' between 0.01 and 0.99): a covariate on another scale, or one predicted
#' above 1 everywhere, gets the quantile knots alone, with a message, which is
#' how an annotation column is handled. With a single column the covariate is
#' also the natural \code{surrogate} for \code{\link{sffdr}}; with several,
#' use the estimated fpi0 (the default of \code{\link{sffdr}}).
#'
#' @return List containing:
#' \item{fmod}{Model formula (using \code{splines::ns})}
#' \item{zt}{Data frame of rank-transformed covariates}
#' \item{knots}{List of the knots kept, per column, on the rank scale}
#' \item{anchor}{Named numeric vector: the rank of \code{null_value} per
#'   column (NA when the anchor was not used)}
#'
#' @examples
#' set.seed(1)
#' n <- 20000
#' a <- rnorm(n)
#' echisq <- 1 + exp(a) / 4
#' nonnull <- runif(n) < plogis(-4 + a)
#' z <- rnorm(n) + ifelse(nonnull, rnorm(n, 0, 3), 0)
#' p <- 2 * pnorm(-abs(z))
#' mod <- pi0_model_echisq(echisq, verbose = FALSE)
#' fpi0 <- fpi0est(p, pi0_model_obj = mod, verbose = FALSE)$fpi0
#' cor(fpi0, echisq) < 0
#'
#' @references
#' Bass AJ, Park J, Mukherjee S, Whittaker JC, Wallace C (2026). Genomic
#' foundation models do not yet boost common variant discovery beyond standard
#' functional annotations. bioRxiv.
#'
#' @seealso \code{\link{pi0_model}}, \code{\link{fpi0est}}, \code{\link{sffdr}}
#' @export
pi0_model_echisq <- function(
  x,
  indep_snps = NULL,
  weights = NULL,
  null_value = 1,
  knots = c(0.5, 0.75, 0.9, 0.95, 0.99, 0.999),
  min_snps_per_knot = 100,
  verbose = TRUE,
  seed = 2026
) {
  if (!is.null(seed)) {
    withr::local_seed(seed)
  }
  if (is.null(dim(x))) {
    x <- matrix(x, ncol = 1, dimnames = list(names(x), "echisq"))
  }
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a numeric vector, matrix or data frame.")
  }
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("'x' must be numeric: expected chi-squared statistics, one row per test.")
  }
  if (!is.null(indep_snps) && !is.null(weights)) {
    stop(
      "Ambiguous LD correction: You have provided 'indep_snps' (pruning) and 'weights' (Inverse-LD). ",
      "Please use only one input."
    )
  }
  if (!is.null(weights) && length(weights) != nrow(x)) {
    stop("'weights' must have one entry per row of 'x'.")
  }
  if (!is.numeric(null_value) || length(null_value) != 1 || is.na(null_value)) {
    stop("'null_value' must be a single number.")
  }
  if (
    !is.numeric(knots) || length(knots) == 0 || anyNA(knots) ||
      any(knots <= 0 | knots >= 1)
  ) {
    stop("'knots' must be quantiles in (0, 1).")
  }
  knots <- sort(unique(knots))
  if (!is.numeric(min_snps_per_knot) || min_snps_per_knot < 1) {
    stop("'min_snps_per_knot' must be at least 1.")
  }

  if (is.null(colnames(x))) {
    colnames(x) <- if (ncol(x) == 1) "echisq" else paste0("echisq", seq_len(ncol(x)))
  }
  colnames(x) <- make.names(colnames(x), unique = TRUE)
  train <- .train_set(indep_snps, nrow(x), rownames(x))

  # Rank transform: the top of the covariate is the signal end (rank -> 1)
  x_ranked <- apply(x, 2, function(col) {
    n_valid <- sum(!is.na(col))
    if (n_valid == 0) {
      return(col)
    }
    rank(col, ties.method = "random", na.last = "keep") / n_valid
  })
  x_ranked <- matrix(x_ranked, ncol = ncol(x), dimnames = dimnames(x))

  if (verbose) {
    message("==================================================")
    message("Building Pi0 Model (Expected chi-squared: Signal -> 1.0)")
    message(sprintf("  Variables        : %d", ncol(x)))
    message(sprintf("  Null value       : %g", null_value))
    message(sprintf("  Min SNPs per bin : %d", as.integer(min_snps_per_knot)))
    message("--------------------------------------------------")
  }

  terms <- character(0)
  kept <- vector("list", ncol(x))
  names(kept) <- colnames(x)
  anchors <- rep(NA_real_, ncol(x))
  names(anchors) <- colnames(x)

  for (i in seq_len(ncol(x))) {
    var_name <- colnames(x)[i]
    raw_vals <- x[, i]
    rank_vals <- x_ranked[, i]
    valid <- !is.na(raw_vals)
    if (!is.null(weights)) {
      valid <- valid & !is.na(weights)
    }
    if (sum(valid) < 2) {
      stop(sprintf("[%s]: fewer than two non-missing values.", var_name))
    }

    # The support of a bin: independent SNPs, total weight, or every SNP
    support_ranks <- rank_vals[valid & train]
    support_w <- if (!is.null(weights)) {
      pmin(pmax(weights[valid & train], 1e-6), 1.0)
    } else {
      rep(1.0, length(support_ranks))
    }

    # The anchor: the share of SNPs predicted at or below the null value
    anchor <- mean(raw_vals[valid] <= null_value)
    use_anchor <- anchor > 0.01 && anchor < 0.99
    current <- if (use_anchor) {
      sort(c(anchor, knots[knots > anchor + 0.01]))
    } else {
      knots
    }
    if (verbose) {
      if (use_anchor) {
        message(sprintf(
          "  [%s]: anchor at rank %.4f (share of SNPs with x <= %g); candidates below it dropped.",
          var_name, anchor, null_value
        ))
      } else {
        message(sprintf(
          "  [%s]: %.1f%% of SNPs at or below %g, so no anchor: quantile knots only.",
          var_name, 100 * anchor, null_value
        ))
      }
    }

    # Remove knots, sparsest bin first, until every bin has enough support;
    # a bin next to the anchor loses its other knot, so the anchor stays
    final <- current
    repeat {
      if (length(final) == 0) {
        break
      }
      bins <- cut(support_ranks, breaks = c(0, final, 1), include.lowest = TRUE)
      counts <- as.numeric(tapply(support_w, bins, sum))
      counts[is.na(counts)] <- 0
      if (all(counts >= min_snps_per_knot)) {
        break
      }
      bad <- which.min(counts)                  # bins are 1..(k + 1); knot j closes bin j
      a_idx <- if (use_anchor) match(anchor, final) else NA_integer_
      drop <- NULL
      if (!is.na(a_idx)) {
        if (bad == a_idx && a_idx > 1) {
          drop <- a_idx - 1                      # the bin below the anchor: drop its lower knot
        } else if (bad == a_idx + 1 && a_idx < length(final)) {
          drop <- a_idx + 1                      # the bin above the anchor: drop its upper knot
        }
      }
      if (is.null(drop)) {
        drop <- min(bad, length(final))
      }
      if (verbose) {
        message(sprintf(
          "  [%s]: bin %d sparse (%.1f SNPs). Removing knot %.4f.",
          var_name, bad, counts[bad], final[drop]
        ))
      }
      if (!is.na(a_idx) && drop == a_idx) {
        use_anchor <- FALSE
      }
      final <- final[-drop]
    }

    if (length(final) > 0) {
      term <- sprintf(
        "splines::ns(%s, knots = c(%s), Boundary.knots = c(0, 1))",
        var_name, paste(final, collapse = ", ")
      )
      if (verbose) {
        message(sprintf(
          "  [%s]: Final Model -> Spline (%d knots: %s)",
          var_name, length(final), paste(round(final, 5), collapse = ", ")
        ))
      }
    } else {
      term <- var_name
      if (verbose) {
        message(sprintf("  [%s]: Final Model -> Linear", var_name))
      }
    }
    terms <- c(terms, term)
    kept[[i]] <- final
    anchors[i] <- if (use_anchor) anchor else NA_real_
  }

  fmod <- formula(paste("~", paste(terms, collapse = " + ")))
  if (verbose) {
    message("--------------------------------------------------")
    message(deparse1(fmod))
    message("==================================================\n")
  }

  list(
    fmod = fmod,
    zt = as.data.frame(x_ranked),
    knots = kept,
    anchor = anchors
  )
}
