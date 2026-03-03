#' Build Model for Functional Proportion of Null Tests (fpi0)
#'
#' @description
#' Generates a natural spline model for the functional proportion of null tests (fpi0) by
#' adaptively selecting knots based on an FDR threshold.
#'
#' @param z Matrix/data.frame of p-values (rows=tests, cols=traits).
#' @param indep_snps Logical vector indicating independent SNPs (training subset).
#'   If NULL, uses all tests. Default: NULL.
#' @param weights Optional numeric vector of weights (e.g., inverse LD scores) for weighted spline fitting. Default: NULL.
#' @param fdr_threshold FDR threshold for signal definition. Default: 0.25.
#' @param min_discoveries Min (LD-independent) significant hits required to include trait. Default: 100.
#' @param n_knots Target knot count. Default: 5. Automatically reduced if insufficient
#'   discoveries (via `min_snps_per_knot`) or capped by `max_knots`.
#' @param min_snps_per_knot Min significant SNPs per knot interval. Default: 100. Note that this is the minimum independent SNPs required per knot interval if `indep_snps` is provided.
#' @param verbose Print selection details. Default: TRUE.
#' @param seed Optional seed for reproducibility. Default: 2026.
#'
#' @details
#' Independent SNPs determine knot placement; all SNPs train the model to capture signal shape.
#' For smaller datasets (<100K), consider reducing `min_snps_per_knot` and `min_discoveries` to improve knot placement.
#'
#' @return List containing:
#' \item{fmod}{Model formula (using \code{splines::ns})}
#' \item{zt}{Data frame of globally rank-transformed p-values}
#'
#' @export
pi0_model <- function(
  z,
  indep_snps = NULL,
  weights = NULL,
  fdr_threshold = 0.25,
  min_discoveries = 150,
  min_snps_per_knot = 100,
  n_knots = 5,
  verbose = TRUE,
  seed = 2026
) {
  if (!is.null(seed)) {
    withr::local_seed(seed)
  }

  if (!is.matrix(z) && !is.data.frame(z)) {
    stop("'z' must be a matrix/data.frame")
  }

  if (!is.null(indep_snps) && !is.null(weights)) {
    stop(
      "Ambiguous LD correction: You have provided 'indep_snps' (pruning) and 'weights' (Inverse-LD). ",
      "Please use only one input."
    )
  }

  z <- as.matrix(z)
  z_ranked <- apply(z, 2, function(col) {
    n_valid <- sum(!is.na(col))
    if (n_valid == 0) {
      return(col)
    }
    (rank(col, ties.method = "random", na.last = "keep")) / n_valid
  })

  n_vars <- ncol(z)
  if (is.null(colnames(z))) {
    colnames(z) <- paste0("z", seq_len(n_vars))
  }
  colnames(z_ranked) <- make.names(colnames(z), unique = TRUE)
  var_names <- colnames(z_ranked)

  z_train_raw <- z
  z_train_ranked <- z_ranked

  use_indep <- !is.null(indep_snps)
  if (use_indep) {
    if (is.logical(indep_snps) && length(indep_snps) != nrow(z)) {
      warning(sprintf(
        "Length mismatch: indep_snps (%d) vs z rows (%d). Ignoring indep_snps.",
        length(indep_snps),
        nrow(z)
      ))
      use_indep <- FALSE
    } else if (is.numeric(indep_snps)) {
      if (
        any(indep_snps < 1, na.rm = TRUE) ||
          any(indep_snps > nrow(z), na.rm = TRUE)
      ) {
        warning(
          "Numeric indep_snps contains out-of-bounds indices. Ignoring indep_snps."
        )
        use_indep <- FALSE
      }
    } else if (is.character(indep_snps) && is.null(rownames(z))) {
      warning(
        "Character indep_snps provided, but input data 'z' has no rownames. Ignoring indep_snps."
      )
      use_indep <- FALSE
    }
  }

  if (verbose) {
    message("==================================================")
    message("Building Pi0 Model (Rank: Signal -> 0.0)")
    message(sprintf("  Variables        : %d", n_vars))

    discovery_type <- if (use_indep || !is.null(weights)) {
      "(Effective Independent)"
    } else {
      "(Raw SNPs)"
    }
    message(sprintf(
      "  Min Discoveries  : %d %s",
      min_discoveries,
      discovery_type
    ))
    message("--------------------------------------------------")
  }

  n_total_snps <- nrow(z)

  if (verbose) {
    message(sprintf(
      "  [Auto-Scale]: N = %d. Bandwidth set to %d independent SNPs per knot.",
      n_total_snps,
      min_snps_per_knot
    ))
  }

  selected_terms <- character(0)

  for (i in seq_len(n_vars)) {
    var_name <- var_names[i]

    raw_vals <- z_train_raw[, i]
    rank_vals <- z_train_ranked[, i]

    valid_idx <- !is.na(raw_vals)
    if (!is.null(weights)) {
      valid_idx <- valid_idx & !is.na(weights)
    }

    raw_vals <- raw_vals[valid_idx]
    rank_vals <- rank_vals[valid_idx]

    w_vals <- if (!is.null(weights)) {
      pmin(pmax(weights[valid_idx], 1e-6), 1.0)
    } else {
      rep(1.0, length(raw_vals))
    }

    valid_indices_map <- which(valid_idx)

    qobj <- tryCatch(qvalue::qvalue(raw_vals), error = function(e) NULL)
    knot_locs <- numeric(0)

    if (!is.null(qobj)) {
      is_sig <- qobj$qvalues <= fdr_threshold

      # Compute effective discovery count
      if (use_indep) {
        is_indep_check <- logical(length(valid_indices_map))
        if (is.logical(indep_snps)) {
          is_indep_check <- indep_snps[valid_indices_map]
        } else if (is.character(indep_snps)) {
          is_indep_check <- rownames(z)[valid_indices_map] %in% indep_snps
        } else {
          is_indep_check <- valid_indices_map %in% indep_snps
        }
        n_sig_effective <- sum(is_sig & is_indep_check)
      } else if (!is.null(weights)) {
        n_sig_effective <- sum(w_vals[is_sig])
      } else {
        n_sig_effective <- sum(is_sig)
      }

      if (n_sig_effective >= min_discoveries) {
        # Safety check: warn if large LD-inflated discovery set with no correction
        if (sum(is_sig) > 10000 && !use_indep && is.null(weights)) {
          warning(sprintf(
            "High number of raw discoveries (%d) detected for %s without LD correction. Spline knots may overfit to massive LD blocks. Consider using 'indep_snps' or 'weights'.",
            sum(is_sig),
            var_name
          ))
        }

        if (verbose) {
          message(sprintf(
            "  [%s]: Found %d %sdiscoveries (FDR < %.2f). Calculating knots...",
            var_name,
            round(n_sig_effective),
            if (use_indep || !is.null(weights)) "independent " else "",
            fdr_threshold
          ))
        }

        basis_ranks <- rank_vals[is_sig]
        basis_weights <- w_vals[is_sig]

        if (use_indep) {
          sig_indep_ranks <- rank_vals[is_sig & is_indep_check]
          sig_indep_weights <- w_vals[is_sig & is_indep_check]

          if (length(sig_indep_ranks) >= 30) {
            basis_ranks <- sig_indep_ranks
            basis_weights <- sig_indep_weights
          }
        }

        max_safe_knots <- floor(
          sum(basis_weights, na.rm = TRUE) / min_snps_per_knot
        )
        n_knots_actual <- max(1, min(max_safe_knots, n_knots))

        # WEIGHTED KNOT PLACEMENT: place knots at equal independent statistical mass
        ord <- order(basis_ranks)
        sorted_ranks <- basis_ranks[ord]
        sorted_weights <- basis_weights[ord]

        cum_w <- cumsum(sorted_weights)
        if (cum_w[length(cum_w)] > 0) {
          cum_w <- cum_w / cum_w[length(cum_w)]
          probs <- seq(0, 1, length.out = n_knots_actual + 1)

          tail_knots <- sapply(probs, function(p) {
            idx <- which(cum_w >= p)[1]
            if (is.na(idx)) {
              sorted_ranks[length(sorted_ranks)]
            } else {
              sorted_ranks[idx]
            }
          })
          knot_locs <- tail_knots[2:length(tail_knots)]
        }
      }
    }

    if (length(knot_locs) > 0) {
      knot_locs <- c(knot_locs, 0.1) # null anchor
    }

    knot_locs <- unique(round(knot_locs, 7))
    knot_locs <- knot_locs[knot_locs < 0.9999]
    knot_locs <- sort(knot_locs)

    if (length(knot_locs) > 0 && (use_indep || !is.null(weights))) {
      is_indep_subset <- rep(TRUE, length(valid_indices_map))

      if (use_indep) {
        if (is.logical(indep_snps)) {
          is_indep_subset <- indep_snps[valid_indices_map]
        } else if (is.character(indep_snps)) {
          is_indep_subset <- rownames(z)[valid_indices_map] %in% indep_snps
        } else {
          is_indep_subset <- valid_indices_map %in% indep_snps
        }
      }

      indep_ranks_subset <- rank_vals[is_indep_subset]
      indep_weights_subset <- w_vals[is_indep_subset]

      min_indep_support <- min_snps_per_knot / 2

      repeat {
        if (length(knot_locs) == 0) {
          break
        }

        bins <- c(0, knot_locs, 1)
        bin_cuts <- cut(
          indep_ranks_subset,
          breaks = bins,
          include.lowest = TRUE
        )

        counts <- as.numeric(tapply(
          indep_weights_subset,
          bin_cuts,
          sum,
          na.rm = TRUE
        ))
        counts[is.na(counts)] <- 0

        if (all(counts >= min_indep_support)) {
          break
        }

        problem_bin_idx <- which.min(counts)
        knot_to_remove <- if (problem_bin_idx > length(knot_locs)) {
          length(knot_locs)
        } else {
          problem_bin_idx
        }

        if (verbose) {
          message(sprintf(
            "  [%s]: Indep Bin %d sparse (%.1f effective SNPs). Collapsing knot %.4f.",
            var_name,
            problem_bin_idx,
            counts[problem_bin_idx],
            knot_locs[knot_to_remove]
          ))
        }
        knot_locs <- knot_locs[-knot_to_remove]
      }
    }

    calc_cutoff <- 10000 / length(raw_vals)
    if (length(knot_locs) < 1) {
      tail_rank_cutoff <- max(calc_cutoff, 0.005)
      tail_rank_cutoff <- min(tail_rank_cutoff, 0.1)
      if (verbose) {
        message(sprintf(
          "  [%s]: Fallback -> 1 Knot at %.4f",
          var_name,
          tail_rank_cutoff
        ))
      }
      knot_locs <- c(tail_rank_cutoff)
    }

    knot_str <- paste(knot_locs, collapse = ", ")
    term <- sprintf(
      "splines::ns(%s, knots = c(%s), Boundary.knots = c(0, 1))",
      var_name,
      knot_str
    )

    if (verbose) {
      model_type <- "Discovered"
      if (length(knot_locs) == 1) {
        chk_cutoff <- min(max(calc_cutoff, 0.005), 0.1)
        if (abs(knot_locs[1] - chk_cutoff) < 0.0001) {
          model_type <- "Fallback"
        }
      }
      if (model_type == "Discovered") {
        message(sprintf(
          "  [%s]: Final Model -> Spline (%d knots: %s)",
          var_name,
          length(knot_locs),
          knot_str
        ))
      }
    }

    selected_terms <- c(selected_terms, term)
  }

  if (length(selected_terms) == 0) {
    message("No variables met criteria.")
    fmod <- formula("~ 1")
  } else {
    fmod <- formula(paste("~", paste(selected_terms, collapse = " + ")))
  }

  if (verbose) {
    message("--------------------------------------------------")
    message(deparse1(fmod))
    message("==================================================\n")
  }

  list(fmod = fmod, zt = as.data.frame(z_ranked))
}
