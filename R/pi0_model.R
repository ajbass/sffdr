#' Build Model for Functional Proportion of Null Tests (fpi0)
#'
#' @description
#' Generates a natural spline model for the functional proportion of null tests (fpi0) by
#' adaptively selecting knots based on an FDR threshold.
#'
#' @param z Matrix/data.frame of p-values (rows=tests, cols=traits).
#' @param indep_snps Logical vector indicating independent SNPs (training subset).
#'   If NULL, uses all tests. Default: NULL.
#' @param fdr_threshold FDR threshold for signal definition. Default: 0.25.
#' @param min_discoveries Min significant hits required to include trait. Default: 2500.
#' @param n_knots Target knot count. Default: 5. Automatically reduced if insufficient
#'   discoveries (via `min_snps_per_knot`) or capped by `max_knots`.
#' @param min_snps_per_knot Min significant SNPs per knot interval. Default: 2500.
#' @param verbose Print selection details. Default: TRUE.
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
  fdr_threshold = 0.25,
  min_discoveries = 2500,
  min_snps_per_knot = 2500,
  n_knots = 5,
  verbose = TRUE
) {
  if (!is.matrix(z) && !is.data.frame(z)) {
    stop("'z' must be a matrix/data.frame")
  }

  # --- 1. Rank Transform (Signal -> 0.0) ---
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

  # Use ALL SNPs for Training (Shape Discovery)
  z_train_raw <- z
  z_train_ranked <- z_ranked

  if (verbose) {
    message("==================================================")
    message("Building Pi0 Model (Rank: Signal -> 0.0)")
    message(sprintf("  Variables        : %d", n_vars))
    message(sprintf("  Min Discoveries  : %d (All SNPs)", min_discoveries))
    message("--------------------------------------------------")
  }
  n_total_snps <- nrow(z)

  # Default bandwidth safety
  dynamic_min_snps <- min_snps_per_knot

  if (verbose) {
    message(sprintf(
      "  [Auto-Scale]: N = %d. Bandwidth set to %d SNPs per knot.",
      n_total_snps,
      dynamic_min_snps
    ))
  }

  selected_terms <- character(0)

  for (i in seq_len(n_vars)) {
    var_name <- var_names[i]

    raw_vals <- z_train_raw[, i]
    rank_vals <- z_train_ranked[, i]

    valid_idx <- !is.na(raw_vals)
    raw_vals <- raw_vals[valid_idx]
    rank_vals <- rank_vals[valid_idx]

    valid_indices_map <- which(valid_idx)

    qobj <- tryCatch(qvalue::qvalue(raw_vals), error = function(e) NULL)
    knot_locs <- numeric(0)

    if (!is.null(qobj)) {
      # Identify discoveries in the FULL set
      is_sig <- qobj$qvalues <= fdr_threshold
      n_sig <- sum(is_sig)

      if (n_sig >= min_discoveries) {
        if (verbose) {
          message(sprintf(
            "  [%s]: Found %d discoveries (FDR < %.2f). Calculating knots...",
            var_name,
            n_sig,
            fdr_threshold
          ))
        }

        basis_ranks <- rank_vals[is_sig]

        if (!is.null(indep_snps)) {
          is_indep <- logical(length(valid_indices_map))

          if (is.logical(indep_snps)) {
            if (length(indep_snps) != nrow(z)) {
              warning(sprintf(
                "Length mismatch: indep_snps (%d) vs z rows (%d). Ignoring indep_snps.",
                length(indep_snps),
                nrow(z)
              ))
            } else {
              is_indep <- indep_snps[valid_indices_map]
            }
          } else if (is.character(indep_snps) && !is.null(rownames(z))) {
            current_names <- rownames(z)[valid_indices_map]
            is_indep <- current_names %in% indep_snps
          } else {
            is_indep <- valid_indices_map %in% indep_snps
          }

          # Intersection: Significant AND Independent
          sig_indep_ranks <- rank_vals[is_sig & is_indep]

          # Only switch to robust basis if we have enough data points (e.g. > 30)
          if (length(sig_indep_ranks) >= 30) {
            basis_ranks <- sig_indep_ranks
          }
        }

        # Calculate max knots based on TOTAL signal (n_sig)
        max_safe_knots <- floor(length(basis_ranks) / min_snps_per_knot)
        n_knots_actual <- max(1, min(max_safe_knots, n_knots))

        probs <- seq(0, 1, length.out = n_knots_actual + 1)
        tail_knots <- quantile(basis_ranks, probs = probs)
        knot_locs <- tail_knots[2:length(tail_knots)]
      }
    }

    if (length(knot_locs) > 0) {
      knot_locs <- c(knot_locs, 0.1) # null anchor
    }

    knot_locs <- unique(round(knot_locs, 7))
    knot_locs <- knot_locs[knot_locs < 0.9999]
    knot_locs <- sort(knot_locs)

    if (length(knot_locs) > 0 && !is.null(indep_snps)) {
      is_indep_subset <- logical(length(valid_indices_map))

      if (is.logical(indep_snps)) {
        if (length(indep_snps) == nrow(z)) {
          is_indep_subset <- indep_snps[valid_indices_map]
        }
      } else if (is.character(indep_snps) && !is.null(rownames(z))) {
        curr_names <- rownames(z)[valid_indices_map]
        is_indep_subset <- curr_names %in% indep_snps
      } else {
        is_indep_subset <- valid_indices_map %in% indep_snps
      }

      indep_ranks_subset <- rank_vals[is_indep_subset]
      # minimum # of independent SNPs required in each knot interval
      min_indep_support <- 50

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
        counts <- as.numeric(table(bin_cuts))

        if (all(counts >= min_indep_support)) {
          break
        }

        problem_bin_idx <- which.min(counts)

        if (problem_bin_idx > length(knot_locs)) {
          knot_to_remove <- length(knot_locs)
        } else {
          knot_to_remove <- problem_bin_idx
        }

        if (verbose) {
          message(sprintf(
            "  [%s]: Indep Bin %d sparse (%d SNPs). Collapsing knot %.4f.",
            var_name,
            problem_bin_idx,
            counts[problem_bin_idx],
            knot_locs[knot_to_remove]
          ))
        }
        knot_locs <- knot_locs[-knot_to_remove]
      }
    }

    if (length(knot_locs) < 1) {
      calc_cutoff <- 10000 / length(raw_vals)
      tail_rank_cutoff <- max(calc_cutoff, 0.005) # Never go below 0.5% rank
      tail_rank_cutoff <- min(max(calc_cutoff, 0.005), 0.1) # Cap at 10% rank for fallback
      if (verbose) {
        message(sprintf(
          "  [%s]: Fallback -> 1 Knot at %.4f",
          var_name,
          tail_rank_cutoff
        ))
      }
      # Assign the fallback knot!
      knot_locs <- c(tail_rank_cutoff)
    }

    knot_str <- paste(knot_locs, collapse = ", ")
    term <- sprintf(
      "splines::ns(%s, knots = c(%s), Boundary.knots = c(0, 1))",
      var_name,
      knot_str
    )

    if (verbose) {
      # Detect if we are in fallback mode for messaging
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
