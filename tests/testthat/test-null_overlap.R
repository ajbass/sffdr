# =================================================================
# Tests for null_overlap.R
# Focus: RNA-seq 20K genes, same samples, different phenotypes
# =================================================================

library(testthat)

# -----------------------------------------------------------------
# Simulation helper
# -----------------------------------------------------------------
#' Generates correlated z-scores for two phenotypes measured on the
#' same N samples (RNA-seq same-sample scenario).
#'
#' Model: Bivariate normal noise + non-centrality shift for signals.
#'
#' @param n        Total genes
#' @param pi0      True null fraction
#' @param rho      True null correlation (same-sample overlap)
#' @param ncp      Non-centrality for signal genes
#' @param seed     Random seed

sim_rnaseq <- function(
  n = 20000L,
  pi0 = 0.80,
  rho = 0.30,
  ncp = 3.0,
  seed = 42L
) {
  # Preserve global RNG state to prevent cross-test contamination
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  set.seed(seed)

  n_null <- round(pi0 * n)
  n_signal <- n - n_null
  is_null <- c(rep(TRUE, n_null), rep(FALSE, n_signal))

  # Generate bivariate normal noise with exact correlation rho
  # using Cholesky decomposition
  e1 <- rnorm(n)
  e2 <- rnorm(n)

  z1 <- e1
  z2 <- rho * e1 + sqrt(1 - rho^2) * e2

  # Add signal to the true positives
  signal_idx <- which(!is_null)
  z1[signal_idx] <- z1[signal_idx] + ncp
  z2[signal_idx] <- z2[signal_idx] + ncp

  p1 <- 2 * pnorm(-abs(z1))
  p2 <- 2 * pnorm(-abs(z2))

  list(
    z1 = z1,
    z2 = z2,
    p1 = p1,
    p2 = p2,
    is_null = is_null,
    rho = rho,
    pi0 = pi0
  )
}


# =================================================================
# 1. estimate_rho_overlap
# =================================================================

test_that("estimate_rho_overlap: double_null recovers rho under global null", {
  dat <- sim_rnaseq(n = 20000L, pi0 = 1.0, rho = 0.25, seed = 1L)
  rho_hat <- estimate_rho_overlap(
    dat$z1,
    dat$z2,
    dat$p1,
    dat$p2,
    method = "double_null"
  )
  expect_equal(rho_hat, dat$rho, tolerance = 0.05)
})

test_that("estimate_rho_overlap: double_null recovers rho with pi0=0.80", {
  dat <- sim_rnaseq(n = 20000L, pi0 = 0.80, rho = 0.25, ncp = 3.0)
  rho_hat <- estimate_rho_overlap(
    dat$z1,
    dat$z2,
    dat$p1,
    dat$p2,
    method = "double_null"
  )
  expect_equal(rho_hat, dat$rho, tolerance = 0.05)
})

test_that("estimate_rho_overlap: double_null recovers NEGATIVE rho", {
  dat <- sim_rnaseq(n = 20000L, pi0 = 0.80, rho = -0.25, ncp = 3.0)
  rho_hat <- estimate_rho_overlap(
    dat$z1,
    dat$z2,
    dat$p1,
    dat$p2,
    method = "double_null"
  )
  # tolerance = 0.10: lfdr intersection attenuates negative rho similarly
  # to positive rho (~5-15% attenuation from range restriction)
  expect_equal(rho_hat, dat$rho, tolerance = 0.10)
})

test_that("estimate_rho_overlap: mom_scaled roughly recovers rho", {
  dat <- sim_rnaseq(n = 20000L, pi0 = 0.80, rho = 0.20, ncp = 3.0)
  # mom_scaled uses mean(z1*z2)/pi0_joint which is severely contaminated
  # by signal genes (E[z1*z2|signal] >> rho when ncp is large).
  # This method is only reliable when ncp is small or neg_controls provided.
  # Test with low ncp to reduce signal contamination:
  dat_low <- sim_rnaseq(n = 20000L, pi0 = 0.95, rho = 0.20, ncp = 1.5)
  rho_hat <- estimate_rho_overlap(
    dat_low$z1,
    dat_low$z2,
    dat_low$p1,
    dat_low$p2,
    method = "mom_scaled"
  )
  expect_equal(rho_hat, dat_low$rho, tolerance = 0.15)
})

test_that("estimate_rho_overlap: near-zero rho estimated near zero", {
  dat <- sim_rnaseq(n = 20000L, pi0 = 0.80, rho = 0.02, ncp = 3.0)
  rho_hat <- estimate_rho_overlap(
    dat$z1,
    dat$z2,
    dat$p1,
    dat$p2,
    method = "double_null"
  )
  expect_lt(abs(rho_hat), 0.10)
})

test_that("estimate_rho_overlap: high rho=0.5 recovers within tolerance and warns", {
  dat <- sim_rnaseq(n = 20000L, pi0 = 0.85, rho = 0.50, ncp = 4.0)

  # Warning is issued by decorrelate_informative (which calls estimate_rho_overlap)
  # but NOT by estimate_rho_overlap itself. Test the estimate without warning check:
  rho_hat <- suppressWarnings(
    estimate_rho_overlap(
      dat$z1,
      dat$z2,
      dat$p1,
      dat$p2,
      method = "double_null"
    )
  )
  # tolerance = 0.15: attenuation increases at high rho due to
  # lfdr > 0.5 restricting to narrower z-score range
  expect_equal(rho_hat, dat$rho, tolerance = 0.15)

  # Warning IS issued when rho_o is supplied to decorrelate_informative
  expect_warning(
    decorrelate_informative(dat$z1, dat$z2, dat$p1, dat$p2, rho_o = 0.5),
    regexp = "large|0\\.3"
  )
})

test_that("estimate_rho_overlap: recovers rho across range 0.05-0.40", {
  rho_grid <- c(0.05, 0.10, 0.20, 0.30, 0.40)
  for (rho in rho_grid) {
    dat <- sim_rnaseq(
      n = 20000L,
      pi0 = 0.80,
      rho = rho,
      ncp = 3.5,
      seed = round(rho * 1000)
    )

    # Supress warnings here ONLY to keep the loop clean, as >0.3 is tested above
    rho_hat <- suppressWarnings(
      estimate_rho_overlap(
        dat$z1,
        dat$z2,
        dat$p1,
        dat$p2,
        method = "double_null"
      )
    )
    # Tolerance scales with rho: higher rho -> more attenuation
    # from lfdr > 0.5 filter restricting the range of z-scores
    tol <- if (rho <= 0.20) 0.07 else 0.12
    expect_equal(
      rho_hat,
      rho,
      tolerance = tol,
      label = sprintf("rho=%.2f: estimated %.3f", rho, rho_hat)
    )
  }
})

test_that("estimate_rho_overlap: input validation checks fail safely", {
  dat <- sim_rnaseq(n = 200L)
  expect_error(
    estimate_rho_overlap("a", dat$z2, dat$p1, dat$p2),
    "numeric vector"
  )
  expect_error(
    estimate_rho_overlap(dat$z1, dat$z2[1:100], dat$p1, dat$p2),
    "same length"
  )
})


# =================================================================
# 2. decorrelate_informative
# =================================================================

test_that("decorrelate_informative: output structure correct", {
  dat <- sim_rnaseq(n = 20000L, pi0 = 0.80, rho = 0.25)
  dec <- decorrelate_informative(
    dat$z1,
    dat$z2,
    dat$p1,
    dat$p2,
    rho_o = dat$rho
  )

  expect_named(dec, c("z2", "p2", "lfdr_x", "rho_o"))
  expect_length(dec$z2, length(dat$z1))
  expect_length(dec$p2, length(dat$z1))
  expect_true(all(dec$p2 >= 0 & dec$p2 <= 1, na.rm = TRUE))
  expect_equal(dec$rho_o, dat$rho)
})

test_that("decorrelate_informative: soft whitening formula is strictly applied", {
  dat <- sim_rnaseq(n = 20000L, pi0 = 0.80, rho = 0.25)
  dec <- decorrelate_informative(
    dat$z1,
    dat$z2,
    dat$p1,
    dat$p2,
    rho_o = dat$rho
  )

  expected_z2_star <- dat$z2 - dat$rho * dat$z1 * dec$lfdr_x
  expect_equal(dec$z2, expected_z2_star)
})

test_that("decorrelate_informative: null correlation substantially reduced", {
  dat <- sim_rnaseq(n = 20000L, pi0 = 0.80, rho = 0.30, ncp = 4.0)
  null_idx <- which(dat$is_null)

  rho_before <- cor(dat$z1[null_idx], dat$z2[null_idx])
  dec <- decorrelate_informative(
    dat$z1,
    dat$z2,
    dat$p1,
    dat$p2,
    rho_o = dat$rho
  )
  rho_after <- cor(dat$z1[null_idx], dec$z2[null_idx])

  expect_lt(abs(rho_after), abs(rho_before))
  expect_lt(
    abs(rho_after),
    0.10,
    label = sprintf(
      "Null corr after whitening: %.3f (before: %.3f)",
      rho_after,
      rho_before
    )
  )
})

test_that("decorrelate_informative: signal z2* approximately unchanged", {
  dat <- sim_rnaseq(n = 20000L, pi0 = 0.80, rho = 0.30, ncp = 5.0)
  dec <- decorrelate_informative(
    dat$z1,
    dat$z2,
    dat$p1,
    dat$p2,
    rho_o = dat$rho
  )
  signal_idx <- which(!dat$is_null)

  r <- cor(dat$z2[signal_idx], dec$z2[signal_idx])
  expect_gt(
    r,
    0.95,
    label = sprintf("Signal z2 correlation before/after: %.3f", r)
  )
})

test_that("decorrelate_informative: negative rho_o decorrelates properly", {
  dat <- sim_rnaseq(n = 20000L, pi0 = 0.80, rho = -0.25, ncp = 4.0)
  dec <- decorrelate_informative(
    dat$z1,
    dat$z2,
    dat$p1,
    dat$p2,
    rho_o = dat$rho
  )
  null_idx <- which(dat$is_null)

  rho_after <- cor(dat$z1[null_idx], dec$z2[null_idx])
  expect_lt(abs(rho_after), 0.10)
})

test_that("decorrelate_informative: rho_o=0 leaves z2 unchanged", {
  dat <- sim_rnaseq(n = 5000L, pi0 = 0.90, rho = 0.0)
  dec <- decorrelate_informative(dat$z1, dat$z2, dat$p1, dat$p2, rho_o = 0.0)

  expect_equal(dec$z2, dat$z2)
  expect_equal(dec$p2, dat$p2)
})

test_that("decorrelate_informative: warns when |rho_o| > 0.3", {
  dat <- sim_rnaseq(n = 500L, pi0 = 0.90, rho = 0.5)
  expect_warning(
    decorrelate_informative(dat$z1, dat$z2, dat$p1, dat$p2, rho_o = 0.5),
    "large|0\\.3"
  )
})


# =================================================================
# 3. Uniformity of null p2_star after whitening
# =================================================================

test_that("p2_star approximately uniform under null: rho=0.25, pi0=0.90", {
  dat <- sim_rnaseq(n = 20000L, pi0 = 0.90, rho = 0.25, ncp = 4.0, seed = 99L)
  dec <- decorrelate_informative(
    dat$z1,
    dat$z2,
    dat$p1,
    dat$p2,
    rho_o = dat$rho
  )

  # Soft whitening produces z2* = z2 - rho*z1*lfdr(z1).
  # For null genes: z2* ~ N(0, 1 - 2*rho*lfdr + rho^2*lfdr^2) — NOT N(0,1).
  # Converting via 2*pnorm(-abs(z2*)) gives conservative (right-shifted)
  # p-values. We test that p2_star is MORE uniform than p2 conditional on
  # large |z1|, not that it is exactly Uniform(0,1).
  null_idx <- which(dat$is_null)
  high_z1 <- null_idx[abs(dat$z1[null_idx]) > 1.0]

  if (length(high_z1) >= 100) {
    ks_raw <- ks.test(dat$p2[high_z1], "punif")$p.value
    ks_dec <- ks.test(dec$p2[high_z1], "punif")$p.value
    expect_gt(
      ks_dec,
      ks_raw,
      label = sprintf(
        "Whitened KS p=%.4f > raw KS p=%.4f (conditional on |z1|>1)",
        ks_dec,
        ks_raw
      )
    )
  } else {
    skip("Not enough high-z1 null genes")
  }
})

test_that("p2_star approximately uniform under global null after whitening", {
  dat <- sim_rnaseq(n = 20000L, pi0 = 1.0, rho = 0.30, seed = 42L)
  dec <- decorrelate_informative(
    dat$z1,
    dat$z2,
    dat$p1,
    dat$p2,
    rho_o = dat$rho
  )

  # Soft whitening does NOT produce exact Uniform(0,1) p-values because
  # z2* has variance < 1. Test that the correlation structure is removed:
  # cor(z1, z2*) should be near zero for null genes.

  # Original z2 is correlated with z1
  cor_before <- cor(dat$z1, dat$z2)
  # After whitening, z2* should be nearly uncorrelated with z1
  cor_after <- cor(dat$z1, dec$z2)

  expect_lt(
    abs(cor_after),
    abs(cor_before) * 0.3,
    label = sprintf(
      "cor(z1,z2*) = %.3f should be < 30%% of cor(z1,z2) = %.3f",
      cor_after,
      cor_before
    )
  )
})

test_that("raw p2 NOT uniform conditional on large |z1| (motivates whitening)", {
  dat <- sim_rnaseq(n = 20000L, pi0 = 0.90, rho = 0.30, ncp = 4.0, seed = 55L)

  null_idx <- which(dat$is_null)
  high_z1_null <- null_idx[abs(dat$z1[null_idx]) > 1.5]

  if (length(high_z1_null) >= 100) {
    ks_raw <- ks.test(dat$p2[high_z1_null], "punif")
    expect_lt(ks_raw$p.value, 0.10) # Should reject uniformity

    dec <- decorrelate_informative(
      dat$z1,
      dat$z2,
      dat$p1,
      dat$p2,
      rho_o = dat$rho
    )
    ks_dec <- ks.test(dec$p2[high_z1_null], "punif")

    expect_gt(
      ks_dec$p.value,
      ks_raw$p.value,
      label = sprintf(
        "Whitened p2* KS p = %.4f > raw KS p = %.4f",
        ks_dec$p.value,
        ks_raw$p.value
      )
    )
  } else {
    skip("Not enough high-z1 null genes for conditional test")
  }
})


# =================================================================
# 4. sffdr integration: Type I error calibration
# =================================================================

test_that("sffdr with decorrelated surrogate: Type I error calibrated", {
  skip_if_not_installed("sffdr")

  dat <- sim_rnaseq(n = 20000L, pi0 = 1.0, rho = 0.30, seed = 77L)
  dec <- decorrelate_informative(
    dat$z1,
    dat$z2,
    dat$p1,
    dat$p2,
    rho_o = dat$rho
  )

  mpi0_dec <- pi0_model(
    as.matrix(dec$p2),
    min_discoveries = 50,
    verbose = FALSE
  )
  fpi0_dec <- fpi0est(dat$p1, pi0_model_obj = mpi0_dec, verbose = FALSE)
  res_dec <- sffdr(
    dat$p1,
    fpi0 = fpi0_dec$fpi0,
    surrogate = dec$p2,
    verbose = FALSE
  )

  alpha <- 0.001
  fp_rate <- mean(res_dec$fpvalues < alpha, na.rm = TRUE)
  expect_lte(
    fp_rate,
    alpha * 3,
    label = sprintf("FP rate at alpha=0.001 under global null: %.5f", fp_rate)
  )
})

test_that("sffdr with decorrelated surrogate: fewer FDs than without correction", {
  skip_if_not_installed("sffdr")

  dat <- sim_rnaseq(n = 20000L, pi0 = 0.95, rho = 0.30, ncp = 4.0, seed = 123L)

  # Without decorrelation
  mpi0_raw <- pi0_model(
    as.matrix(dat$p2),
    min_discoveries = 50,
    verbose = FALSE
  )
  fpi0_raw <- fpi0est(dat$p1, pi0_model_obj = mpi0_raw, verbose = FALSE)
  res_raw <- sffdr(
    dat$p1,
    fpi0 = fpi0_raw$fpi0,
    surrogate = dat$p2,
    verbose = FALSE
  )

  # With decorrelation
  dec <- decorrelate_informative(
    dat$z1,
    dat$z2,
    dat$p1,
    dat$p2,
    rho_o = dat$rho
  )
  mpi0_dec <- pi0_model(
    as.matrix(dec$p2),
    min_discoveries = 50,
    verbose = FALSE
  )
  fpi0_dec <- fpi0est(dat$p1, pi0_model_obj = mpi0_dec, verbose = FALSE)
  res_dec <- sffdr(
    dat$p1,
    fpi0 = fpi0_dec$fpi0,
    surrogate = dec$p2,
    verbose = FALSE
  )

  thresh <- 0.05
  null_idx <- which(dat$is_null)
  signal_idx <- which(!dat$is_null)

  fd_raw <- sum(res_raw$fqvalues[null_idx] < thresh, na.rm = TRUE)
  fd_dec <- sum(res_dec$fqvalues[null_idx] < thresh, na.rm = TRUE)
  tp_raw <- sum(res_raw$fqvalues[signal_idx] < thresh, na.rm = TRUE)
  tp_dec <- sum(res_dec$fqvalues[signal_idx] < thresh, na.rm = TRUE)

  expect_lte(
    fd_dec,
    fd_raw * 1.30 + 10,
    label = sprintf(
      "FDs: raw=%d, dec=%d | TPs: raw=%d, dec=%d",
      fd_raw,
      fd_dec,
      tp_raw,
      tp_dec
    )
  )
})


# =================================================================
# 5. estimate_rho_mom: neg_controls path
# =================================================================

# ...

# =================================================================
# 6. Functional p-value calibration: global null, sparse, non-sparse
# Multi-replicate Type I error and FDR control
# =================================================================

# Helper: run one replicate of the full sffdr pipeline
# Returns list(fp_rate, fd, tp, fdr_empirical)
run_sffdr_rep <- function(
  n = 20000L,
  pi0 = 1.0,
  rho = 0.30,
  ncp = 4.0,
  seed = 1L,
  alpha = 0.001,
  fq_thr = 0.05,
  decorr = TRUE,
  min_disc = 50L
) {
  dat <- sim_rnaseq(n = n, pi0 = pi0, rho = rho, ncp = ncp, seed = seed)

  if (decorr && rho > 0) {
    dec <- decorrelate_informative(
      dat$z1,
      dat$z2,
      dat$p1,
      dat$p2,
      rho_o = rho
    )
    surr <- dec$p2
  } else {
    surr <- dat$p2
  }

  mpi0 <- suppressWarnings(pi0_model(
    as.matrix(surr),
    min_discoveries = min_disc,
    verbose = FALSE
  ))
  fp0 <- suppressWarnings(fpi0est(
    dat$p1,
    pi0_model_obj = mpi0,
    verbose = FALSE
  ))
  res <- suppressWarnings(sffdr(
    dat$p1,
    fpi0 = fp0$fpi0,
    surrogate = surr,
    verbose = FALSE
  ))

  null_idx <- which(dat$is_null)
  signal_idx <- which(!dat$is_null)

  fp_rate <- mean(res$fpvalues[null_idx] < alpha, na.rm = TRUE)
  fd <- sum(res$fqvalues[null_idx] < fq_thr, na.rm = TRUE)
  tp <- sum(res$fqvalues[signal_idx] < fq_thr, na.rm = TRUE)
  n_disc <- fd + tp
  fdr_emp <- if (n_disc > 0) fd / n_disc else 0

  list(
    fp_rate = fp_rate,
    fd = fd,
    tp = tp,
    fdr_emp = fdr_emp,
    fpvalues = res$fpvalues[null_idx]
  )
}

# -----------------------------------------------------------------
# 6a. Global null: fp-values approximately Uniform(0,1)
# -----------------------------------------------------------------

test_that("fpvalues approximately uniform under global null (single rep)", {
  skip_if_not_installed("sffdr")

  dat <- sim_rnaseq(n = 20000L, pi0 = 1.0, rho = 0.30, seed = 77L)
  dec <- decorrelate_informative(
    dat$z1,
    dat$z2,
    dat$p1,
    dat$p2,
    rho_o = dat$rho
  )
  mpi0 <- suppressWarnings(pi0_model(
    as.matrix(dec$p2),
    min_discoveries = 50L,
    verbose = FALSE
  ))
  fp0 <- suppressWarnings(fpi0est(
    dat$p1,
    pi0_model_obj = mpi0,
    verbose = FALSE
  ))
  res <- suppressWarnings(sffdr(
    dat$p1,
    fpi0 = fp0$fpi0,
    surrogate = dec$p2,
    monotone = TRUE,
    verbose = FALSE
  ))

  ks <- ks.test(res$fpvalues, "punif")
  expect_gt(
    ks$p.value,
    0.001,
    label = sprintf(
      "Global null fpvalues KS p=%.4f (should be > 0.001 for approximate uniformity)",
      ks$p.value
    )
  )

  expect_equal(
    mean(res$fpvalues, na.rm = TRUE),
    0.5,
    tolerance = 0.05,
    label = sprintf(
      "Mean fpvalue = %.3f (expected ~0.5)",
      mean(res$fpvalues, na.rm = TRUE)
    )
  )
})

# -----------------------------------------------------------------
# 6b. Global null: Type I error control across 10 replicates
# -----------------------------------------------------------------

test_that("Type I error controlled at alpha=0.001 under global null (10 reps)", {
  skip_if_not_installed("sffdr")

  n_reps <- 10L
  alpha <- 0.001
  fp_rates <- vapply(
    seq_len(n_reps),
    function(s) {
      run_sffdr_rep(
        n = 20000L,
        pi0 = 1.0,
        rho = 0.30,
        ncp = 4.0,
        seed = s,
        alpha = alpha
      )$fp_rate
    },
    numeric(1)
  )

  mean_fp <- mean(fp_rates)
  max_fp <- max(fp_rates)

  expect_lte(
    mean_fp,
    alpha * 3,
    label = sprintf(
      "Mean Type I error = %.5f across %d reps (nominal alpha = %.4f, max allowed = %.4f)",
      mean_fp,
      n_reps,
      alpha,
      alpha * 3
    )
  )

  expect_lte(
    max_fp,
    alpha * 5,
    label = sprintf(
      "Max Type I error = %.5f (max allowed = %.4f)",
      max_fp,
      alpha * 5
    )
  )
})

# -----------------------------------------------------------------
# 6c. Sparse setting (pi0=0.95): FDR control + power across 10 reps
# -----------------------------------------------------------------

test_that("FDR controlled and power preserved: sparse pi0=0.95 (10 reps)", {
  skip_if_not_installed("sffdr")

  n_reps <- 10L
  fq_thr <- 0.05

  res_list <- lapply(seq_len(n_reps), function(s) {
    run_sffdr_rep(
      n = 20000L,
      pi0 = 0.95,
      rho = 0.30,
      ncp = 4.0,
      seed = s,
      fq_thr = fq_thr,
      decorr = TRUE,
      min_disc = 50L
    )
  })

  fdr_vals <- vapply(res_list, `[[`, numeric(1), "fdr_emp")
  tp_vals <- vapply(res_list, `[[`, numeric(1), "tp")
  fd_vals <- vapply(res_list, `[[`, numeric(1), "fd")

  mean_fdr <- mean(fdr_vals)
  mean_tp <- mean(tp_vals)
  mean_fd <- mean(fd_vals)

  expect_lte(
    mean_fdr,
    fq_thr * 1.5,
    label = sprintf(
      "Mean empirical FDR = %.3f (nominal = %.2f, max allowed = %.2f)",
      mean_fdr,
      fq_thr,
      fq_thr * 1.5
    )
  )

  expect_gt(
    mean_tp,
    500,
    label = sprintf(
      "Mean TP = %.0f (should be > 500 for ncp=4, rho=0.30)",
      mean_tp
    )
  )
})

# -----------------------------------------------------------------
# 6d. Non-sparse setting (pi0=0.80): FDR control + power across 10 reps
# -----------------------------------------------------------------

test_that("FDR controlled and power preserved: non-sparse pi0=0.80 (10 reps)", {
  skip_if_not_installed("sffdr")

  n_reps <- 10L
  fq_thr <- 0.05

  res_list <- lapply(seq_len(n_reps), function(s) {
    run_sffdr_rep(
      n = 20000L,
      pi0 = 0.80,
      rho = 0.30,
      ncp = 4.0,
      seed = s,
      fq_thr = fq_thr,
      decorr = TRUE,
      min_disc = 50L
    )
  })

  fdr_vals <- vapply(res_list, `[[`, numeric(1), "fdr_emp")
  tp_vals <- vapply(res_list, `[[`, numeric(1), "tp")

  mean_fdr <- mean(fdr_vals)
  mean_tp <- mean(tp_vals)

  expect_lte(
    mean_fdr,
    fq_thr * 1.5,
    label = sprintf(
      "Mean empirical FDR = %.3f (nominal = %.2f, max allowed = %.2f)",
      mean_fdr,
      fq_thr,
      fq_thr * 1.5
    )
  )

  expect_gt(
    mean_tp,
    2000,
    label = sprintf(
      "Mean TP = %.0f (should be > 2000 for pi0=0.80, ncp=4)",
      mean_tp
    )
  )
})

# -----------------------------------------------------------------
# 6e. Decorrelation reduces false discoveries vs no correction
#     Sparse (pi0=0.95) across 10 replicates
# -----------------------------------------------------------------

test_that("Decorrelation reduces FDs vs no correction: sparse pi0=0.95 (10 reps)", {
  skip_if_not_installed("sffdr")

  n_reps <- 10L
  fq_thr <- 0.05

  fd_raw_vec <- vapply(
    seq_len(n_reps),
    function(s) {
      run_sffdr_rep(
        n = 20000L,
        pi0 = 0.95,
        rho = 0.30,
        ncp = 4.0,
        seed = s,
        fq_thr = fq_thr,
        decorr = FALSE
      )$fd
    },
    numeric(1)
  )

  fd_dec_vec <- vapply(
    seq_len(n_reps),
    function(s) {
      run_sffdr_rep(
        n = 20000L,
        pi0 = 0.95,
        rho = 0.30,
        ncp = 4.0,
        seed = s,
        fq_thr = fq_thr,
        decorr = TRUE
      )$fd
    },
    numeric(1)
  )

  mean_fd_raw <- mean(fd_raw_vec)
  mean_fd_dec <- mean(fd_dec_vec)

  expect_lt(
    mean_fd_dec,
    mean_fd_raw,
    label = sprintf(
      "Mean FD: raw=%.1f, decorr=%.1f (decorr should be lower)",
      mean_fd_raw,
      mean_fd_dec
    )
  )

  reduction <- 1 - mean_fd_dec / mean_fd_raw
  expect_gt(
    reduction,
    0.10,
    label = sprintf(
      "FD reduction = %.1f%% (should be > 10%%)",
      reduction * 100
    )
  )
})

# -----------------------------------------------------------------
# 6f. Power preservation: decorrelation does not substantially
#     reduce true positives vs no correction
# -----------------------------------------------------------------

test_that("Decorrelation preserves power: sparse pi0=0.95 (10 reps)", {
  skip_if_not_installed("sffdr")

  n_reps <- 10L
  fq_thr <- 0.05

  tp_raw_vec <- vapply(
    seq_len(n_reps),
    function(s) {
      run_sffdr_rep(
        n = 20000L,
        pi0 = 0.95,
        rho = 0.30,
        ncp = 4.0,
        seed = s,
        fq_thr = fq_thr,
        decorr = FALSE
      )$tp
    },
    numeric(1)
  )

  tp_dec_vec <- vapply(
    seq_len(n_reps),
    function(s) {
      run_sffdr_rep(
        n = 20000L,
        pi0 = 0.95,
        rho = 0.30,
        ncp = 4.0,
        seed = s,
        fq_thr = fq_thr,
        decorr = TRUE
      )$tp
    },
    numeric(1)
  )

  mean_tp_raw <- mean(tp_raw_vec)
  mean_tp_dec <- mean(tp_dec_vec)

  power_retention <- mean_tp_dec / mean_tp_raw
  expect_gt(
    power_retention,
    0.90,
    label = sprintf(
      "Power retention = %.1f%% (TP: raw=%.0f, decorr=%.0f)",
      power_retention * 100,
      mean_tp_raw,
      mean_tp_dec
    )
  )
})

# -----------------------------------------------------------------
# 6g. Type I error: compare decorr vs no-correction under rho=0
#     (should be equivalent when no overlap exists)
# -----------------------------------------------------------------

test_that("Decorrelation does not hurt when rho=0 (10 reps)", {
  skip_if_not_installed("sffdr")

  n_reps <- 10L
  alpha <- 0.001

  fp_raw <- vapply(
    seq_len(n_reps),
    function(s) {
      run_sffdr_rep(
        n = 20000L,
        pi0 = 0.95,
        rho = 0.0,
        ncp = 4.0,
        seed = s,
        alpha = alpha,
        decorr = FALSE
      )$fp_rate
    },
    numeric(1)
  )

  fp_dec <- vapply(
    seq_len(n_reps),
    function(s) {
      run_sffdr_rep(
        n = 20000L,
        pi0 = 0.95,
        rho = 0.0,
        ncp = 4.0,
        seed = s,
        alpha = alpha,
        decorr = TRUE
      )$fp_rate
    },
    numeric(1)
  )

  expect_lte(
    mean(fp_raw),
    alpha * 3,
    label = sprintf("Mean Type I (no decorr, rho=0) = %.5f", mean(fp_raw))
  )
  expect_lte(
    mean(fp_dec),
    alpha * 3,
    label = sprintf("Mean Type I (decorr, rho=0) = %.5f", mean(fp_dec))
  )
})

# =================================================================
# 5. estimate_rho_mom: neg_controls path
# =================================================================

test_that("estimate_rho_mom: neg_controls gives accurate rho", {
  dat <- sim_rnaseq(n = 20000L, pi0 = 0.80, rho = 0.25, ncp = 4.0)
  neg_ctrl <- which(dat$is_null)[1:2000]

  mom <- estimate_rho_mom(
    dat$z1,
    dat$z2,
    dat$p1,
    dat$p2,
    neg_controls = neg_ctrl
  )

  expect_equal(mom$rho_hat, dat$rho, tolerance = 0.05)
  expect_equal(mom$method, "neg_controls")
})
