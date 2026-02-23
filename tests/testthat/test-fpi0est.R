test_that("fpi0est returns valid pi0 estimates", {
  data(bmi)
  p <- sumstats$bmi
  z <- as.matrix(sumstats[, -1])

  mpi0 <- pi0_model(
    z,
    min_discoveries = 100,
    min_snps_per_knot = 50,
    verbose = FALSE
  )
  result <- fpi0est(p, z = mpi0$zt, pi0_model = mpi0$fmod, verbose = FALSE)

  # Check return structure matches actual output
  expect_named(result, c("fpi0", "tableLambda", "MISE", "lambda"))
  expect_true(all(result$fpi0 >= 0 & result$fpi0 <= 1))
  expect_true(result$lambda > 0 & result$lambda < 1)
})

test_that("fpi0est works with custom lambda values", {
  data(bmi)
  p <- sumstats$bmi
  z <- as.matrix(sumstats[, -1])

  mpi0 <- pi0_model(
    z,
    min_discoveries = 100,
    min_snps_per_knot = 50,
    verbose = FALSE
  )
  custom_lambda <- seq(0.1, 0.9, 0.1)
  result <- fpi0est(
    p,
    z = mpi0$zt,
    pi0_model = mpi0$fmod,
    lambda = custom_lambda,
    verbose = FALSE
  )

  expect_true(result$lambda %in% custom_lambda)
})

test_that("fpi0est uses constrained binomial correctly", {
  data(bmi)
  p <- sumstats$bmi
  z <- as.matrix(sumstats[, -1])

  mpi0 <- pi0_model(
    z,
    min_discoveries = 100,
    min_snps_per_knot = 50,
    verbose = FALSE
  )
  result <- fpi0est(
    p,
    z = mpi0$zt,
    pi0_model = mpi0$fmod,
    constrained.p = TRUE,
    verbose = FALSE
  )

  # All fpi0 estimates should be strictly between 0 and 1
  expect_true(all(result$fpi0 > 0 & result$fpi0 < 1))
})

test_that("fpi0est handles indep_snps correctly", {
  data(bmi)
  # Independent SNPs
  p <- sumstats$bmi
  z <- as.matrix(sumstats[, -1])
  mpi0 <- pi0_model(
    z,
    min_discoveries = 100,
    min_snps_per_knot = 50,
    verbose = FALSE
  )
  result1 <- fpi0est(
    p,
    z = mpi0$zt,
    pi0_model = mpi0$fmod,
    constrained.p = TRUE,
    verbose = FALSE
  )

  # Add dependent SNPs
  p <- c(sumstats$bmi, runif(100))
  z <- rbind(as.matrix(sumstats[, -1]), as.matrix(sumstats[1:100, -1]))
  indep_snps <- c(rep(TRUE, 10000), rep(FALSE, 100))
  mpi0 <- pi0_model(
    z,
    min_discoveries = 100,
    min_snps_per_knot = 50,
    indep_snps = indep_snps,
    verbose = FALSE
  )
  result <- fpi0est(
    p,
    z = mpi0$zt,
    pi0_model = mpi0$fmod,
    indep_snps = indep_snps,
    verbose = FALSE
  )

  # All fpi0 estimates should be strictly between 0 and 1
  expect_equal(result$fpi0[1:10000], result1$fpi0[1:10000], tolerance = 1e-3)
})

test_that("fpi0est handles small datasets", {
  set.seed(123)
  n <- 1000
  p <- runif(n)
  z <- data.frame(
    var1 = runif(n),
    var2 = runif(n)
  )

  # Simple formula without splines for small test
  result <- fpi0est(p, z = z, pi0_model = ~ var1 + var2, verbose = FALSE)

  expect_true(all(result$fpi0 >= 0 & result$fpi0 <= 1))
})

test_weighted_pi0_fixed <- function() {
  set.seed(42)

  n_indep <- 10000
  is_signal <- runif(n_indep) < 0.10
  p_indep <- numeric(n_indep)
  p_indep[!is_signal] <- runif(sum(!is_signal), 0, 1) # Proper uniform nulls
  p_indep[is_signal] <- runif(sum(is_signal), 1e-10, 1e-3)
  w_indep <- rep(1.0, n_indep)

  n_ld <- 40000
  p_ld <- runif(n_ld, 0, 1)
  w_ld <- rep(1 / n_ld, n_ld)

  p_all <- c(p_indep, p_ld)
  w_all <- c(w_indep, w_ld)

  res_unweighted <- pi0est_weighted(
    p_all,
    weights = NULL,
    pi0.method = "bootstrap"
  )
  res_weighted <- pi0est_weighted(
    p_all,
    weights = w_all,
    pi0.method = "bootstrap"
  )

  message(
    "Unweighted pi0 Estimate: ",
    round(res_unweighted$pi0, 4),
    " (Inflated by the 40k LD Nulls)"
  )
  message(
    "Weighted pi0 Estimate:   ",
    round(res_weighted$pi0, 4),
    " (Recovers the 0.90 background!)"
  )

  if (abs(res_weighted$pi0 - 0.90) < 0.05 && res_unweighted$pi0 > 0.97) {
    message(
      "PASS: The inverse-LD weights successfully neutralized the LD block!"
    )
  } else {
    message("FAIL")
  }
}

test_fpi0est_weights <- function() {
  message("Generating synthetic GWAS data with a covariate-linked LD block...")
  set.seed(123)

  # --- 1. Simulate the Independent Biological Background ---
  n_indep <- 10000
  z_indep <- runif(n_indep, 0, 1) # Covariate from 0 to 1

  # True biology: Higher Z means more signal (pi0 drops from 0.9 to 0.4)
  true_pi0_indep <- 0.90 - (0.50 * z_indep)
  is_null_indep <- runif(n_indep) < true_pi0_indep

  p_indep <- numeric(n_indep)
  p_indep[is_null_indep] <- runif(sum(is_null_indep), 0, 1) # Nulls
  p_indep[!is_null_indep] <- rbeta(sum(!is_null_indep), 1, 20) # Strong Signals
  w_indep <- rep(1.0, n_indep)

  # --- 2. Simulate the Massive "LD Bomb" ---
  n_ld <- 40000
  z_ld <- rnorm(n_ld, mean = 0.2, sd = 0.01)
  p_ld <- runif(n_ld, 0, 1) # 100% Null
  w_ld <- rep(1 / n_ld, n_ld) # Inverse LD weights sum to 1

  # --- 3. Combine and Format ---
  p_all <- c(p_indep, p_ld)
  z_raw <- c(z_indep, z_ld)
  w_all <- c(w_indep, w_ld)

  z_ranked <- rank(z_raw, ties.method = "random") / length(z_raw)

  z_df <- data.frame(Z = z_ranked)

  fm <- formula("~ splines::ns(Z, df=3)")

  # --- 4. Run UNWEIGHTED fpi0est ---
  message("Running UNWEIGHTED fpi0est...")
  res_unw <- fpi0est(
    p = p_all,
    z = z_df,
    pi0_model = fm,
    weights = NULL,
    verbose = FALSE
  )

  # --- 5. Run WEIGHTED fpi0est ---
  message("Running WEIGHTED fpi0est...")
  res_w <- fpi0est(
    p = p_all,
    z = z_df,
    pi0_model = fm,
    weights = w_all,
    verbose = FALSE
  )

  # --- 6. Evaluate Accuracy ---
  fpi0_est_unw <- res_unw$fpi0[1:10000]
  fpi0_est_w <- res_w$fpi0[1:10000]

  mae_unw <- mean(abs(fpi0_est_unw - true_pi0_indep))
  mae_w <- mean(abs(fpi0_est_w - true_pi0_indep))

  message("\n--- Test Results ---")
  message(sprintf(
    "Unweighted Mean Abs Error: %.4f (Spline warped by LD block)",
    mae_unw
  ))
  message(sprintf(
    "Weighted Mean Abs Error:   %.4f (Spline recovered true biology)",
    mae_w
  ))

  if (mae_w < 0.05 && mae_unw > 0.075) {
    message(
      "\nPASS: fastglm correctly used the weights in the likelihood to ignore the LD bomb!"
    )
  } else {
    message("\nFAIL: The weights did not protect the spline.")
  }
}
