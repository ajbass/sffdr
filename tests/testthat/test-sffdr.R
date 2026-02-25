test_that("sffdr produces valid functional p-values", {
  set.seed(123)
  n <- 1000
  p <- runif(n, 0, 1)
  fpi0 <- runif(n, 0.7, 0.95)

  result <- sffdr(p, fpi0)

  expect_s3_class(result, "sffdr")
  expect_named(
    result,
    c("call", "pvalues", "fpvalues", "fqvalues", "flfdr", "fpi0", "fx")
  )
  expect_true(all(result$fpvalues >= 0 & result$fpvalues <= 1))
  expect_true(all(result$fqvalues >= 0 & result$fqvalues <= 1))
  expect_true(all(result$flfdr >= 0 & result$flfdr <= 1))
  expect_equal(length(result$fpvalues), n)
})

test_that("sffdr handles custom surrogate variable", {
  set.seed(123)
  n <- 500
  p <- runif(n, 0, 1)
  fpi0 <- runif(n, 0.8, 0.95)
  surrogate <- rnorm(n)

  result <- sffdr(p, fpi0, surrogate = surrogate)

  expect_s3_class(result, "sffdr")
  expect_true(all(result$fpvalues >= 0 & result$fpvalues <= 1))
})

test_that("sffdr validates input lengths", {
  p <- runif(100, 0, 1)
  fpi0_wrong <- runif(50, 0.8, 0.95)

  expect_error(sffdr(p, fpi0_wrong), "fpi0")
})

test_that("sffdr validates fpi0 range", {
  p <- runif(100, 0, 1)
  fpi0_invalid <- runif(100, 1.1, 1.5)

  expect_error(sffdr(p, fpi0_invalid), "\\[0, 1\\]")
})

test_that("sffdr handles zero p-values", {
  set.seed(123)
  n <- 100
  p <- c(0, 0, runif(98, 0, 1))
  fpi0 <- runif(n, 0.8, 0.95)

  expect_warning(result <- sffdr(p, fpi0), "p-values <= 0")
  expect_s3_class(result, "sffdr")
  expect_true(all(result$fpvalues >= 0 & result$fpvalues <= 1))
})

test_that("sffdr with seed produces reproducible results", {
  set.seed(456)
  n <- 100
  p <- runif(n, 0, 1)
  fpi0 <- runif(n, 0.8, 0.95)

  result1 <- sffdr(p, fpi0, seed = 42)
  result2 <- sffdr(p, fpi0, seed = 42)

  expect_equal(result1$fpvalues, result2$fpvalues)
  expect_equal(result1$fqvalues, result2$fqvalues)
})

test_that("sffdr custom nn parameter", {
  set.seed(123)
  n <- 500
  p <- runif(n, 0, 1)
  fpi0 <- runif(n, 0.8, 0.95)

  result <- sffdr(p, fpi0, nn = 0.05)

  expect_s3_class(result, "sffdr")
  expect_true(all(result$fpvalues >= 0 & result$fpvalues <= 1))
})

# Helper functions
simulate_mesa_data <- function(n, n_mesa = floor(0.2 * n)) {
  set.seed(2026)

  p_vals <- runif(n)

  n_signal <- max(10, floor(n * 0.05))
  p_vals[1:n_signal] <- rbeta(n_signal, 1, 150)

  mesa_idx <- sample((n_signal + 1):n, n_mesa)
  p_vals[mesa_idx] <- runif(n_mesa, 0.75, 0.85)

  # Clamp bounds safely
  p_vals <- pmax(pmin(p_vals, 0.999), 0.0001)
  surrogate <- runif(n)
  fpi0 <- rep(0.95, n)

  list(
    p = p_vals,
    z = surrogate,
    fpi0 = fpi0,
    mesa_idx = mesa_idx
  )
}

simulate_mountain_data <- function(n, mountain_p = 0.8, n_mountain = 0.1 * n) {
  set.seed(2026)
  p_vals <- runif(n)

  n_signal <- max(10, floor(n * 0.01))
  p_vals[1:n_signal] <- rbeta(n_signal, 1, 50)

  mountain_idx <- sample((n_signal + 1):n, n_mountain)
  p_vals[mountain_idx] <- rnorm(n_mountain, mean = mountain_p, sd = 0.05)

  p_vals <- pmax(pmin(p_vals, 1), .Machine$double.xmin)

  surrogate <- runif(n)
  fpi0 <- rep(0.95, n)

  list(
    p = p_vals,
    z = surrogate,
    fpi0 = fpi0,
    mountain_idx = mountain_idx,
    clean_idx = setdiff(1:n, mountain_idx)
  )
}

test_that("C++ Core Engine: Conditional smoothing enforces strict running minimums within groups", {
  p_cond <- c(0.1, 0.3, 0.2, 0.4)
  den_cond <- c(10, 12, 8, 4) # Peak of 12 in group 1
  grp_cond <- c(1, 1, 2, 2)

  res_cond <- monoSmooth_conditional(p_cond, den_cond, grp_cond)

  expect_equal(as.numeric(res_cond), c(10, 10, 8, 4))
})

test_that("sffdr Conditional Smoothing flattens local LD mountains without cross-contamination", {
  dat <- simulate_mountain_data(n = 5000)

  # Run without monotonicity
  res_std <- sffdr(
    dat$p,
    dat$fpi0,
    dat$z,
    nn = 0.5,
    monotone = FALSE,
    verbose = FALSE
  )

  # Run with monotonicity
  res_mono <- sffdr(
    dat$p,
    dat$fpi0,
    dat$z,
    nn = 0.5,
    monotone = TRUE,
    verbose = FALSE
  )

  idx_mountain <- dat$mountain_idx
  idx_clean <- which(abs(dat$p - 0.5) < 0.05)

  # Prove the standard density estimator gets tricked by the mountain
  expect_gt(mean(res_std$fx[idx_mountain]), mean(res_std$fx[idx_clean]) * 1.1)

  # Prove the monotonicity flag conditionally crushes the mountain
  expect_lte(
    mean(res_mono$fx[idx_mountain]),
    mean(res_mono$fx[idx_clean]) + 1e-6
  )
})

test_that("sffdr Conditional Smoothing scales efficiently to large datasets (N > 50000)", {
  dat <- simulate_mesa_data(n = 55000)

  res_mono <- sffdr(
    dat$p,
    dat$fpi0,
    dat$z,
    monotone = TRUE,
    verbose = FALSE
  )

  idx_mountain <- dat$mesa_idx
  idx_clean <- which(abs(dat$p - 0.5) < 0.05)

  # Assert the mesa was crushed
  expect_lte(
    mean(res_mono$fx[idx_mountain]),
    mean(res_mono$fx[idx_clean]) + 1e-6
  )

  # Assert the true biological signal was preserved (not crushed by a global minimum)
  idx_signal <- which(dat$p < 0.001)
  expect_gt(mean(res_mono$fx[idx_signal]), mean(res_mono$fx[idx_clean]))

  # Assert no NAs were introduced during mapping
  expect_false(any(is.na(res_mono$fx)))
})

test_that("sffdr gracefully handles data with ZERO true signals (Pure Null)", {
  set.seed(2026)
  n <- 55000

  # Simulate a purely flat, disappointing GWAS with NO signals near 0
  p_vals <- runif(n, 0.0, 1.0)
  surrogate <- runif(n)
  fpi0 <- rep(1.0, n)

  # 1. Assert that the function runs completely without throwing an error
  res_mono <- expect_error(
    sffdr(
      p_vals,
      fpi0,
      surrogate,
      monotone = TRUE,
      verbose = FALSE
    ),
    NA
  )

  # 2. Assert that the density wasn't accidentally crushed to 0
  expect_gt(mean(res_mono$fx), 0.5)

  # 3. Assert no missing values
  expect_false(any(is.na(res_mono$fx)))
})
