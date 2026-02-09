test_that("kernelEstimator works with 1D p-values", {
  set.seed(123)
  p <- runif(1000, 0, 1)

  result <- kernelEstimator(p)

  expect_s3_class(result, "data.frame")
  expect_named(result, c("x", "fx", "s", "fs"))
  expect_equal(nrow(result), length(p))
  expect_true(all(result$fx >= 0))
  expect_true(all(is.finite(result$fx)))
})

test_that("kernelEstimator works with 2D (p-values + covariate)", {
  set.seed(123)
  n <- 5000
  # Create realistic GWAS-like data with signal
  z <- c(rnorm(4000, 0, 1), rnorm(1000, -3, 0.5)) # 20% signal in tail
  p <- c(runif(4000, 0.1, 1), runif(1000, 0, 0.01)) # Signal has small p-values
  x_mat <- cbind(p, z)

  result <- kernelEstimator(x_mat)

  expect_s3_class(result, "data.frame")
  expect_true("fx" %in% names(result))
  expect_true("fs" %in% names(result))
  expect_equal(nrow(result), n)
  expect_true(all(result$fx >= 0))
  expect_true(all(is.finite(result$fx)))
})

test_that("kernelEstimator handles small p-values correctly", {
  set.seed(123)
  p <- c(c(1e-200, 1e-50, 1e-10, 0.001, 0.5, 0.9), runif(50))

  # Suppress expected locfit warnings for extreme values
  suppressWarnings({
    result <- kernelEstimator(p)
  })

  expect_s3_class(result, "data.frame")
  expect_true(all(is.finite(result$fx)))
  expect_true(all(result$fx >= 0))
})

test_that("kernelEstimator respects epsilon bounds", {
  set.seed(123)
  p <- c(c(0, 1e-320, 0.5, 0.999999, 1.0), runif(50))

  suppressWarnings({
    result <- kernelEstimator(p, epsilon = 1e-15, epsilon.max = 1 - 1e-4)
  })

  expect_true(all(is.finite(result$s)))
  expect_true(all(is.finite(result$fx)))
})

test_that("kernelEstimator handles custom evaluation points", {
  set.seed(123)
  p <- runif(1000, 0, 1)
  eval_pts <- seq(0.01, 0.99, by = 0.01)

  result <- kernelEstimator(p, eval.points = eval_pts)

  expect_equal(nrow(result), length(eval_pts))
})

test_that("kernelEstimator handles custom nn parameter", {
  set.seed(123)
  p <- runif(1000, 0, 1)

  result1 <- kernelEstimator(p, nn = 0.05)
  result2 <- kernelEstimator(p, nn = 0.1)

  # Different nn should give different results
  expect_false(identical(result1$fx, result2$fx))
})

test_that("kernelEstimator downsampling works correctly", {
  set.seed(123)
  n <- 10000
  # Create data with clear tail structure
  z <- c(rnorm(8000, 0, 1), rnorm(2000, -3, 0.5))
  p <- c(runif(8000, 0.1, 1), runif(2000, 0, 0.01))
  x_mat <- cbind(p, z)

  # Use smaller target to force downsampling
  result <- kernelEstimator(x_mat, target_null = 5000)

  expect_s3_class(result, "data.frame")
  expect_true(all(is.finite(result$fx)))
})

test_that("kernelEstimator tail prioritization works", {
  set.seed(123)
  n <- 10000
  # Create data with strong signal in tails
  z <- c(rnorm(8000, 0, 1), rnorm(2000, -3, 0.5))
  p <- c(runif(8000, 0.1, 1), runif(2000, 0, 0.01))
  x_mat <- cbind(p, z)

  result <- kernelEstimator(x_mat, tail_threshold = -2)

  expect_s3_class(result, "data.frame")
  expect_true(all(is.finite(result$fx)))
})

test_that("kernelEstimator returns lfit attribute", {
  set.seed(123)
  p <- runif(1000, 0, 1)

  result <- kernelEstimator(p)

  expect_true(!is.null(attr(result, "lfit")))
  expect_s3_class(attr(result, "lfit"), "locfit")
})

test_that("kernelEstimator trim parameter works", {
  set.seed(123)
  p <- runif(1000, 0, 1)

  result <- kernelEstimator(p, trim = 0.05)

  # Values below 0.05 should be constant
  fx_below <- result$fx[result$x < 0.05]
  expect_true(length(unique(fx_below)) <= 2) # Allow for some numerical noise

  # Values above 0.95 should be constant
  fx_above <- result$fx[result$x > 0.95]
  expect_true(length(unique(fx_above)) <= 2)
})

test_that("kernelEstimator handles matrix with wrong dimensions", {
  set.seed(123)
  x_wrong <- matrix(runif(100), ncol = 1)

  expect_error(kernelEstimator(x_wrong), "at least 2 columns")
})

test_that("kernelEstimator handles all identical values", {
  p <- rep(0.5, 100)

  suppressWarnings({
    result <- kernelEstimator(p)
  })

  expect_s3_class(result, "data.frame")
  expect_true(all(is.finite(result$fx)))
})

test_that("kernelEstimator handles very small sample sizes", {
  set.seed(123)
  p <- runif(50, 0, 1)

  result <- kernelEstimator(p)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 50)
  expect_true(all(is.finite(result$fx)))
})

test_that("kernelEstimator reproducibility with set.seed", {
  p <- runif(1000, 0, 1)

  set.seed(42)
  result1 <- kernelEstimator(p)

  set.seed(42)
  result2 <- kernelEstimator(p)

  expect_equal(result1$fx, result2$fx)
})

test_that("kernelEstimator handles maxk parameter", {
  set.seed(123)
  p <- runif(1000, 0, 1)

  # Should work with custom maxk
  result <- kernelEstimator(p, maxk = 100000)

  expect_s3_class(result, "data.frame")
  expect_true(all(is.finite(result$fx)))
})

test_that("kernelEstimator adaptive jitter works", {
  set.seed(123)
  # Create p-values with known resolution
  p <- seq(0.001, 0.999, by = 0.001)

  result <- kernelEstimator(p)

  expect_s3_class(result, "data.frame")
  expect_true(all(is.finite(result$fx)))
})

test_that("kernelEstimator cascade fallback works", {
  set.seed(123)
  # Create realistic GWAS-like data
  n <- 5000
  z <- c(rnorm(4000, 0, 1), rnorm(1000, -3, 0.5))
  p <- c(runif(4000, 0.1, 1), runif(1000, 0, 0.01))
  x_mat <- cbind(p, z)

  # Should fall back through strategies if needed
  result <- kernelEstimator(x_mat)

  expect_s3_class(result, "data.frame")
  expect_true(all(is.finite(result$fx)))
})
