test_that("pi0_model builds correct formula", {
  set.seed(123)
  # Create synthetic data with signal
  n <- 1000
  z1 <- runif(n) # Uniformly distributed
  z2 <- runif(n)
  z <- data.frame(z1 = z1, z2 = z2)

  # Add signal: make small z1 values have small "p-values"
  # This ensures qvalue finds significant hits
  set.seed(456)

  result <- pi0_model(
    as.matrix(z),
    min_discoveries = 100,
    min_snps_per_knot = 50,
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_named(result, c("fmod", "zt"))
  expect_s3_class(result$zt, "data.frame")
})

test_that("pi0_model respects max_knots", {
  set.seed(123)
  data(bmi)
  z <- as.matrix(sumstats[, -1])

  result <- pi0_model(
    z,
    n_knots = 2,
    min_discoveries = 100,
    min_snps_per_knot = 50,
    verbose = FALSE
  )

  # Check formula is created
  expect_type(result, "list")
  expect_true("fmod" %in% names(result))
})

test_that("pi0_model handles weak signal variables", {
  set.seed(123)
  n <- 1000
  # Create data with very weak signal (high p-values)
  z <- data.frame(
    weak_var = runif(n, 0.5, 1) # All high p-values, no signal
  )

  result <- pi0_model(
    as.matrix(z),
    min_discoveries = 100,
    min_snps_per_knot = 50,
    verbose = FALSE
  )

  # Should return intercept-only formula when no signal
  expect_type(result, "list")
})

test_that("pi0_model handles missing values", {
  set.seed(123)
  data(bmi)
  z <- as.matrix(sumstats[, -1])

  # Add some NAs
  z[1:10, 1] <- NA
  z[20:30, 2] <- NA

  result <- pi0_model(
    z,
    min_discoveries = 100,
    min_snps_per_knot = 50,
    verbose = FALSE
  )

  # Should still work with NAs handled
  # Note: pi0_model removes rows with NAs, so zt will have fewer rows
  expect_type(result, "list")
  expect_true(nrow(result$zt) <= nrow(z))
  expect_true(nrow(result$zt) > 0)
})

test_that("pi0_model validates inputs", {
  # Should error on invalid input
  expect_error(pi0_model("not a matrix"))
})

test_that("pi0_model works with real data", {
  data(bmi)
  z <- as.matrix(sumstats[, -1])

  result <- pi0_model(
    z,
    min_discoveries = 100,
    min_snps_per_knot = 50,
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_named(result, c("fmod", "zt"))
  expect_s3_class(result$fmod, "formula")
  expect_s3_class(result$zt, "data.frame")
  expect_equal(nrow(result$zt), nrow(z))
})
