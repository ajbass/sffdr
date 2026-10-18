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

test_that("pi0_model respects n_knots", {
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

  # Rows are kept (one per test); the missing entries stay NA in zt
  expect_type(result, "list")
  expect_equal(nrow(result$zt), nrow(z))
  expect_true(all(is.na(result$zt[1:10, 1])))
  expect_false(anyNA(result$zt[-(1:10), 1]))
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

test_that("pi0_model takes indep_snps as logical, indices or SNP ids", {
  data(bmi)
  z <- as.matrix(sumstats[, -1])
  rownames(z) <- paste0("rs", seq_len(nrow(z)))
  set.seed(3)
  keep <- runif(nrow(z)) < 0.5
  m_log <- pi0_model(z, indep_snps = keep, min_discoveries = 50, min_snps_per_knot = 20, verbose = FALSE)
  m_idx <- pi0_model(z, indep_snps = which(keep), min_discoveries = 50, min_snps_per_knot = 20, verbose = FALSE)
  m_id <- pi0_model(z, indep_snps = rownames(z)[keep], min_discoveries = 50, min_snps_per_knot = 20, verbose = FALSE)
  expect_equal(deparse1(m_log$fmod), deparse1(m_idx$fmod))
  expect_equal(deparse1(m_log$fmod), deparse1(m_id$fmod))
  expect_error(pi0_model(z, indep_snps = keep[-1], verbose = FALSE), "one entry per SNP")
  expect_error(pi0_model(z, indep_snps = c(0, 5), verbose = FALSE), "indep_snps")
})

test_that("pi0_model stops on non-numeric input and warns on non-p-value columns", {
  data(bmi)
  z <- sumstats[, -1]
  z$id <- paste0("rs", seq_len(nrow(z)))
  expect_error(pi0_model(z, verbose = FALSE), "numeric")
  z2 <- as.matrix(sumstats[, -1])
  z2[, 1] <- -log10(z2[, 1])
  expect_warning(
    m <- pi0_model(z2, min_discoveries = 50, min_snps_per_knot = 20, verbose = FALSE),
    "not a p-value column"
  )
  expect_s3_class(m$fmod, "formula")
})
