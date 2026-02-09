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
