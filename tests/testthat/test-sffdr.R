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
