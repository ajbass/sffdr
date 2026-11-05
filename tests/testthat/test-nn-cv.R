## sffdr(folds = , nn = "cv"): the window is one of the candidates and the result is the cross-fit at that window
test_that("nn = \"cv\" picks a candidate window and returns the cross-fit at it", {
  expect_error(sffdr(runif(100), rep(.9, 100), nn = "cv", verbose = FALSE), "folds")
  skip_on_cran()
  set.seed(5); m <- 20000; u <- runif(m)
  z <- ifelse(runif(m) < .1 * exp(-3 * u), rnorm(m, 0, 3), 0) + rnorm(m); p <- 2 * pnorm(-abs(z))
  f0 <- 1 - .1 * exp(-3 * u); chr <- rep(1:5, length.out = m)
  a <- sffdr(p, f0, folds = chr, nn = "cv", verbose = FALSE); nn <- attr(a, "nn")
  n_tail <- sum(p < pnorm(-2)); nn0 <- max(.005, min(.05, 5000 / (n_tail + min(m - n_tail, 1e5))))
  expect_true(any(abs(nn / nn0 - c(.25, .5, 1, 2, 4)) < 1e-12))
  expect_identical(a$fpvalues, sffdr(p, f0, folds = chr, nn = nn, verbose = FALSE)$fpvalues)
  expect_null(attr(sffdr(p, f0, folds = chr, verbose = FALSE), "nn"))
})
