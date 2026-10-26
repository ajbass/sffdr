sim_echisq <- function(n = 20000, seed = 1) {
  set.seed(seed)
  a <- rnorm(n)
  pi1 <- plogis(-4 + a)
  nonnull <- runif(n) < pi1
  z <- rnorm(n) + ifelse(nonnull, rnorm(n, 0, 3), 0)
  list(p = 2 * pnorm(-abs(z)), echisq = 1 + pi1 * 9 + rnorm(n, 0, 0.5), pi0 = 1 - pi1)
}

test_that("pi0_model_echisq returns a model with the anchor and the upper-quantile knots", {
  d <- sim_echisq()
  mod <- pi0_model_echisq(d$echisq, verbose = FALSE)
  expect_type(mod, "list")
  expect_named(mod, c("fmod", "zt", "knots", "anchor"))
  expect_s3_class(mod$zt, "data.frame")
  expect_equal(names(mod$zt), "echisq")
  expect_true(all(mod$zt$echisq > 0 & mod$zt$echisq <= 1))
  # the anchor is the share of SNPs at or below the null value, and it is a knot
  expect_equal(unname(mod$anchor), mean(d$echisq <= 1))
  expect_true(mod$anchor %in% mod$knots$echisq)
  # candidates at or below the anchor are dropped, the rest kept in order
  expect_true(all(diff(mod$knots$echisq) > 0))
  expect_true(all(mod$knots$echisq[-1] > mod$anchor + 0.01))
  expect_match(deparse1(mod$fmod), "splines::ns\\(echisq")
})

test_that("no anchor when the null value is outside the covariate's range", {
  d <- sim_echisq()
  expect_message(mod <- pi0_model_echisq(d$echisq + 5), "no anchor")
  expect_true(is.na(mod$anchor))
  # 20,000 SNPs: the 0.999 bin (20 SNPs) is pruned, the other quantile knots stay
  expect_equal(mod$knots$echisq, c(0.5, 0.75, 0.9, 0.95, 0.99))
  mod2 <- pi0_model_echisq(d$echisq, null_value = -Inf, verbose = FALSE)
  expect_true(is.na(mod2$anchor))
})

test_that("sparse bins lose knots and the anchor survives while a neighbour can go", {
  d <- sim_echisq(n = 2000)
  mod <- pi0_model_echisq(d$echisq, min_snps_per_knot = 100, verbose = FALSE)
  # 2,000 SNPs: the 0.99 and 0.999 bins (20 and 2 SNPs) cannot stand
  expect_false(any(mod$knots$echisq >= 0.99))
  expect_true(mod$anchor %in% mod$knots$echisq)
  # with a huge requirement every knot goes and the term is linear
  mod_lin <- pi0_model_echisq(d$echisq, min_snps_per_knot = 5000, verbose = FALSE)
  expect_length(mod_lin$knots$echisq, 0)
  expect_equal(deparse1(mod_lin$fmod), "~echisq")
})

test_that("indep_snps and weights define the support of a bin", {
  d <- sim_echisq()
  n <- length(d$p)
  few <- rep(FALSE, n)
  few[sample(n, 1500)] <- TRUE
  mod_i <- pi0_model_echisq(d$echisq, indep_snps = few, verbose = FALSE)
  expect_false(any(mod_i$knots$echisq >= 0.99))
  mod_w <- pi0_model_echisq(d$echisq, weights = rep(0.05, n), verbose = FALSE)
  expect_false(any(mod_w$knots$echisq >= 0.99))
  expect_error(pi0_model_echisq(d$echisq, indep_snps = few, weights = rep(0.5, n)), "Ambiguous")
  expect_error(pi0_model_echisq(d$echisq, weights = rep(0.5, 10)), "one entry per row")
})

test_that("several columns give an additive model", {
  d <- sim_echisq()
  x <- cbind(gfm = d$echisq, base = d$echisq + rnorm(length(d$p), 0, 0.3))
  mod <- pi0_model_echisq(x, verbose = FALSE)
  expect_equal(names(mod$zt), c("gfm", "base"))
  expect_length(mod$knots, 2)
  expect_match(deparse1(mod$fmod), "ns\\(gfm.*\\+.*ns\\(base")
})

test_that("the model feeds fpi0est and the estimate falls with the expected chi-squared", {
  d <- sim_echisq()
  mod <- pi0_model_echisq(d$echisq, verbose = FALSE)
  fpi0 <- fpi0est(d$p, pi0_model_obj = mod, verbose = FALSE)$fpi0
  expect_length(fpi0, length(d$p))
  expect_true(all(fpi0 >= 0 & fpi0 <= 1))
  expect_lt(cor(fpi0, d$echisq), 0)
  top <- d$echisq > quantile(d$echisq, 0.99)
  expect_lt(mean(fpi0[top]), mean(fpi0[!top]))
})

test_that("input checks", {
  expect_error(pi0_model_echisq(letters), "numeric")
  expect_error(pi0_model_echisq(c(1, 2, 3), knots = c(0, 0.5)), "quantiles in \\(0, 1\\)")
  expect_error(pi0_model_echisq(c(1, 2, 3), null_value = c(1, 2)), "single number")
  expect_error(pi0_model_echisq(c(1, 2, 3), min_snps_per_knot = 0), "at least 1")
})

test_that("annotation columns alongside the expected chi-squared get their own spline, unanchored", {
  d <- sim_echisq()
  ld <- exp(rnorm(length(d$p), 3, 0.6))                      # an LD score: null value 1 below its range
  expect_message(mod <- pi0_model_echisq(cbind(echisq = d$echisq, ld = ld)), "no anchor")
  expect_equal(names(mod$zt), c("echisq", "ld"))
  expect_false(is.na(mod$anchor[["echisq"]]))
  expect_true(is.na(mod$anchor[["ld"]]))
  expect_true(all(mod$knots$ld %in% c(0.5, 0.75, 0.9, 0.95, 0.99, 0.999)))
  expect_match(deparse1(mod$fmod), "ns\\(echisq.*ns\\(ld")
  fpi0 <- fpi0est(d$p, pi0_model_obj = mod, verbose = FALSE)$fpi0
  expect_true(all(fpi0 > 0 & fpi0 <= 1))
})
