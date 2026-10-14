# =================================================================
# Tests for the overlap correction: decorrelate_informative() and estimate_overlap()
# =================================================================

library(testthat)

# -----------------------------------------------------------------
# Simulation helper: a primary study and K informative studies on
# shared samples (null noise correlation Sigma, primary first).
# Primary non-null at rate pi1; informative non-null at rate pik,
# enriched among primary non-nulls (shared) by the factor `share`.
# -----------------------------------------------------------------
sim_overlap <- function(
  n = 20000L,
  Sigma = matrix(c(1, 0.4, 0.4, 1), 2),
  pi1 = 0.02,
  pik = 0.05,
  v1 = 4,
  vk = 16,
  share = 0.5,
  seed = 1L
) {
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  set.seed(seed)
  K <- ncol(Sigma) - 1
  c1 <- stats::runif(n) < pi1
  ck <- matrix(FALSE, n, K)
  for (k in seq_len(K)) {
    ck[, k] <- ifelse(c1, stats::runif(n) < share, stats::runif(n) < pik)
  }
  mu1 <- ifelse(c1, stats::rnorm(n, 0, sqrt(v1)), 0)
  muk <- ifelse(ck, stats::rnorm(n * K, 0, sqrt(vk)), 0)
  e <- matrix(stats::rnorm(n * (K + 1)), n) %*% chol(Sigma)
  list(z1 = mu1 + e[, 1], zi = muk + e[, -1, drop = FALSE], c1 = c1)
}

tail_ratio <- function(p, t) sum(p < t, na.rm = TRUE) / (t * sum(!is.na(p)))

# =================================================================
# 1. estimate_overlap
# =================================================================

test_that("estimate_overlap: recovers rho with single-study signal, negative rho, K = 3", {
  dat <- sim_overlap(n = 20000L, Sigma = matrix(c(1, 0.3, 0.3, 1), 2), pi1 = 0.02, pik = 0.05, share = 0.1)
  expect_lt(abs(estimate_overlap(dat$z1, dat$zi)[1, 2] - 0.3), 0.03)
  dat <- sim_overlap(n = 20000L, Sigma = matrix(c(1, -0.6, -0.6, 1), 2), pi1 = 0.02, pik = 0.05, share = 0.1, seed = 4L)
  expect_lt(abs(estimate_overlap(dat$z1, dat$zi)[1, 2] + 0.6), 0.03)
  S <- matrix(c(1, .7, .1, .25, .7, 1, .15, .15, .1, .15, 1, .35, .25, .15, .35, 1), 4)
  dat3 <- sim_overlap(n = 20000L, Sigma = S, pi1 = 0.02, pik = 0.05, seed = 5L)
  R <- estimate_overlap(dat3$z1, dat3$zi)
  expect_equal(dim(R), c(4L, 4L))
  expect_true(min(eigen(R, symmetric = TRUE)$values) > 0)
  expect_lt(max(abs(unname(R) - S)), 0.04)
})

test_that("estimate_overlap: dense signal biases toward zero by a few hundredths", {
  dat <- sim_overlap(n = 20000L, Sigma = matrix(c(1, -0.6, -0.6, 1), 2), pi1 = 0.05, pik = 0.1, share = 0.1, seed = 4L)
  r <- estimate_overlap(dat$z1, dat$zi)[1, 2]
  expect_lt(r, -0.55)
  expect_gt(r, -0.62)
})

test_that("estimate_overlap: near-zero rho and the box moment", {
  dat <- sim_overlap(n = 20000L, Sigma = diag(2), seed = 6L)
  expect_lt(abs(estimate_overlap(dat$z1, dat$zi)[1, 2]), 0.03)
  m <- sffdr:::.box_moment
  expect_equal(m(0), 0)
  expect_equal(m(0.5), -m(-0.5), tolerance = 1e-12)
  expect_true(all(diff(sapply(seq(-0.9, 0.9, by = 0.1), m)) > 0))
})

# =================================================================
# 2. The replaced overlap functions
# =================================================================

test_that("the replaced overlap functions are defunct and point to the new ones", {
  expect_error(overlap_null_density(), "decorrelate_informative")
  expect_error(run_empirical_overlap_correction(1, 2, 3, 4), "decorrelate_informative")
  expect_error(estimate_rho_overlap(1, 2), "estimate_overlap")
  expect_error(discover_empirical_nulls(1, 2), "estimate_overlap")
})

# =================================================================
# 3. decorrelate_informative
# =================================================================

test_that("decorrelate_informative: the formula, vector and matrix input", {
  dat <- sim_overlap(n = 5000L, Sigma = matrix(c(1, .5, .3, .5, 1, .2, .3, .2, 1), 3))
  d1 <- decorrelate_informative(dat$z1, dat$zi[, 1], rho_o = 0.5)
  expect_named(d1, c("z2", "p2", "rho_o"))
  expect_equal(d1$z2, (dat$zi[, 1] - 0.5 * dat$z1) / sqrt(1 - 0.25))
  expect_equal(d1$p2, 2 * pnorm(-abs(d1$z2)))
  d2 <- decorrelate_informative(dat$z1, dat$zi, rho_o = c(0.5, 0.3))
  expect_equal(dim(d2$z2), c(5000L, 2L))
  expect_equal(d2$z2[, 1], d1$z2)
  # sffdr 1.1.2 calls (p-values passed positionally) still work; lfdr_x is ignored
  p <- 2 * pnorm(-abs(cbind(dat$z1, dat$zi)))
  expect_equal(decorrelate_informative(dat$z1, dat$zi[, 1], p[, 1], p[, 2], rho_o = 0.5)$z2, d1$z2)
  expect_warning(decorrelate_informative(dat$z1, dat$zi[, 1], rho_o = 0.5, lfdr_x = rep(1, 5000)), "no longer used")
  expect_error(decorrelate_informative(dat$z1, dat$zi, rho_o = 0.5), "one correlation")
  expect_error(decorrelate_informative(dat$z1, dat$zi[, 1], rho_o = 1), "one correlation")
  expect_error(decorrelate_informative(abs(dat$z1), dat$zi[, 1], rho_o = 0.5), "signed")
  expect_lt(abs(decorrelate_informative(dat$z1, dat$zi[, 1])$rho_o - 0.5), 0.05)
})

test_that("decorrelate_informative: the covariate is independent of the primary's null p-values", {
  # primary null; informative effects at 10%, all shifted the same way
  dat <- sim_overlap(n = 50000L, Sigma = matrix(c(1, .7, .7, 1), 2), pi1 = 0, pik = .1, vk = 16, seed = 11L)
  p1 <- 2 * pnorm(-abs(dat$z1))
  raw <- 2 * pnorm(-abs(dat$zi[, 1]))
  dec <- decorrelate_informative(dat$z1, dat$zi[, 1], rho_o = 0.7)
  expect_gt(stats::cor(-log(p1), -log(raw), method = "spearman"), 0.2)
  expect_lt(abs(stats::cor(-log(p1), -log(dec$p2), method = "spearman")), 0.01)
  strong <- dec$p2 < 1e-3
  expect_equal(tail_ratio(p1[strong], 0.05), 1, tolerance = 0.25)
})

test_that("sffdr on decorrelated covariates: no false discoveries when the primary is null", {
  skip_on_cran()
  dat <- sim_overlap(n = 20000L, Sigma = matrix(c(1, .9, .9, 1), 2), pi1 = 0, pik = .1, vk = 25, seed = 12L)
  p1 <- 2 * pnorm(-abs(dat$z1))
  dec <- decorrelate_informative(dat$z1, dat$zi[, 1], rho_o = 0.9)
  m <- pi0_model(as.matrix(dec$p2), verbose = FALSE)
  f <- fpi0est(p1, pi0_model_obj = m, verbose = FALSE)
  res <- sffdr(p1, fpi0 = f$fpi0, verbose = FALSE)
  expect_equal(sum(res$fqvalues <= 0.05), 0)
})
