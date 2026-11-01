## sffdr(folds = ...): the cross-fitted functional p-value (the v-value construction) and local FDRs. keep_fit leaves
## the default result unchanged; the model-form functional p-value returns the p-value under an uninformative
## surrogate; cross-fitting by fold is complete, reproducible and parallel-safe; the fit-only training path is the
## surface sffdr() keeps and the two-pass evaluation is exact in the tail; the self-fit version (internal) is not it.
sim <- function(m = 20000, k = 10, seed = 3) {
  set.seed(seed); nb <- m %/% k
  u <- rep(runif(nb), each = k); pi0 <- 1 - 0.25 * exp(-8 * u)
  nonnull <- rep(runif(nb), each = k) < 1 - pi0
  mu <- ifelse(nonnull, rep(rnorm(nb, 0, 4), each = k), 0)
  z <- mu + sqrt(.9) * rep(rnorm(nb), each = k) + sqrt(.1) * rnorm(m)
  list(p = 2 * pnorm(-abs(z)), u = u, chr = rep(1:5, length.out = nb)[rep(seq_len(nb), each = k)],
       indep = rep(c(TRUE, rep(FALSE, k - 1)), nb), nonnull = nonnull)
}
fit_pi0 <- function(d) { mp <- pi0_model(z = data.frame(u = d$u), indep_snps = d$indep, verbose = FALSE)
  fpi0est(d$p, mp, indep_snps = d$indep, verbose = FALSE)$fpi0 }

test_that("keep_fit keeps the fit and changes nothing else", {
  d <- sim(); f0 <- fit_pi0(d)
  a <- sffdr(d$p, f0, verbose = FALSE); b <- sffdr(d$p, f0, keep_fit = TRUE, verbose = FALSE)
  expect_null(a$fit)
  expect_identical(a[setdiff(names(a), "call")], b[setdiff(names(a), "call")])
  expect_true(all(c("lfit", "z", "groups", "n_bins", "p", "fpi0", "surrogate", "fx", "lfdr", "valid") %in% names(b$fit)))
})

test_that("the model-form functional p-value is the p-value under an uninformative surrogate", {
  d <- sim()
  fit <- sffdr(d$p, rep(0.9, length(d$p)), keep_fit = TRUE, verbose = FALSE)
  v <- sffdr:::vvalues(fit)
  expect_length(v, length(d$p)); expect_false(anyNA(v))
  top <- d$p <= 1e-3
  expect_lt(median(abs(log10(v[top] / d$p[top]))), log10(1.1))
  expect_lt(max(abs(log10(v[top] / d$p[top]))), log10(2))
  expect_gt(cor(rank(v), rank(fit$flfdr), method = "spearman"), 0.999)
  expect_error(sffdr:::vvalues(sffdr(d$p, rep(0.9, length(d$p)), verbose = FALSE)), "keep_fit")
  # null-heavy scan: the lfdr reaches 1 early in the tail; those SNPs have v = 1, the steep tail is still the p-value
  set.seed(3); m <- 50000; nonnull <- runif(m) < 0.01
  p <- 2 * pnorm(-abs(ifelse(nonnull, rnorm(m, 0, 4), 0) + rnorm(m)))
  fit <- sffdr(p, rep(0.99, m), keep_fit = TRUE, verbose = FALSE); v <- sffdr:::vvalues(fit)
  sel <- p <= 1e-2 & fit$flfdr >= 1
  expect_gt(sum(sel), 20); expect_true(all(v[sel] >= 1 - 1e-3))
  tail <- p <= 1e-4
  expect_lt(median(abs(log10(v[tail] / p[tail]))), log10(1.1))
})

test_that("sffdr(folds = ) scores every SNP from the other folds, reproducibly, with the same elements", {
  d <- sim(); f0 <- fit_pi0(d)
  a <- sffdr(d$p, f0, verbose = FALSE)
  cv <- sffdr(d$p, f0, folds = d$chr, verbose = FALSE)
  expect_s3_class(cv, "sffdr")
  expect_identical(names(cv), c("call", "pvalues", "fpvalues", "fqvalues", "flfdr", "fpi0", "fx", "folds"))
  for (el in c("fpvalues", "fqvalues", "flfdr", "fx")) { expect_length(cv[[el]], length(d$p)); expect_false(anyNA(cv[[el]])) }
  expect_true(all(cv$fpvalues >= 0 & cv$fpvalues <= 1)); expect_true(all(cv$flfdr >= 0 & cv$flfdr <= 1))
  expect_identical(cv$fpi0, f0)
  cv2 <- sffdr(d$p, f0, folds = d$chr, verbose = FALSE); expect_identical(cv2$fpvalues, cv$fpvalues)
  top <- d$p <= 1e-3                                                    # the two estimators agree in the tail
  expect_gt(cor(log(cv$fpvalues[top]), log(a$fpvalues[top])), 0.95)
  expect_gt(cor(cv$flfdr[top], a$flfdr[top]), 0.95)
  expect_error(sffdr(d$p, f0, folds = rep(1, length(d$p)), verbose = FALSE), "two distinct")
  expect_error(sffdr(d$p, f0, folds = 1:3, verbose = FALSE), "one label per SNP")
  pc <- sffdr:::.sffdr_cross(d$p, f0, f0, d$chr, folds_run = 1:2, verbose = FALSE)   # job-array pieces pool to the same result
  pd <- sffdr:::.sffdr_cross(d$p, f0, f0, d$chr, folds_run = 3:5, verbose = FALSE)
  expect_true(anyNA(pc$fpvalues))
  mg <- sffdr:::.merge_folds(list(pc, pd), d$p)
  expect_identical(mg$fpvalues, cv$fpvalues); expect_identical(mg$fqvalues, cv$fqvalues)
  if (.Platform$OS.type != "windows") {
    cp <- sffdr(d$p, f0, folds = d$chr, ncores = 2, verbose = FALSE); expect_identical(cp$fpvalues, cv$fpvalues)
  }
})

test_that("the fit-only surface is the surface sffdr(keep_fit = TRUE) keeps, and the two-pass evaluation is exact in the tail", {
  d <- sim(); f0 <- fit_pi0(d)
  tr <- d$chr != 1; te <- which(!tr)
  a <- sffdr(d$p[tr], f0[tr], keep_fit = TRUE, verbose = FALSE)$fit; b <- sffdr:::.fit_surface(d$p[tr], f0[tr])
  expect_identical(a$z, b$z); expect_identical(a$groups, b$groups)
  ca <- sffdr:::.bin_curves(a, min(d$p), 600L); cb <- sffdr:::.bin_curves(b, min(d$p), 600L)
  expect_identical(ca$lfdr, cb$lfdr); expect_identical(ca$W, cb$W)
  pl <- sffdr:::.place_on_fit(b, cb, d$p[te], f0[te])
  v_exact <- sffdr:::.v_from_curves(cb, pl$lfdr, d$p[te], pl$grp, exact_below = Inf)
  v_two <- sffdr:::.v_from_curves(cb, pl$lfdr, d$p[te], pl$grp, n_grid = 100L)      # forces the grid pass at this size
  sel <- v_exact <= 1e-2
  expect_gt(sum(sel), 10); expect_equal(v_two[sel], v_exact[sel], tolerance = 1e-12)
  expect_true(all(v_two >= 0 & v_two <= 1))
})

test_that("a tied surrogate (few levels) is spread over the training bins: null genomes stay at or below nominal", {
  set.seed(9); m <- 60000; p0 <- runif(m); lev <- sample(c(.95, .97, .99), m, replace = TRUE); chr <- rep(1:6, length.out = m)
  a <- sffdr(p0, lev, folds = chr, verbose = FALSE)
  expect_lte(sum(a$fpvalues <= 1e-3), 1.5 * 60); expect_lte(sum(a$fpvalues <= 1e-4), 1.5 * 6 + 3)
  b <- sffdr(p0, lev, folds = chr, verbose = FALSE); expect_identical(a$fpvalues, b$fpvalues)        # reproducible placement
  cst <- sffdr(p0, rep(.98, m), folds = chr, verbose = FALSE); expect_lte(sum(cst$fpvalues <= 1e-3), 1.5 * 60)
})
