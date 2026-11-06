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
  c(fpi0est(d$p, mp, indep_snps = d$indep, verbose = FALSE)$fpi0) }                 # c(): fpi0 the surrogate, not the covariate

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
  lex <- integer(length(v)); lex[order(fit$flfdr, d$p)] <- seq_along(v)              # the order of the regions: lfdr, then p
  expect_gt(cor(rank(v), lex), 0.999)
  expect_error(sffdr:::vvalues(sffdr(d$p, rep(0.9, length(d$p)), verbose = FALSE)), "keep_fit")
  # null-heavy scan: the lfdr reaches 1 early in the tail; within that tie v is ordered by p and at least p, the steep
  # tail is still the p-value
  set.seed(3); m <- 50000; nonnull <- runif(m) < 0.01
  p <- 2 * pnorm(-abs(ifelse(nonnull, rnorm(m, 0, 4), 0) + rnorm(m)))
  fit <- sffdr(p, rep(0.99, m), keep_fit = TRUE, verbose = FALSE); v <- sffdr:::vvalues(fit)
  sel <- p <= 1e-2 & fit$flfdr >= 1
  expect_gt(sum(sel), 20); expect_true(all(v[sel] >= p[sel] * (1 - 1e-12)))
  expect_true(all(diff(v[sel][order(p[sel])]) >= -1e-15))
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
  cap <- pl$lfdr >= 1                                                                # the cap is exact in the grid pass
  expect_gt(sum(cap), 100); expect_equal(v_two[cap], v_exact[cap], tolerance = 1e-12)
  expect_true(all(v_two >= 0 & v_two <= 1))
})

test_that("a fold refits the tail only, on the same rows and window: the functional p-values of a full refit", {
  d <- sim(); f0 <- fit_pi0(d); ns <- asNamespace("sffdr"); e <- min(d$p); tr <- d$chr != 1
  ft <- ns$.fit_surface(d$p[tr], f0[tr], tail = .1); ff <- ns$.fit_surface(d$p[tr], f0[tr])
  expect_identical(ft[c("z", "groups", "n_bins")], ff[c("z", "groups", "n_bins")]); expect_identical(ff$scale, 1)
  expect_gt(ft$scale, .4); expect_lt(ft$scale, .6)                                  # the weight below p = .5
  skip_on_cran()
  cv <- sffdr(d$p, f0, folds = d$chr, verbose = FALSE); v_full <- rep(NA_real_, length(d$p))
  for (f in 1:5) {                                                                   # leave-one-fold-out full refits (1.2.3)
    tr <- d$chr != f; te <- d$chr == f
    withr::with_seed(2026, { fit <- ns$.fit_surface(d$p[tr], f0[tr], seed = 2026); cur <- ns$.bin_curves(fit, p_min = e)
      pl <- ns$.place_on_fit(fit, cur, d$p[te], f0[te]) })
    rk <- ns$.continue_curves(fit, cur)
    v_full[te] <- ns$.v_from_curves(rk, ns$.on_curves(rk, pl$grp, d$p[te], f0[te]), d$p[te], pl$grp, credit = min(d$p[tr]) > e)
  }
  top <- pmin(cv$fpvalues, v_full) <= 1e-3
  expect_gt(sum(top), 100); expect_gt(cor(log(cv$fpvalues[top]), log(v_full[top])), .99)
  expect_lt(median(abs(log10(cv$fpvalues[top] / v_full[top]))), log10(1.1))
  for (t in c(1e-4, 1e-3)) expect_lte(abs(sum(cv$fpvalues <= t) - sum(v_full <= t)), .1 * sum(v_full <= t) + 5)
})

test_that("pi0 in folds: each fold is scored on its own pi0 model's ranks and pi0, the bulk on the global model's", {
  d <- sim(); ns <- asNamespace("sffdr"); mp <- pi0_model(z = data.frame(u = d$u), verbose = FALSE)
  f0 <- fpi0est(d$p, mp, folds = d$chr, verbose = FALSE)$fpi0; pf <- attr(f0, "pi0_folds")
  attr(f0, "surrogate") <- NULL                                                          # the surrogate of several covariates
  expect_identical(attr(pf$terms, ".Environment"), ns)
  expect_equal(ns$.fpi0_of_fold(pf, NULL), c(fpi0est(d$p, mp, verbose = FALSE)$fpi0))
  cv <- sffdr(d$p, f0, folds = d$chr, verbose = FALSE)
  expect_null(attributes(cv$fpi0)); expect_equal(cv$flfdr, pmin(c(f0) / cv$fx, 1))      # each SNP's lfdr on its held-out pi0
  tr <- d$chr != 2; F2 <- ns$.fpi0_of_fold(pf, 2); e <- min(d$p)                          # fold 2 by hand
  st <- ns$.cross_setup(d$p, c(f0), ns$.fpi0_of_fold(pf, NULL), p_min = e)
  k <- ns$.cross_fold(d$p, F2, F2, tr, st, p_min = e, own_s = TRUE)
  expect_identical(cv$fpvalues[!tr], ns$.v_from_curves(k$rk, c = k$c, p = d$p[!tr], grp = k$pl$grp, credit = min(d$p[tr]) > e))
  cu <- sffdr(d$p, f0, surrogate = d$u, folds = d$chr, verbose = FALSE)                 # a given surrogate: the fold pi0 only
  expect_equal(cu$flfdr, pmin(c(f0) / cu$fx, 1))
  expect_lte(sum(cv$fpvalues[!d$nonnull] <= 1e-3), 1.5 * 1e-3 * sum(!d$nonnull) + 5)
  expect_error(sffdr(d$p, f0, folds = d$chr, nfolds = 2, verbose = FALSE), "other folds")
  expect_identical(sffdr(d$p, c(f0), verbose = FALSE)$fpvalues, sffdr(d$p, f0, verbose = FALSE)$fpvalues)   # no folds: the values
})

test_that("a pi0 model with one covariate: the covariate is the default surrogate, SNPs without it at the null end", {
  d <- sim(); mp <- pi0_model(z = data.frame(u = d$u), indep_snps = d$indep, verbose = FALSE)
  f0 <- fpi0est(d$p, mp, indep_snps = d$indep, verbose = FALSE)$fpi0
  expect_identical(attr(f0, "surrogate"), mp$zt$u)
  expect_identical(sffdr(d$p, f0, verbose = FALSE)$fpvalues, sffdr(d$p, c(f0), surrogate = mp$zt$u, verbose = FALSE)$fpvalues)
  ff <- fpi0est(d$p, mp, folds = d$chr, verbose = FALSE)$fpi0; fs <- ff; attr(fs, "surrogate") <- NULL
  expect_identical(sffdr(d$p, ff, folds = d$chr, verbose = FALSE)$fpvalues, sffdr(d$p, fs, surrogate = mp$zt$u, folds = d$chr, verbose = FALSE)$fpvalues)
  zt <- mp$zt; zt$u[1:50] <- NA
  fm <- fpi0est(d$p, z = zt, pi0_model = mp$fmod, indep_snps = d$indep, verbose = FALSE)$fpi0
  expect_equal(c(fm)[1:50], rep(1, 50)); expect_identical(attr(fm, "surrogate")[1:50], rep(Inf, 50))
  expect_false(anyNA(sffdr(d$p, fm, verbose = FALSE)$fpvalues))
  two <- fpi0est(d$p, z = data.frame(u = mp$zt$u, w = runif(length(d$p))), pi0_model = "~ u + w", verbose = FALSE)$fpi0
  expect_null(attr(two, "surrogate"))
})

test_that("a tied surrogate (few levels) is spread over the training bins: null genomes stay at or below nominal", {
  set.seed(9); m <- 60000; p0 <- runif(m); lev <- sample(c(.95, .97, .99), m, replace = TRUE); chr <- rep(1:6, length.out = m)
  a <- sffdr(p0, lev, folds = chr, verbose = FALSE)
  expect_lte(sum(a$fpvalues <= 1e-3), 1.5 * 60); expect_lte(sum(a$fpvalues <= 1e-4), 1.5 * 6 + 3)
  b <- sffdr(p0, lev, folds = chr, verbose = FALSE); expect_identical(a$fpvalues, b$fpvalues)        # reproducible placement
  cst <- sffdr(p0, rep(.98, m), folds = chr, verbose = FALSE); expect_lte(sum(cst$fpvalues <= 1e-3), 1.5 * 60)
})

test_that("nfolds groups whole labels into folds of similar size and equals passing the grouped labels", {
  d <- sim(); f0 <- fit_pi0(d)
  g <- sffdr:::.group_folds(d$chr, 2)
  expect_equal(sort(unique(g)), 1:2)
  expect_true(all(tapply(g, d$chr, function(x) length(unique(x))) == 1))       # a chromosome is never split
  expect_lte(abs(diff(table(g))), max(table(d$chr)))                           # balanced to within one label
  a <- sffdr(d$p, f0, folds = d$chr, nfolds = 2, verbose = FALSE)
  b <- sffdr(d$p, f0, folds = g, verbose = FALSE)
  expect_identical(a$fpvalues, b$fpvalues); expect_identical(a$folds, g)
  expect_error(sffdr(d$p, f0, folds = d$chr, nfolds = 1, verbose = FALSE), "nfolds")
  expect_error(sffdr(d$p, f0, folds = d$chr, nfolds = 6, verbose = FALSE), "nfolds")
})

test_that("a bin that has not reached the lfdr by the smallest p is credited the range below it", {
  for (s in c(18, 19)) {                                   # weak scans, uninformative surrogate: the top SNP is not lifted
    set.seed(s); m <- 20000; u <- runif(m)
    z <- ifelse(runif(m) < .005, rnorm(m, 0, 1.5), 0) + rnorm(m); p <- 2 * pnorm(-abs(z))
    a <- sffdr(p, rep(.995, m), surrogate = u, folds = rep(1:5, length.out = m), verbose = FALSE)
    i <- which.min(p)
    expect_gte(a$fpvalues[i], p[i] * (1 - 1e-12))
  }
})

test_that("only the fold whose own SNPs hold the smallest p is credited the range below it", {
  s_p <- qnorm(c(1e-6, 1e-4, 1e-2, .1, .5)); L <- c(.2, .5, .5, .8, 1)
  expect_equal(sffdr:::.t_of_c(L, s_p, .1), 1e-6); expect_equal(sffdr:::.t_of_c(L, s_p, .1, credit = FALSE), 0)
  expect_equal(sffdr:::.t_of_c(L, s_p, .35, credit = FALSE), sffdr:::.t_of_c(L, s_p, .35))
  skip_on_cran()
  ns <- asNamespace("sffdr")
  one <- function(p, f0, folds, f, credit) {                                 # fold f as .sffdr_cross computes it
    st <- ns$.cross_setup(p, f0, f0, p_min = min(p)); k <- ns$.cross_fold(p, f0, f0, folds != f, st, p_min = min(p))
    ns$.v_from_curves(k$rk, c = k$c, p = p[folds == f], grp = k$pl$grp, credit = credit)
  }
  set.seed(19); m <- 20000; u <- runif(m); pi1 <- .005 * 8 * exp(-8 * u) / (1 - exp(-8))
  z <- ifelse(runif(m) < pi1, rnorm(m, 0, 1.5), 0) + rnorm(m); p <- 2 * pnorm(-abs(z)); f0 <- 1 - pi1
  folds <- rep(1:5, length.out = m); fm <- folds[which.min(p)]; fo <- setdiff(1:5, fm)[1]
  a <- sffdr(p, f0, surrogate = f0, folds = folds, verbose = FALSE)
  expect_identical(a$fpvalues[folds == fm], one(p, f0, folds, fm, TRUE))     # holds the smallest p
  nc <- one(p, f0, folds, fo, FALSE)
  expect_identical(a$fpvalues[folds == fo], nc)                              # its training holds it
  expect_true(any(nc < one(p, f0, folds, fo, TRUE)))
  p[which(folds == fo)[1]] <- min(p)                                         # tied across two folds: neither credited
  b <- sffdr(p, f0, surrogate = f0, folds = folds, verbose = FALSE)
  expect_identical(b$fpvalues[folds == fm], one(p, f0, folds, fm, FALSE))
})

test_that("a tie in the lfdr is broken by p inside the region: no atom at 1, uniform under the null, copies share values", {
  s_p <- qnorm(c(1e-6, 1e-4, 1e-2, .1, .5)); L1 <- c(.2, .5, .5, .8, 1); L2 <- c(.1, .3, .6, .9, 1)
  expect_equal(sffdr:::.a_of_c(L1, s_p, .5), 1e-4); expect_equal(sffdr:::.t_of_c(L1, s_p, .5), 1e-2)  # the stretch at .5
  expect_equal(sffdr:::.a_of_c(L1, s_p, .35), sffdr:::.t_of_c(L1, s_p, .35))                       # a crossing: one point
  expect_equal(sffdr:::.a_of_c(L1, s_p, .1), 1e-6)
  cur <- list(s_p = s_p, lfdr = cbind(L1, L2), W = c(1, 1), ids = 1:2)
  a2 <- sffdr:::.a_of_c(L2, s_p, .5)
  expect_equal(sffdr:::.v_from_curves(cur, c(.5, .5), c(1e-3, 5e-2), c(2, 1)), c((1e-3 + max(a2, 1e-3)) / 2, (5e-2 + a2) / 2))
  set.seed(4); m <- 60000; nb <- m / 5                         # null genome in blocks of 5 exact copies (same p and fpi0)
  p <- rep(runif(nb), each = 5); f0 <- rep(1 - .02 * runif(nb), each = 5); chr <- rep(rep(1:6, length.out = nb), each = 5)
  a <- sffdr(p, f0, folds = chr, verbose = FALSE)
  cap <- a$flfdr >= 1
  expect_gt(mean(cap), .5); expect_lt(mean(a$fpvalues > .999), .01)
  expect_true(all(a$fpvalues[cap] >= pmin(p[cap], 1 - 1e-4) * (1 - 1e-12)))      # at least p below the grid top
  for (t in c(.05, .2, .5)) expect_lt(abs(mean(a$fpvalues <= t) / t - 1), .1)
  expect_true(all(tapply(a$fpvalues, rep(seq_len(nb), each = 5), function(x) diff(range(x)) == 0)))
})

test_that("below a fold's smallest training p the ranking falls with p: a lone locus is called, nothing else moves", {
  skip_on_cran()
  set.seed(1); m <- 20000; p0 <- 2 * pnorm(-abs(ifelse(runif(m) < .01, rnorm(m, 0, 3), 0) + rnorm(m))); f00 <- runif(m, .97, 1)
  fit <- sffdr:::.fit_surface(p0, f00); e <- min(p0)
  cur <- sffdr:::.bin_curves(fit, p_min = e / 1e6); rk <- sffdr:::.continue_curves(fit, cur); lo <- cur$s_p < qnorm(e)
  expect_identical(rk$fx[!lo, ], cur$fx[!lo, ]); expect_identical(rk$lfdr[!lo, ], cur$lfdr[!lo, ])
  a <- sweep(rk$fx[lo, , drop = FALSE], 1, e / pnorm(cur$s_p[lo]), "/")             # below e: a value at e times e / p
  expect_lt(max(abs(sweep(a, 2, a[1, ], "/") - 1)), 1e-12); expect_true(all(a[1, ] >= cur$fx[max(which(lo)) + 1L, ] * (1 - 1e-12)))
  expect_true(all(apply(rk$fx, 2, function(x) all(diff(x) <= 0))))
  expect_identical(sffdr:::.continue_curves(fit, sffdr:::.bin_curves(fit, p_min = e)), sffdr:::.bin_curves(fit, p_min = e))
  ns <- asNamespace("sffdr")
  set.seed(2); n_chr <- 8; m <- 10000; n <- n_chr * m; chr <- rep(seq_len(n_chr), each = m)
  z <- rnorm(n); enr <- sample(n, n * .02); bg <- sample(enr, 120L); z[bg] <- rnorm(120L, 0, 2)
  hit <- setdiff(which(chr == 1), enr)[1:20]; z[hit] <- 7.9 + rnorm(20, 0, .15)
  p <- 2 * pnorm(-abs(z)); f0 <- ifelse(seq_len(n) %in% enr, .8, .995) * runif(n, .995, 1)
  a <- sffdr(p, f0, folds = chr, verbose = FALSE)
  st <- ns$.cross_setup(p, f0, f0, p_min = min(p))
  fit_121 <- function(f) {                                                   # fold f ranked by the fit's curves
    tr <- chr != f; te <- chr == f; k <- ns$.cross_fold(p, f0, f0, tr, st, p_min = min(p)); cur <- k$cur; pl <- k$pl
    list(v = ns$.v_from_curves(cur, c = pl$lfdr, p = p[te], grp = pl$grp, credit = min(p[tr]) > min(p)), lfdr = pl$lfdr,
         edge = pnorm(cur$s_p[min(which(cur$s_p >= qnorm(min(p[tr]))))]), e = min(p[tr]))
  }
  o1 <- fit_121(1); o2 <- fit_121(2); te <- chr == 1
  expect_identical(a$flfdr[te], o1$lfdr)
  expect_identical(a$fpvalues[chr == 2], o2$v)
  expect_true(all(a$fpvalues[hit] <= 5e-8))
  expect_lt(sum(o1$v[match(hit, which(te))] <= 5e-8), 15)
  inr <- p[te] >= o1$e & o1$v <= 1e-3 & a$fpvalues[te] <= 1e-3
  expect_gt(sum(inr), 10); expect_lte(max(abs(a$fpvalues[te][inr] - o1$v[inr])), o1$edge * (1 + 1e-9))
})

test_that("where the continued ranking passes the largest double the read is a number, and only the NaN reads move", {
  s <- seq(qnorm(.Machine$double.xmin), qnorm(1 - 1e-4), length.out = 600)       # the grid of .bin_curves at the floor of p
  fx <- outer(1 / pnorm(s), c(10, 1))                                             # ratio e / p as continued: bin 1 passes the largest double at the floor
  cur <- list(s_p = s, fx = fx, lfdr = sweep(fx, 2, c(.9, .95), function(f, p0) pmin(p0 / f, 1)), W = c(500, 500), ids = 1:2)
  expect_true(is.infinite(fx[1, 1])); expect_true(all(is.finite(fx[-1, ])))
  pt <- c(.Machine$double.xmin, pnorm(mean(s[1:2])), 1e-300, 1e-50, 1e-6, .5); grp <- rep(1L, length(pt)); f1 <- rep(.9, length(pt))
  q0 <- pmin(f1 / pmax(approx(s, fx[, 1], xout = qnorm(pt), rule = 2)$y, .Machine$double.xmin), 1)    # the plain read, as before 1.2.5
  expect_identical(is.nan(q0), c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE))
  q <- sffdr:::.on_curves(cur, grp, pt, f1); v <- sffdr:::.v_from_curves(cur, q, pt, grp)
  expect_identical(q[-2], q0[-2]); expect_true(is.finite(q[2]))
  expect_false(anyNA(v)); expect_true(all(diff(q) >= 0)); expect_true(all(diff(v) >= 0))
})

test_that("the v sums in C++ are the R loops over the bins, to the bit", {
  set.seed(3); G <- 60; B <- 30; s_p <- seq(qnorm(1e-10), qnorm(1 - 1e-4), length.out = G)
  L <- apply(matrix(runif(G * B), G), 2, function(x) pmin(sort(x) * 1.3, 1)); Wn <- runif(B); Wn <- Wn / sum(Wn)
  cg <- sort(unique(c(exp(seq(log(1e-6), 0, length.out = 500)), seq(.5, 1, length.out = 500))))
  c <- c(runif(150), rep(1, 20), L[cbind(sample(G, 30, TRUE), sample(B, 30, TRUE))]); p <- runif(200, 0, .01)
  grp <- sample(c(seq_len(B), 99L), 200, TRUE); own <- match(grp, seq_len(B)) - 1L; own[is.na(own)] <- -1L
  for (cr in c(TRUE, FALSE)) {
    V <- numeric(length(cg)); for (b in seq_len(B)) V <- V + Wn[b] * sffdr:::.t_of_c(L[, b], s_p, cg, cr)
    expect_identical(sffdr:::.v_grid(L, s_p, cg, Wn, cr), V)
    v <- numeric(length(c))
    for (b in seq_len(B)) {
      a <- sffdr:::.a_of_c(L[, b], s_p, c); tb <- pmin(sffdr:::.t_of_c(L[, b], s_p, c, cr), pmax(a, p))
      tb[grp == b] <- pmax(a[grp == b], p[grp == b]); v <- v + Wn[b] * tb
    }
    expect_identical(sffdr:::.v_exact(L, s_p, c, p, own, Wn, cr), v)
  }
})
