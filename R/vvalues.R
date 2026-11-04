## ============================================================================
## vvalues.R -- cross-fitting behind sffdr(folds = ...): the density is fitted
## leave-one-fold-out and each SNP's functional p-value is computed under the
## fit that did not see its fold (the v-value of Liley and Wallace 2021):
## v_i = sum_b W_b t_b(lfdr_i) / sum_b W_b over the surrogate bins b, with W_b
## the null mass of the bin and t_b(c) the largest p at which its lfdr curve
## reaches c. A tie in the lfdr (the cap at 1, a flat stretch of a curve) is
## broken by p inside the region. The self-fit version (vvalues()) is internal,
## for tests.
## ============================================================================

## Per fold: fit the surface on the other folds, take its bin curves, place the
## held-out SNPs on it, read their functional p-values off the curves.
## folds_run restricts the run to some folds (job arrays); .merge_folds pools.
#' @importFrom stats approx predict qnorm pnorm dnorm
#' @noRd
.sffdr_cross <- function(p, fpi0, surrogate, folds, ncores = 1L, folds_run = NULL, grid_size = 600L, seed = 2026, verbose = TRUE, ...) {
  n <- length(p)
  fold_ids <- sort(unique(folds[!is.na(folds)]))
  if (length(fold_ids) < 2L) stop("'folds' must have at least two distinct labels.")
  run_ids <- if (is.null(folds_run)) fold_ids else intersect(fold_ids, folds_run)
  if (!length(run_ids)) stop("'folds_run' matches no fold label.")
  ok <- !is.na(p) & !is.na(fpi0) & !is.na(surrogate) & !is.na(folds)
  p_min <- min(p[ok])
  one_fold <- function(f) {
    if (!is.null(seed)) withr::local_seed(seed)     # the fit and the tied-surrogate placement reproducible per fold
    te <- ok & folds == f; tr <- ok & folds != f
    if (verbose) message(sprintf("  Fold %s: %d SNPs held out, fitting on %d", as.character(f), sum(te), sum(tr)))
    fit <- .fit_surface(p[tr], fpi0[tr], surrogate[tr], seed = seed, ...)
    cur <- .bin_curves(fit, p_min = p_min, grid_size = grid_size)
    pl <- .place_on_fit(fit, cur, p_new = p[te], fpi0_new = fpi0[te], surrogate_new = surrogate[te])
    v <- .v_from_curves(cur, c = pl$lfdr, p = p[te], grp = pl$grp, credit = min(p[tr]) > p_min)
    list(idx = which(te), lfdr = pl$lfdr, v = v, fx = pl$fx)
  }
  if (verbose) message(sprintf("  Cross-fitting the density: %d folds%s", length(fold_ids), if (ncores > 1L) sprintf(", %d cores", as.integer(ncores)) else ""))
  res <- if (ncores > 1L && .Platform$OS.type != "windows") parallel::mclapply(run_ids, one_fold, mc.cores = as.integer(ncores)) else lapply(run_ids, one_fold)
  bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r$idx), logical(1))
  if (any(bad)) stop(sprintf("fold(s) %s failed", paste(run_ids[bad], collapse = ", ")))
  out <- list(flfdr = rep(NA_real_, n), fpvalues = rep(NA_real_, n), fx = rep(NA_real_, n))
  for (r in res) { out$flfdr[r$idx] <- r$lfdr; out$fpvalues[r$idx] <- r$v; out$fx[r$idx] <- r$fx }
  out
}

## Group fold labels into nfolds folds of similar size, whole labels kept
## together (largest label first, each to the smallest fold so far).
.group_folds <- function(folds, nfolds) {
  nfolds <- as.integer(nfolds)
  tab <- sort(table(folds[!is.na(folds)]), decreasing = TRUE)
  if (nfolds < 2L || nfolds > length(tab)) stop("'nfolds' must be between 2 and the number of distinct labels.")
  size <- numeric(nfolds); grp <- integer(length(tab)); names(grp) <- names(tab)
  for (i in seq_along(tab)) { k <- which.min(size); grp[i] <- k; size[k] <- size[k] + tab[i] }
  out <- grp[as.character(folds)]; out[is.na(folds)] <- NA_integer_
  as.integer(unname(out))
}

## Pool pieces of .sffdr_cross run on subsets of folds, then the functional
## q-values from the pooled held-out local FDRs.
.merge_folds <- function(pieces, p, fp_ties = TRUE) {
  out <- pieces[[1]]
  for (pc in pieces[-1]) for (el in c("flfdr", "fpvalues", "fx")) { have <- !is.na(pc[[el]]); out[[el]][have] <- pc[[el]][have] }
  ok <- !is.na(out$flfdr) & !is.na(p)
  out$fqvalues <- rep(NA_real_, length(p))
  out$fqvalues[ok] <- if (fp_ties) fpvalues(out$flfdr[ok], p[ok])$fq else fpvalues(out$flfdr[ok])$fq
  out
}

## Self-fit v-values from a kept fit (sffdr(keep_fit = TRUE)). Internal: the
## surface has seen the SNPs it scores, so this is anti-conservative. Tests only.
vvalues <- function(object, grid_size = 600L) {
  fit <- object$fit
  if (is.null(fit)) stop("'object' has no kept fit: run sffdr() with keep_fit = TRUE.")
  cur <- .bin_curves(fit, p_min = min(fit$p), grid_size = grid_size)
  v_valid <- .v_from_curves(cur, c = fit$lfdr, p = fit$p, grp = fit$groups)
  out <- rep(NA_real_, length(fit$valid)); out[fit$valid] <- v_valid
  out
}

## -- internals ---------------------------------------------------------------

## The bins' density and lfdr curves on a probit-scale grid of p, from the
## kept fit: the surface at each bin's surrogate centre, the running minimum
## along p (as the fit), then L_b(p) = min(pi0_b / f_b(p), 1).
.bin_curves <- function(fit, p_min, grid_size = 600L) {
  s_p <- seq(qnorm(max(p_min, .Machine$double.xmin)), qnorm(1 - 1e-4), length.out = grid_size)
  zb <- tapply(fit$z, fit$groups, mean)
  W <- tapply(fit$fpi0, fit$groups, sum)
  pi0b <- tapply(fit$fpi0, fit$groups, mean)
  ids <- as.integer(names(zb))
  nd <- data.frame(V1 = rep(qnorm(as.numeric(zb)), each = grid_size), V2 = rep(s_p, times = length(zb)))
  fs <- predict(fit$lfit, newdata = nd, log = FALSE)
  fx <- fs / pmax(dnorm(nd$V1) * dnorm(nd$V2), .Machine$double.xmin)
  fx <- matrix(pmax(fx, .Machine$double.xmin), nrow = grid_size)
  fx <- apply(fx, 2, cummin)
  L <- sweep(fx, 2, as.numeric(pi0b), FUN = function(f, p0) pmin(p0 / f, 1))
  rz <- range(fit$z); dz <- diff(rz)
  list(s_p = s_p, fx = fx, lfdr = L, W = as.numeric(W), ids = ids,
       breaks = seq(rz[1] - dz / 1000, rz[2] + dz / 1000, length.out = fit$n_bins + 1))
}

## t_b(c): the largest p with L_b(p) <= c, linear in probit(p) between grid
## points; the grid top when the bin always reaches c. A bin that has not
## reached c by the grid bottom gets the range below it with credit (the fold
## holding the smallest p), else 0
.t_of_c <- function(L, s_p, c, credit = TRUE) {
  G <- length(s_p)
  i <- findInterval(c, L)
  out <- rep(if (credit) pnorm(s_p[1]) else 0, length(c))
  top <- i >= G; out[top] <- pnorm(s_p[G])
  m <- i >= 1 & !top
  if (any(m)) {
    i1 <- i[m]; L1 <- L[i1]; L2 <- L[i1 + 1]
    w <- ifelse(L2 > L1, (c[m] - L1) / (L2 - L1), 0)
    out[m] <- pnorm(s_p[i1] + w * (s_p[i1 + 1] - s_p[i1]))
  }
  out
}

## a_b(c): the smallest p with L_b(p) >= c, where the bin's stretch at level c
## starts (= t_b(c) where the curve crosses c); the grid bottom when the curve
## starts at or above c, the grid top when it never reaches c.
.a_of_c <- function(L, s_p, c) {
  G <- length(s_p)
  i <- findInterval(c, L, left.open = TRUE)
  out <- rep(pnorm(s_p[1]), length(c))
  top <- i >= G; out[top] <- pnorm(s_p[G])
  m <- i >= 1 & !top
  if (any(m)) {
    i1 <- i[m]; L1 <- L[i1]; L2 <- L[i1 + 1]
    out[m] <- pnorm(s_p[i1] + (c[m] - L1) / (L2 - L1) * (s_p[i1 + 1] - s_p[i1]))
  }
  out
}

## sum_b w_b (top - max(a_b, p))_+ for each p: the null mass of the bins' stretches at lfdr 1 above p
.cap_share <- function(a, w, top, p) {
  k <- a < top; a <- a[k]; w <- w[k]
  if (!length(a)) return(numeric(length(p)))
  o <- order(a); a <- a[o]; w <- w[o]
  j <- findInterval(p, a)                                       # bins whose stretch starts at or below p
  w_lo <- c(0, cumsum(w))[j + 1]; wa_hi <- sum(w * a) - c(0, cumsum(w * a))[j + 1]
  pmax(w_lo * (top - p) + (sum(w) - w_lo) * top - wa_hi, 0)
}

## v_i = sum_b W_b t_b(c_i) / sum_b W_b over the region {lfdr < c_i} U {lfdr = c_i, P <= p_i}: a bin whose curve sits
## at c_i on [a_b, t_b] contributes p_i clamped to that stretch, the SNP's own bin at least its own p
.v_from_curves <- function(cur, c, p, grp, exact_below = 1e-2, n_grid = 5000L, credit = TRUE) {
  Wn <- cur$W / sum(cur$W); B <- length(Wn); n <- length(c)
  exact <- function(idx) {                                      # the per-bin evaluation, on a subset of SNPs
    v <- numeric(length(idx)); ci <- c[idx]; pi <- p[idx]; gi <- grp[idx]
    for (b in seq_len(B)) {
      a <- .a_of_c(cur$lfdr[, b], cur$s_p, ci)
      tb <- pmin(.t_of_c(cur$lfdr[, b], cur$s_p, ci, credit), pmax(a, pi))
      own <- gi == cur$ids[b]
      tb[own] <- pmax(a[own], pi[own])
      v <- v + Wn[b] * tb
    }
    v
  }
  if (n <= n_grid) return(exact(seq_len(n)))
  # pass 1: V(c) = sum_b Wn_b t_b(c) tabulated on a log grid in c and interpolated (the bulk), the own bin's exact term
  # in place of its share; at the cap every bin's stretch clamped at p (exact). pass 2: the exact per-bin evaluation
  # for every other SNP with v at most exact_below
  lo <- max(min(c[c > 0], na.rm = TRUE), 1e-300)                # log grid to 1, and a linear one on [.5, 1] where V is steep
  cg <- sort(unique(c(exp(seq(log(lo) - 1e-9, 0, length.out = n_grid)), seq(.5, 1, length.out = n_grid)))); G <- length(cg)
  V <- numeric(G)
  for (b in seq_len(B)) V <- V + Wn[b] * .t_of_c(cur$lfdr[, b], cur$s_p, cg, credit)
  a1 <- vapply(seq_len(B), function(b) .a_of_c(cur$lfdr[, b], cur$s_p, 1), numeric(1))
  V1 <- V[G]; V[G] <- sum(Wn * a1)                              # V just below the cap: its jump is not spread over the last interval
  lc <- log(pmin(pmax(c, lo), 1)); i <- pmin(pmax(findInterval(lc, log(cg)), 1L), G - 1L)
  w <- pmin(pmax((lc - log(cg)[i]) / (log(cg)[i + 1L] - log(cg)[i]), 0), 1)
  vlo <- V[i]                                                   # V at the grid point below c: a lower bound on v
  v <- V[i] * (1 - w) + V[i + 1L] * w
  top <- pnorm(cur$s_p[length(cur$s_p)])
  cap <- which(c >= 1)
  if (length(cap)) v[cap] <- V1 - .cap_share(a1, Wn, top, p[cap])
  idx <- split(seq_along(grp), factor(grp, levels = cur$ids))  # each bin's SNPs, in one pass
  for (b in seq_len(B)) {                                       # the own bin: its exact term in place of its interpolated share
    own <- idx[[b]]; if (!length(own)) next
    Lb <- cur$lfdr[, b]; ib <- i[own]
    t_hi <- .t_of_c(Lb, cur$s_p, cg[ib + 1L], credit); t_hi[ib + 1L == G] <- a1[b]
    ti <- (1 - w[own]) * .t_of_c(Lb, cur$s_p, cg[ib], credit) + w[own] * t_hi
    v[own] <- v[own] + Wn[b] * ifelse(c[own] < 1, pmax(.a_of_c(Lb, cur$s_p, c[own]), p[own]) - ti, pmax(0, p[own] - top))
  }
  v <- pmin(pmax(v, 0), 1)                                      # the grid value plus the own-bin floor can overshoot at a plateau
  redo <- which(vlo <= exact_below & c < 1)                     # every other SNP whose exact v can be at most exact_below
  if (length(redo)) v[redo] <- exact(redo)
  v
}

## Fit-only surface for a training fold: the steps of sffdr() up to the kernel fit (seed, surrogate ranks with random
## ties, perfect-LD copies, the fit on the tail and the downsampled nulls), the surrogate bins and nothing evaluated on
## the training SNPs. Same RNG order as sffdr(), so the fit is the one sffdr(keep_fit = TRUE) would keep.
.fit_surface <- function(p, fpi0, surrogate = fpi0, seed = 2026, nn = NULL, epsilon = .Machine$double.xmin, ...) {
  if (!is.null(seed)) withr::local_seed(seed)
  n <- length(p)
  z <- rank(surrogate, ties.method = "random") / n
  copies <- .ld_copies(p)
  xp <- cbind(z, p)
  kd <- kernelEstimator(x = if (any(copies)) xp[!copies, , drop = FALSE] else xp, eval.points = xp[1:2, , drop = FALSE],
                        nn = nn, epsilon = epsilon, weights = NULL, verbose = FALSE, ...)
  n_bins <- as.integer(max(10, min(1000, floor(n / 1000))))
  groups <- as.integer(cut(z, breaks = n_bins, labels = FALSE))
  list(lfit = attr(kd, "lfit"), z = z, groups = groups, n_bins = n_bins, p = p, fpi0 = fpi0, surrogate = surrogate)
}

## Held-out SNPs on a kept fit: surrogate rank by mid-rank in the training
## surrogate, bin by the training breaks, lfdr = own fpi0 over the bin's
## monotone density curve at the SNP's p.
.place_on_fit <- function(fit, cur, p_new, fpi0_new, surrogate_new = fpi0_new) {
  s_tr <- sort(fit$surrogate); n_tr <- length(s_tr)
  lo <- findInterval(surrogate_new, s_tr, left.open = TRUE); hi <- findInterval(surrogate_new, s_tr)
  tied <- hi > lo                                     # value tied with training surrogates: a random rank among them,
  z_new <- (lo + hi) / (2 * n_tr)                     # as the training ranks were (ties.method = "random"); else mid-rank
  if (any(tied)) z_new[tied] <- (lo[tied] + ceiling(runif(sum(tied)) * (hi[tied] - lo[tied]))) / n_tr
  z_new <- pmin(pmax(z_new, 1 / n_tr), 1)
  grp <- pmin(pmax(findInterval(z_new, cur$breaks), 1L), fit$n_bins)
  sp <- qnorm(pmin(pmax(p_new, .Machine$double.xmin), 1))
  fx <- rep(NA_real_, length(p_new))
  idx <- split(seq_along(grp), factor(grp, levels = cur$ids))  # each bin's SNPs, in one pass
  for (b in seq_along(cur$ids)) {
    i <- idx[[b]]
    if (length(i)) fx[i] <- approx(cur$s_p, cur$fx[, b], xout = sp[i], rule = 2)$y
  }
  miss <- is.na(fx)                                   # a bin empty in training: nearest bin's curve
  if (any(miss)) {
    for (i in which(miss)) { b <- which.min(abs(cur$ids - grp[i])); fx[i] <- approx(cur$s_p, cur$fx[, b], xout = sp[i], rule = 2)$y; grp[i] <- cur$ids[b] }
  }
  list(grp = grp, fx = fx, lfdr = pmin(fpi0_new / pmax(fx, .Machine$double.xmin), 1))
}
