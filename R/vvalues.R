## ============================================================================
## vvalues.R -- cross-fitting behind sffdr(folds = ...): the density's tail
## (p < .1) refitted leave-one-fold-out, the bulk one fit (cross-fitting as in
## Liley and Wallace 2021):
## v_i = sum_b W_b t_b(lfdr_i) / sum_b W_b over the surrogate bins b, with W_b
## the null mass of the bin and t_b(c) the largest p at which its lfdr curve
## reaches c. A tie in the lfdr (the cap at 1, a flat stretch of a curve) is
## broken by p inside the region. Below the fold's smallest training p the
## ranking falls in proportion to p (the local FDR stays the fit's). The
## self-fit version (vvalues()) is internal, for tests.
## ============================================================================

## Per fold: refit the tail on the other folds, take the bulk from the one fit,
## place the held-out SNPs on the bin curves. pi0_folds: each fold's pi0, and
## its surrogate when none is given.
## folds_run restricts the run to some folds (job arrays); .merge_folds pools.
#' @importFrom stats approx predict qnorm pnorm dnorm
#' @noRd
.sffdr_cross <- function(p, fpi0, surrogate, folds, ncores = 1L, folds_run = NULL, grid_size = 600L, seed = 2026, verbose = TRUE, tail = .1,
                         pi0_folds = NULL, ...) {
  n <- length(p)
  fold_ids <- sort(unique(folds[!is.na(folds)]))
  if (length(fold_ids) < 2L) stop("'folds' must have at least two distinct labels.")
  run_ids <- if (is.null(folds_run)) fold_ids else intersect(fold_ids, folds_run)
  if (!length(run_ids)) stop("'folds_run' matches no fold label.")
  ok <- !is.na(p) & !is.na(fpi0) & !is.na(folds); if (!is.null(surrogate)) ok <- ok & !is.na(surrogate)
  P <- p[ok]; F0 <- fpi0[ok]; FO <- folds[ok]; p_min <- min(P); own_s <- is.null(surrogate)
  S <- if (own_s) .fpi0_of_fold(pi0_folds, NULL, rows = which(ok)) else surrogate[ok]
  if (verbose) message("  Fitting the bulk on every SNP")
  st <- .cross_setup(P, F0, S, p_min = p_min, grid_size = grid_size, seed = seed, tail = tail, ...)
  one_fold <- function(f) {
    tr <- FO != f
    if (!is.null(pi0_folds)) { F0 <- .fpi0_of_fold(pi0_folds, f, rows = which(ok)); if (own_s) S <- F0 }
    k <- .cross_fold(P, F0, S, tr, st, p_min = p_min, grid_size = grid_size, seed = seed, tail = tail, own_s = own_s, ...)
    v <- .v_from_curves(k$rk, c = k$c, p = P[!tr], grp = k$pl$grp, credit = min(P[tr]) > p_min)
    list(idx = which(ok)[!tr], lfdr = k$pl$lfdr, v = v, fx = k$pl$fx)
  }
  if (verbose) message(sprintf("  Cross-fitting the density: %d folds%s", length(fold_ids), if (ncores > 1L) sprintf(", %d cores", as.integer(ncores)) else ""))
  par <- ncores > 1L && .Platform$OS.type != "windows"            # folds in batches of ncores
  batches <- split(run_ids, ceiling(seq_along(run_ids) / if (par) as.integer(ncores) else 1L))
  res <- unlist(.lapply_pb(batches, function(ids) if (par) parallel::mclapply(ids, one_fold, mc.cores = as.integer(ncores)) else lapply(ids, one_fold),
                           verbose = verbose), recursive = FALSE)
  bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r$idx), logical(1))
  if (any(bad)) stop(sprintf("fold(s) %s failed", paste(run_ids[bad], collapse = ", ")))
  out <- list(flfdr = rep(NA_real_, n), fpvalues = rep(NA_real_, n), fx = rep(NA_real_, n))
  for (r in res) { out$flfdr[r$idx] <- r$lfdr; out$fpvalues[r$idx] <- r$v; out$fx[r$idx] <- r$fx }
  out
}

## Once per run: the bulk curves and the surrogate's order (ties at random)
.cross_setup <- function(P, F0, S, p_min, grid_size = 600L, seed = 2026, tail = .1, epsilon = .Machine$double.xmin, ...) {
  bulk <- .bulk_curves(.fit_surface(P, F0, S, seed = seed, epsilon = epsilon, ...), p_min = p_min, grid_size = grid_size, tail = tail)
  o <- if (is.null(seed)) order(S, runif(length(S))) else withr::with_seed(seed, order(S, runif(length(S))))
  list(bulk = bulk, o = o, s2 = qnorm(pmin(pmax(P, epsilon), 1 - 1e-4)), copies = .ld_copies(P))
}

## One fold (tr: training SNPs). own_s: S is the fold's own surrogate, ranked
## afresh in the training SNPs
.cross_fold <- function(P, F0, S, tr, st, p_min, grid_size = 600L, seed = 2026, tail = .1, epsilon = .Machine$double.xmin, own_s = FALSE, ...) {
  if (!is.null(seed)) withr::local_seed(seed)        # the fit reproducible per fold
  n_tr <- sum(tr)
  if (own_s) z <- rank(S[tr], ties.method = "random") / n_tr
  else { pos <- integer(length(tr)); pos[st$o] <- cumsum(tr[st$o]); z <- pos[tr] / n_tr }     # running rank
  fit <- .fit_tail(P[tr], F0[tr], z, st$s2[tr], st$copies[tr], seed = seed, epsilon = epsilon, tail = tail, ...)
  cur <- .bin_curves(fit, p_min = p_min, grid_size = grid_size, bulk = st$bulk)
  if (own_s) { fit$surrogate <- S[tr]; pl <- .place_on_fit(fit, cur, p_new = P[!tr], fpi0_new = F0[!tr], surrogate_new = S[!tr]) }
  else pl <- .place_on_fit(fit, cur, p_new = P[!tr], fpi0_new = F0[!tr], z_new = pos[!tr] / n_tr)
  rk <- .continue_curves(fit, cur)
  list(fit = fit, cur = cur, pl = pl, rk = rk, c = if (identical(rk, cur)) pl$lfdr else .on_curves(rk, pl$grp, P[!tr], F0[!tr]))
}

## A fold's tail fit (rows drawn as in kernelEstimator) and its bins
.fit_tail <- function(p, fpi0, z, s2, copies, seed = 2026, nn = NULL, epsilon = .Machine$double.xmin, tail = .1, target_null = 1e5,
                      maxit = 200, maxk = 5000, ...) {
  if (!is.null(seed)) withr::local_seed(seed)
  tl <- which(!copies & s2 < -2); nl <- which(!copies & s2 >= -2); bl <- if (length(nl) > target_null) sample(nl, target_null) else nl
  rows <- c(tl, bl)
  ts <- cbind(qnorm(pmin(pmax(z[rows], epsilon), 1 - 1e-4)), s2[rows] + runif(length(rows), -1e-6, 1e-6))
  lfit <- fit_bivariate_density(ts, tail_threshold = -2, target_null = length(bl), nn = nn, maxit = maxit, maxk = maxk,
                                weights = c(rep(1, length(tl)), rep(length(nl) / max(1, length(bl)), length(bl))), verbose = FALSE,
                                tail_box = qnorm(tail), ...)
  n_bins <- as.integer(max(10, min(1000, floor(length(z) / 1000)))); groups <- as.integer(cut(z, breaks = n_bins, labels = FALSE))
  cnt <- tabulate(groups, n_bins); k <- which(cnt > 0); sm <- rowsum(cbind(z, fpi0), groups, reorder = TRUE)
  list(lfit = lfit, z = z, groups = groups, n_bins = n_bins, p = p, fpi0 = fpi0, scale = attr(lfit, "scale"),
       bins = list(ids = k, zb = sm[, 1] / cnt[k], W = sm[, 2], pi0b = sm[, 2] / cnt[k]))
}

## nn = "cv": the window among multiples of the default (~5000 fitted rows)
## with the highest held-out log-likelihood of the tail (p < 0.023, the rows
## the kernel fits without downsampling) given the surrogate bin
.cv_nn <- function(p, fpi0, surrogate, folds, ncores = 1L, seed = 2026, epsilon = .Machine$double.xmin, mult = c(.25, .5, 1, 2, 4),
                   grid_size = 600L, verbose = TRUE, tail = .1, pi0_folds = NULL, ...) {
  ok <- !is.na(p) & !is.na(fpi0) & !is.na(folds); if (!is.null(surrogate)) ok <- ok & !is.na(surrogate)
  P <- p[ok]; F0 <- fpi0[ok]; FO <- folds[ok]; p_min <- min(P); fold_ids <- sort(unique(FO)); pt <- pnorm(-2); own_s <- is.null(surrogate)
  S <- if (own_s) .fpi0_of_fold(pi0_folds, NULL, rows = which(ok)) else surrogate[ok]
  n_tail <- sum(P < pt); nn0 <- max(0.005, min(0.05, 5000 / (n_tail + min(length(P) - n_tail, 1e5))))
  st <- .cross_setup(P, F0, S, p_min = p_min, grid_size = grid_size, seed = seed, tail = tail, epsilon = epsilon, ...)
  one_fold <- function(f, nn) {                       # scores the tail only (p < .023)
    tr <- FO != f
    if (own_s) S <- .fpi0_of_fold(pi0_folds, f, rows = which(ok))
    k <- .cross_fold(P, F0, S, tr, st, p_min = p_min, grid_size = grid_size, seed = seed, tail = tail, nn = nn, epsilon = epsilon, own_s = own_s, ...)
    s <- pnorm(k$cur$s_p); tl <- P[!tr] < pt
    mass <- apply(k$cur$fx, 2, function(fx) approx(s, cumsum(c(fx[1] * s[1], diff(s) * (fx[-1] + fx[-length(fx)]) / 2)), xout = pt, rule = 2)$y)
    sum(log(k$pl$fx[tl] / mass[match(k$pl$grp[tl], k$cur$ids)]))
  }
  cand <- nn0 * mult
  if (verbose) message(sprintf("  Choosing nn: %d windows x %d folds", length(cand), length(fold_ids)))
  ll <- vapply(cand, function(nn) {
    r <- if (ncores > 1L && .Platform$OS.type != "windows") parallel::mclapply(fold_ids, function(f) try(one_fold(f, nn), silent = TRUE), mc.cores = as.integer(ncores))
         else lapply(fold_ids, function(f) try(one_fold(f, nn), silent = TRUE))
    if (any(vapply(r, inherits, logical(1), "try-error"))) -Inf else sum(unlist(r))
  }, numeric(1))
  if (!any(is.finite(ll))) stop("nn = \"cv\": the density fit failed at every window.")
  if (verbose) message(sprintf("  nn = %.4g (held-out tail log-likelihood %s)", cand[which.max(ll)], paste(sprintf("%.1f", ll), collapse = ", ")))
  cand[which.max(ll)]
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
.bin_curves <- function(fit, p_min, grid_size = 600L, bulk = NULL) {
  s_p <- seq(qnorm(max(p_min, .Machine$double.xmin)), qnorm(1 - 1e-4), length.out = grid_size)
  if (is.null(fit$bins)) {
    zb <- tapply(fit$z, fit$groups, mean)
    W <- tapply(fit$fpi0, fit$groups, sum)
    pi0b <- tapply(fit$fpi0, fit$groups, mean)
    ids <- as.integer(names(zb))
  } else {
    zb <- fit$bins$zb; W <- fit$bins$W; pi0b <- fit$bins$pi0b; ids <- fit$bins$ids
  }
  if (is.null(bulk)) {
    nd <- data.frame(V1 = rep(qnorm(as.numeric(zb)), each = grid_size), V2 = rep(s_p, times = length(zb)))
    fs <- predict(fit$lfit, newdata = nd, log = FALSE)
    fx <- fs / pmax(dnorm(nd$V1) * dnorm(nd$V2), .Machine$double.xmin)
    fx <- matrix(pmax(fx, .Machine$double.xmin), nrow = grid_size)
  } else {                                            # the tail refit below the splice, the one fit above
    V1 <- qnorm(as.numeric(zb)); lo <- s_p < bulk$s; fs <- matrix(NA_real_, grid_size, length(zb))
    Vq <- if (length(V1) > 250L) seq(min(V1), max(V1), length.out = 250L) else V1      # smooth across the surrogate
    fq <- matrix(fit$scale * predict(fit$lfit, newdata = data.frame(V1 = rep(Vq, each = sum(lo)), V2 = rep(s_p[lo], times = length(Vq))), log = FALSE), sum(lo))
    fs[lo, ] <- if (length(Vq) == length(V1)) fq else t(apply(fq, 1, function(r) approx(Vq, r, xout = V1)$y))
    fs[!lo, ] <- t(apply(bulk$fs, 1, function(r) approx(bulk$V1, r, xout = V1, rule = 2)$y))
    fx <- pmax(fs / pmax(outer(dnorm(s_p), dnorm(V1)), .Machine$double.xmin), .Machine$double.xmin)
  }
  fx <- apply(fx, 2, cummin)
  L <- sweep(fx, 2, as.numeric(pi0b), FUN = function(f, p0) pmin(p0 / f, 1))
  rz <- range(fit$z); dz <- diff(rz)
  list(s_p = s_p, fx = fx, lfdr = L, W = as.numeric(W), ids = ids,
       breaks = seq(rz[1] - dz / 1000, rz[2] + dz / 1000, length.out = fit$n_bins + 1))
}

## The one fit above the splice (p >= tail): every fold's bulk
.bulk_curves <- function(fit, p_min, grid_size = 600L, tail = .1) {
  s_p <- seq(qnorm(max(p_min, .Machine$double.xmin)), qnorm(1 - 1e-4), length.out = grid_size); hi <- s_p >= qnorm(tail)
  V1 <- qnorm(as.numeric(tapply(fit$z, fit$groups, mean)))
  fs <- predict(fit$lfit, newdata = data.frame(V1 = rep(V1, each = sum(hi)), V2 = rep(s_p[hi], times = length(V1))), log = FALSE)
  list(s = qnorm(tail), V1 = V1, fs = matrix(fs, nrow = sum(hi)))
}

## The ranking: below the fold's smallest training p (e) each bin's density
## ratio is the fit's at e (before the running minimum) times e / p, above e
## the fit's curve
.continue_curves <- function(fit, cur) {
  e <- min(fit$p); lo <- cur$s_p < qnorm(e); j <- sum(lo) + 1L
  if (j == 1L || j > length(lo)) return(cur)
  zb <- qnorm(as.numeric(if (is.null(fit$bins)) tapply(fit$z, fit$groups, mean) else fit$bins$zb))
  fe <- (if (is.null(fit$scale)) 1 else fit$scale) * predict(fit$lfit, newdata = data.frame(V1 = zb, V2 = cur$s_p[j]), log = FALSE)
  fe <- fe / pmax(dnorm(zb) * dnorm(cur$s_p[j]), .Machine$double.xmin)
  cur$fx[lo, ] <- outer(e / pnorm(cur$s_p[lo]), pmax(fe, cur$fx[j, ]))
  pi0b <- if (is.null(fit$bins)) tapply(fit$fpi0, fit$groups, mean) else fit$bins$pi0b
  cur$lfdr[lo, ] <- sweep(cur$fx[lo, , drop = FALSE], 2, as.numeric(pi0b), FUN = function(f, p0) pmin(p0 / f, 1))
  cur
}

## lfdr of SNPs placed in bins grp, read off cur
.on_curves <- function(cur, grp, p, fpi0) {
  sp <- qnorm(pmin(pmax(p, .Machine$double.xmin), 1)); fx <- rep(NA_real_, length(p))
  idx <- split(seq_along(grp), factor(grp, levels = cur$ids))
  for (b in seq_along(cur$ids)) {
    i <- idx[[b]]
    if (length(i)) fx[i] <- approx(cur$s_p, cur$fx[, b], xout = sp[i], rule = 2)$y
  }
  pmin(fpi0 / pmax(fx, .Machine$double.xmin), 1)
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
  own_b <- match(grp, cur$ids) - 1L; own_b[is.na(own_b)] <- -1L
  exact <- function(idx) .v_exact(cur$lfdr, cur$s_p, c[idx], p[idx], own_b[idx], Wn, credit)  # src/vsums.cpp
  if (n <= n_grid) return(exact(seq_len(n)))
  # pass 1: V(c) = sum_b Wn_b t_b(c) tabulated on a log grid in c and interpolated (the bulk), the own bin's exact term
  # in place of its share; at the cap every bin's stretch clamped at p (exact). pass 2: the exact per-bin evaluation
  # for every other SNP with v at most exact_below
  lo <- max(min(c[c > 0], na.rm = TRUE), 1e-300)                # log grid to 1, and a linear one on [.5, 1] where V is steep
  cg <- sort(unique(c(exp(seq(log(lo) - 1e-9, 0, length.out = n_grid)), seq(.5, 1, length.out = n_grid)))); G <- length(cg)
  V <- .v_grid(cur$lfdr, cur$s_p, cg, Wn, credit)
  a1 <- vapply(seq_len(B), function(b) .a_of_c(cur$lfdr[, b], cur$s_p, 1), numeric(1))
  V1 <- V[G]; V[G] <- sum(Wn * a1)                              # V just below the cap: its jump is not spread over the last interval
  lc <- log(pmin(pmax(c, lo), 1)); i <- pmin(pmax(findInterval(lc, log(cg)), 1L), G - 1L)
  w <- pmin(pmax((lc - log(cg)[i]) / (log(cg)[i + 1L] - log(cg)[i]), 0), 1)
  vlo <- V[i]                                                   # V at the grid point below c: a lower bound on v
  v <- V[i] * (1 - w) + V[i + 1L] * w
  top <- pnorm(cur$s_p[length(cur$s_p)])
  cap <- which(c >= 1)
  if (length(cap)) v[cap] <- V1 - .cap_share(a1, Wn, top, p[cap])
  v <- .v_own(cur$lfdr, cur$s_p, cg, i, w, c, p, own_b, Wn, a1, top, credit, v)   # the own bin: its exact term in place of its share
  v <- pmin(pmax(v, 0), 1)                                      # the grid value plus the own-bin floor can overshoot at a plateau
  redo <- which(vlo <= exact_below & c < 1)                     # every other SNP whose exact v can be at most exact_below
  if (length(redo)) v[redo] <- exact(redo)
  v
}

## Fit-only surface for a training fold: the steps of sffdr() up to the kernel fit (seed, surrogate ranks with random
## ties, perfect-LD copies, the fit on the tail and the downsampled nulls), the surrogate bins and nothing evaluated on
## the training SNPs. Same RNG order as sffdr(), so the fit is the one sffdr(keep_fit = TRUE) would keep. With tail,
## only the surface below it is fitted.
.fit_surface <- function(p, fpi0, surrogate = fpi0, seed = 2026, nn = NULL, epsilon = .Machine$double.xmin, tail = NULL, z = NULL, ...) {
  if (!is.null(seed)) withr::local_seed(seed)
  n <- length(p)
  if (is.null(z)) z <- rank(surrogate, ties.method = "random") / n
  copies <- .ld_copies(p)
  xp <- cbind(z, p)
  kd <- kernelEstimator(x = if (any(copies)) xp[!copies, , drop = FALSE] else xp, eval.points = xp[1:2, , drop = FALSE],
                        nn = nn, epsilon = epsilon, weights = NULL, verbose = FALSE, tail_box = if (!is.null(tail)) qnorm(tail), ...)
  n_bins <- as.integer(max(10, min(1000, floor(n / 1000))))
  groups <- as.integer(cut(z, breaks = n_bins, labels = FALSE))
  lfit <- attr(kd, "lfit")
  out <- list(lfit = lfit, z = z, groups = groups, n_bins = n_bins, p = p, fpi0 = fpi0, surrogate = surrogate,
              scale = if (is.null(attr(lfit, "scale"))) 1 else attr(lfit, "scale"))
  if (!is.null(tail)) {                               # bin centres, null mass, mean pi0
    cnt <- tabulate(groups, n_bins); k <- which(cnt > 0); sm <- rowsum(cbind(z, fpi0), groups, reorder = TRUE)
    out$bins <- list(ids = k, zb = sm[, 1] / cnt[k], W = sm[, 2], pi0b = sm[, 2] / cnt[k])
  }
  out
}

## Held-out SNPs on a kept fit: surrogate rank by mid-rank in the training
## surrogate, bin by the training breaks, lfdr = own fpi0 over the bin's
## monotone density curve at the SNP's p.
.place_on_fit <- function(fit, cur, p_new, fpi0_new, surrogate_new = fpi0_new, z_new = NULL) {
  n_tr <- length(fit$z)
  if (is.null(z_new)) {
    s_tr <- sort(fit$surrogate)
    lo <- findInterval(surrogate_new, s_tr, left.open = TRUE); hi <- findInterval(surrogate_new, s_tr)
    tied <- hi > lo                                   # value tied with training surrogates: a random rank among them,
    z_new <- (lo + hi) / (2 * n_tr)                   # as the training ranks were (ties.method = "random"); else mid-rank
    if (any(tied)) z_new[tied] <- (lo[tied] + ceiling(runif(sum(tied)) * (hi[tied] - lo[tied]))) / n_tr
  }
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
