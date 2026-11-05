#include <Rcpp.h>
using namespace Rcpp;

// The two loops over the surrogate bins in .v_from_curves (R/vvalues.R). The arithmetic of .t_of_c and .a_of_c is
// kept in the same order and each product is stored before it is added, so the result is the R loop's to the bit.
// L: grid x bins lfdr curves, non-decreasing in p; s: the probit grid.

static inline double t_one(const double* L, const double* s, int G, int i, double c, double base, double top) {
  if (i >= G) return top;
  if (i < 1) return base;
  double L1 = L[i - 1], L2 = L[i];
  double w = L2 > L1 ? (c - L1) / (L2 - L1) : 0;
  double ds = s[i] - s[i - 1];
  volatile double wd = w * ds;
  return R::pnorm(s[i - 1] + wd, 0.0, 1.0, 1, 0);
}

static inline double a_one(const double* L, const double* s, int G, int i, double c, double bot, double top) {
  if (i >= G) return top;
  if (i < 1) return bot;
  double L1 = L[i - 1], L2 = L[i];
  double r = (c - L1) / (L2 - L1);
  double ds = s[i] - s[i - 1];
  volatile double rd = r * ds;
  return R::pnorm(s[i - 1] + rd, 0.0, 1.0, 1, 0);
}

// #{L <= c}, or #{L < c} with strict (findInterval, and with left.open = TRUE)
static inline int n_below(const double* L, int G, double c, bool strict) {
  int lo = 0, hi = G;
  while (lo < hi) {
    int mid = (lo + hi) / 2;
    if (strict ? L[mid] < c : L[mid] <= c) lo = mid + 1; else hi = mid;
  }
  return lo;
}

// V(c) = sum_b W_b t_b(c) on a grid of c in increasing order
// [[Rcpp::export(name = ".v_grid")]]
NumericVector v_grid(NumericMatrix L, NumericVector s, NumericVector cg, NumericVector W, bool credit) {
  int G = L.nrow(), B = L.ncol(), n = cg.size();
  NumericVector V(n);
  double bot = R::pnorm(s[0], 0.0, 1.0, 1, 0), top = R::pnorm(s[G - 1], 0.0, 1.0, 1, 0), base = credit ? bot : 0;
  for (int b = 0; b < B; b++) {
    const double* Lb = &L(0, b);
    int i = 0;
    for (int k = 0; k < n; k++) {
      while (i < G && Lb[i] <= cg[k]) i++;
      volatile double wt = W[b] * t_one(Lb, s.begin(), G, i, cg[k], base, top);
      V[k] = V[k] + wt;
    }
  }
  return V;
}

// v_i = sum_b W_b min(t_b(c_i), max(a_b(c_i), p_i)), the SNP's own bin (column own_i, -1 for none) max(a_b(c_i), p_i)
// [[Rcpp::export(name = ".v_exact")]]
NumericVector v_exact(NumericMatrix L, NumericVector s, NumericVector c, NumericVector p, IntegerVector own, NumericVector W, bool credit) {
  int G = L.nrow(), B = L.ncol(), n = c.size();
  NumericVector v(n);
  double bot = R::pnorm(s[0], 0.0, 1.0, 1, 0), top = R::pnorm(s[G - 1], 0.0, 1.0, 1, 0), base = credit ? bot : 0;
  for (int j = 0; j < n; j++) {
    double acc = 0;
    for (int b = 0; b < B; b++) {
      const double* Lb = &L(0, b);
      double a = a_one(Lb, s.begin(), G, n_below(Lb, G, c[j], true), c[j], bot, top);
      double ap = a > p[j] ? a : p[j], tb = ap;
      if (b != own[j]) {
        double t = t_one(Lb, s.begin(), G, n_below(Lb, G, c[j], false), c[j], base, top);
        if (t < ap) tb = t;
      }
      volatile double wt = W[b] * tb;
      acc = acc + wt;
    }
    v[j] = acc;
  }
  return v;
}

// The own bin's exact term in place of its interpolated share in the c-grid pass: i (R's 1-based) and w the SNP's grid
// interval and weight, a1 each bin's a_b(1), own as in v_exact; adds W_b times the term to v
// [[Rcpp::export(name = ".v_own")]]
NumericVector v_own(NumericMatrix L, NumericVector s, NumericVector cg, IntegerVector i, NumericVector w, NumericVector c,
                    NumericVector p, IntegerVector own, NumericVector W, NumericVector a1, double top, bool credit, NumericVector v) {
  int G = L.nrow(), Gc = cg.size(), n = c.size();
  NumericVector out = clone(v);
  double bot = R::pnorm(s[0], 0.0, 1.0, 1, 0), tg = R::pnorm(s[G - 1], 0.0, 1.0, 1, 0), base = credit ? bot : 0;
  for (int j = 0; j < n; j++) {
    int b = own[j];
    if (b < 0) continue;
    const double* Lb = &L(0, b);
    int ib = i[j];
    double t_lo = t_one(Lb, s.begin(), G, n_below(Lb, G, cg[ib - 1], false), cg[ib - 1], base, tg);
    double t_hi = ib + 1 == Gc ? a1[b] : t_one(Lb, s.begin(), G, n_below(Lb, G, cg[ib], false), cg[ib], base, tg);
    volatile double x1 = (1 - w[j]) * t_lo;
    volatile double x2 = w[j] * t_hi;
    double ti = x1 + x2, term;
    if (c[j] < 1) {
      double a = a_one(Lb, s.begin(), G, n_below(Lb, G, c[j], true), c[j], bot, tg);
      term = (a > p[j] ? a : p[j]) - ti;
    } else {
      double d = p[j] - top;
      term = d > 0 ? d : 0;
    }
    volatile double wt = W[b] * term;
    out[j] = out[j] + wt;
  }
  return out;
}
