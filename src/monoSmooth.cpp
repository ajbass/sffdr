#include <Rcpp.h>
using namespace Rcpp;

//' Fast Conditional Monotonic Smoothing
//' 
//' Enforces strictly non-increasing density conditional on surrogate bins.
//' Note: Data must be sorted by (group, pvalue) before passing.
//' 
//' @param pvalue Numeric vector of p-values.
//' @param density Numeric vector of estimated densities.
//' @param group Integer vector of surrogate bins.
//' @return A numeric vector of smoothed densities.
//' @export
// [[Rcpp::export]]
NumericVector monoSmooth_conditional(NumericVector pvalue, NumericVector density, IntegerVector group) {
  
  // Claude's memory safety fix: clone to prevent in-place mutation of R data
  NumericVector out = clone(density); 
  int n = pvalue.size();
  
  if (n == 0) return out;
  
  double current_min = out[0];
  int current_group = group[0];
  
  // Single O(N) pass
  for (int i = 1; i < n; ++i) {
    if (group[i] != current_group) {
      // We entered a new bin. Reset the running minimum.
      current_group = group[i];
      current_min = out[i];
    } else {
      // Same bin: enforce monotonicity
      if (out[i] > current_min) {
        out[i] = current_min; // Crush the peak
      } else {
        current_min = out[i]; // Update the new lowest seen density
      }
    }
  }
  
  return out;
}