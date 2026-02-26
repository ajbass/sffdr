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
#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

//' Fast Conditional PAVA (Decreasing)
//'  
//' @param pvalue A numeric vector of p-values. Must be sorted ascendingly within each group.
//' @param density A numeric vector of initial density estimates (e.g., from locfit).
//' @param group An integer vector denoting the surrogate bin/stratum for each observation.
//' 
//' @return A numeric vector of smoothed, monotonically non-increasing densities.
//' @export
// [[Rcpp::export]]
NumericVector monoSmooth_pava(NumericVector pvalue, NumericVector density, IntegerVector group) {
  int n = density.size();
  NumericVector out = clone(density);
  if (n <= 1) return out;
  
  std::vector<double> val(n);
  std::vector<int> weight(n);
  
  int current_group = group[0];
  int top = -1; 
  int start_idx = 0;
  
  for (int i = 0; i <= n; ++i) {
    bool group_change = (i == n || group[i] != current_group);
    
    if (i < n && !group_change) {
      // Add new element to the working stack
      top++;
      val[top] = density[i];
      weight[top] = 1;
      
      // PAVA Logic: Enforce decreasing constraint
      // If the previous value is LESS than the current value, pool them.
      while (top > 0 && val[top - 1] < val[top]) {
        double new_val = (val[top - 1] * weight[top - 1] + val[top] * weight[top]) / (double)(weight[top - 1] + weight[top]);
        weight[top - 1] += weight[top];
        val[top - 1] = new_val;
        top--; // Shrink the stack as blocks are merged
      }
    } else {
      // Group boundary reached: unpack the pooled blocks into output array
      int write_idx = start_idx;
      for (int b = 0; b <= top; ++b) {
        for (int w = 0; w < weight[b]; ++w) {
          out[write_idx++] = val[b];
        }
      }
      
      // Reset tracker for the next surrogate bin
      if (i < n) {
        current_group = group[i];
        start_idx = i;
        top = 0;
        val[top] = density[i];
        weight[top] = 1;
      }
    }
  }
  return out;
}