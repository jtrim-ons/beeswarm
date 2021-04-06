#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

int which_min_abs(std::vector<double> & y_best, std::vector<bool> & placed)
{
  int i = 0;
  while (placed[i]) ++i;
  double best_val = std::abs(y_best[i]);
  int result = i;
  for (++i; i<y_best.size(); ++i) {
    if (placed[i]) continue;
    double a = std::abs(y_best[i]);
    if (a < best_val) {
      best_val = a;
      result = i;
    }
  }
  return result;
}


// [[Rcpp::export]]
NumericVector fastSwarm(NumericVector x, LogicalVector compact)
{
  NumericVector y = NumericVector(x.size());
  std::vector<bool> placed(x.size());
  std::vector<double> y_low(x.size());
  std::vector<double> y_high(x.size());
  std::vector<double> y_best(x.size());
  
  LogicalVector x_is_na = is_na(x);
  int na_count = 0;
  for (int i=0; i<x.size(); i++) {
    if (x_is_na[i]) {
      y[i] = NumericVector::get_na();
      placed[i] = true;
      ++na_count;
    }
  }
  
  for (int iter=0; iter<x.size() - na_count; iter++) {
    int i = which_min_abs(y_best, placed);
    double xi = x[i];
    double yi = y_best[i];
    y[i] = yi;
    placed[i] = true;
    for (int j=0; j<x.size(); j++) {
      if (placed[j]) continue;
      double xdiff = std::abs(xi - x[j]);
      if (xdiff >= 1) continue;
      double y_offset = std::sqrt(1 - xdiff * xdiff);
      double y_hi = std::max(y_high[j], yi + y_offset);
      y_high[j] = y_hi;
      double y_lo = std::min(y_low[j], yi - y_offset);
      y_low[j] = y_lo;
      y_best[j] = -y_lo < y_hi ? y_lo : y_hi;
    }
  }
  return y;
}