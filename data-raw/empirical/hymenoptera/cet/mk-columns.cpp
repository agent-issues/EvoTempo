// Hymenoptera CET pipeline: per-column Mk log-likelihoods on a TempoTree, for
// one partition of characters with K states and a symmetric rate matrix
// Q = V diag(ev) V' (uniform stationary distribution, so a uniform root).
// Used by HymMk() in common.R, which builds V and ev for unordered (JC) and
// ordered (adjacent changes only) characters and adds the Mkv correction
// per missing-data mask from extra constant columns.  MkLikelihood() in
// EvoTempo covers only unordered characters with one state or NA per cell;
// this kernel takes any 0/1 tip partials, so ambiguity sets such as {01} and
// ordered characters are scored as RevBayes scores them.
//
// Compiled by Rcpp::sourceCpp() from common.R.  ColumnLogLikR() in common.R
// is the pure-R reference (check-mkv.R checks the two agree).
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

// [[Rcpp::export]]
Rcpp::NumericVector MkColumnLogLik(const Rcpp::IntegerVector& postorder,
                                   const Rcpp::IntegerMatrix& kids, int nTip, int root,
                                   const Rcpp::NumericVector& lengths,
                                   const Rcpp::NumericMatrix& V, const Rcpp::NumericVector& ev,
                                   const Rcpp::NumericVector& tip, int nCol) {
  // tip: array [K, nCol, nTip] of tip partial likelihoods (0/1)
  const int K = V.nrow();
  const int nInternal = kids.nrow();
  if (tip.size() != static_cast<R_xlen_t>(K) * nCol * nTip) {
    Rcpp::stop("`tip` must have K x nCol x nTip elements");
  }
  const size_t block = static_cast<size_t>(nCol) * K;
  std::vector<double> part(static_cast<size_t>(nInternal) * block);
  std::vector<double> logScale(nCol, 0.0);
  std::vector<double> acc(block), P(K * K), e(K);
  for (R_xlen_t ii = 0; ii < postorder.size(); ++ii) {
    const int v = postorder[ii];  // 1-based internal node
    std::fill(acc.begin(), acc.end(), 1.0);
    for (int side = 0; side < 2; ++side) {
      const int kid = kids(v - nTip - 1, side);  // 1-based
      const double t = lengths[kid - 1];
      for (int m = 0; m < K; ++m) e[m] = std::exp(ev[m] * t);
      // P[a, b] = P(a -> b) = sum_m V[a, m] e[m] V[b, m]
      for (int a = 0; a < K; ++a) {
        for (int b = 0; b < K; ++b) {
          double s = 0.0;
          for (int m = 0; m < K; ++m) s += V(a, m) * e[m] * V(b, m);
          P[a * K + b] = s;
        }
      }
      const double* L = kid <= nTip ? &tip[static_cast<size_t>(kid - 1) * block]
                                    : &part[static_cast<size_t>(kid - nTip - 1) * block];
      for (int c = 0; c < nCol; ++c) {
        const double* l = L + static_cast<size_t>(c) * K;
        double* o = &acc[static_cast<size_t>(c) * K];
        for (int a = 0; a < K; ++a) {
          double s = 0.0;
          for (int b = 0; b < K; ++b) s += P[a * K + b] * l[b];
          o[a] *= s;
        }
      }
    }
    double* out = &part[static_cast<size_t>(v - nTip - 1) * block];
    for (int c = 0; c < nCol; ++c) {
      const double* a = &acc[static_cast<size_t>(c) * K];
      double mx = 0.0;
      for (int s = 0; s < K; ++s) mx = std::max(mx, a[s]);
      if (mx > 0.0) {
        for (int s = 0; s < K; ++s) out[c * K + s] = a[s] / mx;
        logScale[c] += std::log(mx);
      } else {
        for (int s = 0; s < K; ++s) out[c * K + s] = 0.0;
        logScale[c] = R_NegInf;
      }
    }
  }
  Rcpp::NumericVector res(nCol);
  const double* r = &part[static_cast<size_t>(root - nTip - 1) * block];
  for (int c = 0; c < nCol; ++c) {
    double s = 0.0;
    for (int a = 0; a < K; ++a) s += r[c * K + a];
    res[c] = std::log(s / K) + logScale[c];
  }
  return res;
}
