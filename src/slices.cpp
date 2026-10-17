// Frozen-E slices of the path proposal (.AddSlices() in R/Paths.R): per
// slice, E at the slice midpoint and the eigendecomposition of the
// symmetrized frozen generator, which is tridiagonal.

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include "frozen.h"
#include <algorithm>
#include <vector>
#ifndef FCONE
# define FCONE
#endif

using namespace Rcpp;

// `mids`: slice midpoints; `off`: off-diagonal of the symmetrized generator.
// Eigenvalues are returned in decreasing order, as by eigen(symmetric = TRUE);
// eigenvector signs may differ from eigen(), which no caller depends on.
// [[Rcpp::export(name = ".SlicesCpp")]]
List SlicesCpp(List set, NumericVector mids, NumericVector off) {
  const EvoTempo::FrozenSet fs(set, false);
  const int K = fs.K;
  const int nSlice = mids.size();
  List out(nSlice);
  std::vector<double> d(K), e(std::max(K - 1, 1)), z(static_cast<std::size_t>(K) * K),
    work(std::max(2 * K - 2, 1));
  for (int s = 0; s < nSlice; ++s) {
    NumericVector eBar(K), values(K);
    NumericMatrix vectors(K, K);
    for (int k = 0; k < K; ++k) {
      eBar[k] = fs.Hermite(mids[s], k, false);
      d[k] = fs.q(k, k) + 2 * fs.lambda[k] * eBar[k] - fs.sigma[k];
    }
    std::copy(off.begin(), off.end(), e.begin());
    int info = 0;
    F77_CALL(dstev)("V", &K, d.data(), e.data(), z.data(), &K, work.data(), &info
                    FCONE);
    if (info != 0) Rcpp::stop("dstev failed (info = %d)", info);  // # nocov
    // Ascending from dstev; reverse
    for (int k = 0; k < K; ++k) {
      const int src = K - 1 - k;
      values[k] = d[src];
      std::copy(z.begin() + static_cast<std::size_t>(src) * K,
                z.begin() + static_cast<std::size_t>(src + 1) * K,
                vectors.begin() + static_cast<std::size_t>(k) * K);
    }
    out[s] = List::create(_["eBar"] = eBar, _["values"] = values,
                          _["vectors"] = vectors);
  }
  return out;
}
