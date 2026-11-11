// Trilobite CET pipeline: pruning of the characters that MkLikelihood() cannot
// score (ambiguous cells; ordered characters) and of the constant columns of
// the Mkv corrections, compiled by common.R with Rcpp::sourceCpp().
//
// Columns have their own number of states k[c].  A column with ordId[c] < 0
// evolves under symmetric Mk (Jukes-Cantor with k[c] states); otherwise under
// the reversible generator whose orthonormal eigenvectors and eigenvalues are
// ordV[[ordId[c]]] and ordLam[[ordId[c]]] (P(t) = V diag(exp(lam t)) V').  Root
// distribution uniform.  Same recursion, per-node per-column rescaling and
// root sum as SymLogLik() and OrdLogLik() in common.R, which check-mkv.R
// compares with this function.
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector CetPrune(IntegerVector postorder, IntegerMatrix kids, int nTip, int root,
                       NumericVector lengths, NumericVector tip, IntegerVector k,
                       IntegerVector ordId, List ordV, List ordLam) {
  const int nCol = k.size();
  std::vector<int> off(nCol + 1, 0);
  for (int c = 0; c < nCol; ++c) off[c + 1] = off[c] + k[c];
  const int S = off[nCol];
  const int nAll = nTip + kids.nrow();
  if (tip.size() != (R_xlen_t) S * nTip) stop("`tip` must hold sum(k) values per tip");
  std::vector<double> part((size_t) S * nAll, 0.0);
  std::copy(tip.begin(), tip.end(), part.begin());
  std::vector<double> logScale(nCol, 0.0), val(S);
  const int nOrd = ordV.size();
  std::vector<NumericMatrix> V;
  std::vector<NumericVector> Lam;
  for (int o = 0; o < nOrd; ++o) {
    V.push_back(as<NumericMatrix>(ordV[o]));
    Lam.push_back(as<NumericVector>(ordLam[o]));
  }
  std::vector<std::vector<double> > P(nOrd);
  for (R_xlen_t pi = 0; pi < postorder.size(); ++pi) {
    const int v = postorder[pi] - 1;
    std::fill(val.begin(), val.end(), 1.0);
    for (int j = 0; j < 2; ++j) {
      const int kid = kids(v - nTip, j) - 1;
      const double t = lengths[kid];
      const double* Lk = &part[(size_t) kid * S];
      for (int o = 0; o < nOrd; ++o) {
        const int K = V[o].nrow();
        P[o].assign((size_t) K * K, 0.0);
        std::vector<double> e(K);
        for (int m = 0; m < K; ++m) e[m] = std::exp(Lam[o][m] * t);
        for (int a = 0; a < K; ++a) {
          for (int b = 0; b < K; ++b) {
            double s = 0.0;
            for (int m = 0; m < K; ++m) s += V[o](a, m) * e[m] * V[o](b, m);
            P[o][a + (size_t) K * b] = s;
          }
        }
      }
      for (int c = 0; c < nCol; ++c) {
        const int K = k[c];
        const double* l = Lk + off[c];
        double* out = &val[off[c]];
        if (ordId[c] < 0) {
          const double stay = 1.0 / K + (K - 1.0) / K * std::exp(-K * t / (K - 1.0));
          const double move = (1.0 - stay) / (K - 1.0);
          double sum = 0.0;
          for (int a = 0; a < K; ++a) sum += l[a];
          for (int a = 0; a < K; ++a) out[a] *= move * sum + (stay - move) * l[a];
        } else {
          const std::vector<double>& Pm = P[ordId[c]];
          for (int a = 0; a < K; ++a) {
            double s = 0.0;
            for (int b = 0; b < K; ++b) s += Pm[a + (size_t) K * b] * l[b];
            out[a] *= s;
          }
        }
      }
    }
    double* dst = &part[(size_t) v * S];
    for (int c = 0; c < nCol; ++c) {
      double mx = 0.0;
      for (int a = 0; a < k[c]; ++a) mx = std::max(mx, val[off[c] + a]);
      for (int a = 0; a < k[c]; ++a) dst[off[c] + a] = mx > 0 ? val[off[c] + a] / mx : 0.0;
      logScale[c] += std::log(mx);
    }
  }
  NumericVector out(nCol);
  const double* r = &part[(size_t) (root - 1) * S];
  for (int c = 0; c < nCol; ++c) {
    double s = 0.0;
    for (int a = 0; a < k[c]; ++a) s += r[off[c] + a];
    out[c] = std::log(s / k[c]) + logScale[c];
  }
  return out;
}
