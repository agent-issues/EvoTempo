#ifndef EVOTEMPO_FROZEN_H
#define EVOTEMPO_FROZEN_H

// C++ view of the frozen-E proposal model (`set` from .Precompute() and
// .AddSlices() in R/Paths.R), and the helpers that R/Paths.R defines on it.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace EvoTempo {

// A proposal that cannot be drawn numerically; R callers reject it
struct InfeasibleProposal : public std::runtime_error {
  explicit InfeasibleProposal(const std::string& m) : std::runtime_error(m) {}
};

// findInterval(x, v) for sorted v: number of elements <= x
inline int FindInterval(double x, const std::vector<double>& v) {
  return static_cast<int>(std::upper_bound(v.begin(), v.end(), x) - v.begin());
}

struct Slice {
  std::vector<double> eBar, values, vectors;  // vectors: K x K column-major
};

struct FrozenSet {
  int K = 0;
  int nSlice = 0;
  int nAge = 0;
  bool joint = false;
  double removal = 0.0, rho = 1.0;
  std::string condition;
  std::vector<double> Q, lambda, sigma, psi, prior, survive, halfPi, bounds;
  std::vector<Slice> slices;
  std::vector<double> egAges, egE, egG, egDE;  // nAge x K column-major

  // withSlices = FALSE for a set from .Precompute() alone
  explicit FrozenSet(Rcpp::List set, bool withSlices = true) {
    const Rcpp::NumericMatrix q = set["Q"];
    K = q.nrow();
    Q.assign(q.begin(), q.end());
    lambda = Rep(set["lambda"]);
    sigma = Rep(set["sigma"]);
    psi = Rep(set["psi"]);
    prior = Rep(set["prior"]);
    survive = Rep(set["survive"]);
    joint = Rcpp::as<bool>(set["joint"]);
    removal = Rcpp::as<double>(set["removal"]);
    rho = Rcpp::as<double>(set["rho"]);
    condition = Rcpp::as<std::string>(set["condition"]);
    Rcpp::List eg = set["eg"];
    const Rcpp::NumericVector ages = eg["ages"];
    egAges.assign(ages.begin(), ages.end());
    nAge = egAges.size();
    const Rcpp::NumericMatrix e = eg["E"], g = eg["G"], de = eg["dE"];
    egE.assign(e.begin(), e.end());
    egG.assign(g.begin(), g.end());
    egDE.assign(de.begin(), de.end());
    if (withSlices) {
      const Rcpp::NumericVector hp = set["halfPi"], b = set["bounds"];
      halfPi.assign(hp.begin(), hp.end());
      bounds.assign(b.begin(), b.end());
      Rcpp::List sl = set["slices"];
      nSlice = sl.size();
      slices.resize(nSlice);
      for (int s = 0; s < nSlice; ++s) {
        Rcpp::List one = sl[s];
        const Rcpp::NumericVector eb = one["eBar"], va = one["values"];
        const Rcpp::NumericMatrix ve = one["vectors"];
        slices[s].eBar.assign(eb.begin(), eb.end());
        slices[s].values.assign(va.begin(), va.end());
        slices[s].vectors.assign(ve.begin(), ve.end());
      }
    }
  }

  std::vector<double> Rep(SEXP x) const {
    const Rcpp::NumericVector v(x);
    std::vector<double> out(K);
    for (int k = 0; k < K; ++k) out[k] = v[v.size() == 1 ? 0 : k];
    return out;
  }

  double q(int i, int j) const { return Q[i + static_cast<std::size_t>(j) * K]; }

  // .Hermite(eg, at, bin, value): G (with derivative E) or E (with dE);
  // bin is zero-based
  double Hermite(double at, int bin, bool G) const {
    int i = FindInterval(at, egAges);
    i = std::max(1, std::min(i, nAge - 1)) - 1;  // all.inside, zero-based
    const double h = egAges[i + 1] - egAges[i];
    const double u = (at - egAges[i]) / h;
    const std::vector<double>& y = G ? egG : egE;
    const std::vector<double>& dy = G ? egE : egDE;
    const std::size_t i0 = i + static_cast<std::size_t>(bin) * nAge, i1 = i0 + 1;
    const double u2 = u * u, u3 = u2 * u;
    return (2 * u3 - 3 * u2 + 1) * y[i0] + (u3 - 2 * u2 + u) * h * dy[i0] +
      (3 * u2 - 2 * u3) * y[i1] + (u3 - u2) * h * dy[i1];
  }

  // .TipValues() for one tip in bin k (zero-based)
  double TipValue(int type, double age, int k) const {
    if (type == 0) return rho;
    if (type == 1) return psi[k] * (removal + (1 - removal) * Hermite(age, k, false));
    return 0.0;
  }

  // .SliceExp(set, s, span), s zero-based; out K x K column-major
  void SliceExp(int s, double span, std::vector<double>& out) const {
    const Slice& sl = slices[s];
    std::vector<double> ex(K);
    for (int k = 0; k < K; ++k) ex[k] = std::exp(sl.values[k] * span);
    out.assign(static_cast<std::size_t>(K) * K, 0.0);
    for (int j = 0; j < K; ++j) {
      for (int i = 0; i < K; ++i) {
        double s0 = 0.0;
        for (int k = 0; k < K; ++k) {
          s0 += sl.vectors[i + k * K] * (ex[k] * sl.vectors[j + k * K]);
        }
        out[i + static_cast<std::size_t>(j) * K] = s0 * ((1 / halfPi[i]) * halfPi[j]);
      }
    }
  }
};

// Pieces of one branch (.Pieces): boundaries, old to young; one-based slice
// indices; frozen transfer matrices
struct Pieces {
  std::vector<double> edges;
  std::vector<int> slice;
  std::vector<std::vector<double> > M;
};

}  // namespace EvoTempo

#endif  // EVOTEMPO_FROZEN_H
