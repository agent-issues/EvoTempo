// Mk likelihood of discrete characters on a TempoTree, by the MkPrime
// header-only pruning kernels (LinkingTo: MkPrime).  MkLikelihood() and
// .MkCharactersLogLik() in R/Characters.R hold the pure-R references.

#include <Rcpp.h>
#include <MkPrime/MkPrime.h>
#include "tree.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace Rcpp;

namespace {

// exp(Q t) by uniformization and squaring, as .MkTransition() in
// R/Characters.R (#98, R98-13): mu = max_i -Q_ii, B = I + Q / mu >= 0;
// x = mu t is halved s times until y = x / 2^s <= 1/2; exp(Q t / 2^s) =
// exp(-y) sum_{n <= K + 17} y^n / n! B^n; then s squarings.  Only
// non-negative terms are added, so small probabilities keep their relative
// accuracy, and P(0) = I exactly.  P is row-major: P[a * K + b] = P(a -> b).
inline std::vector<double> UniformizedTransition(const std::vector<double>& Q, int K,
                                                 double t) {
  const std::size_t KK = static_cast<std::size_t>(K) * K;
  std::vector<double> P(KK, 0.0);
  if (!(t > 0.0)) {
    for (int a = 0; a < K; ++a) P[a * K + a] = 1.0;
    return P;
  }
  double mu = 0.0;
  for (int a = 0; a < K; ++a) mu = std::max(mu, -Q[a * K + a]);
  const double x = mu * t;
  const int s = x > 0.5 ? static_cast<int>(std::ceil(std::log2(x / 0.5))) : 0;
  const double y = x / std::ldexp(1.0, s);
  std::vector<double> B(KK), term(KK, 0.0), next(KK), total(KK, 0.0);
  for (std::size_t i = 0; i < KK; ++i) B[i] = Q[i] / mu;
  for (int a = 0; a < K; ++a) {
    B[a * K + a] += 1.0;
    term[a * K + a] = 1.0;
    total[a * K + a] = 1.0;
  }
  for (int n = 1; n <= K + 17; ++n) {
    const double f = y / n;
    for (int a = 0; a < K; ++a) {
      for (int b = 0; b < K; ++b) {
        double acc = 0.0;
        for (int m = 0; m < K; ++m) acc += term[a * K + m] * B[m * K + b];
        next[a * K + b] = acc * f;
      }
    }
    term.swap(next);
    for (std::size_t i = 0; i < KK; ++i) total[i] += term[i];
  }
  const double e = std::exp(-y);
  for (std::size_t i = 0; i < KK; ++i) P[i] = e * total[i];
  for (int i = 0; i < s; ++i) {
    for (int a = 0; a < K; ++a) {
      for (int b = 0; b < K; ++b) {
        double acc = 0.0;
        for (int m = 0; m < K; ++m) acc += P[a * K + m] * P[m * K + b];
        next[a * K + b] = acc;
      }
    }
    P.swap(next);
  }
  return P;
}

// Ordered characters (or any rate matrix Q with a uniform root distribution)
// for MkPrime::Pruner, with P(t) from UniformizedTransition()
class UniformizedModel {
 public:
  UniformizedModel(int k, const double* Q)  // Q column-major, as in R
    : k_(k), Q_(static_cast<std::size_t>(k) * k) {
    for (int a = 0; a < k; ++a) {
      for (int b = 0; b < k; ++b) Q_[a * k + b] = Q[a + b * k];
    }
  }

  int n_states() const { return k_; }

  void transition(const double* src, double* dst, int nChar, double t) const {
    const std::vector<double> P = UniformizedTransition(Q_, k_, t);
    for (int c = 0; c < nChar; ++c) {
      const double* in = src + static_cast<std::size_t>(c) * k_;
      double* out = dst + static_cast<std::size_t>(c) * k_;
      for (int a = 0; a < k_; ++a) {
        double s = 0.0;
        for (int b = 0; b < k_; ++b) s += P[a * k_ + b] * in[b];
        out[a] = s;
      }
    }
  }

  // Q is symmetric, so P(t) is too
  void transition_transpose(const double* src, double* dst, int nChar,
                            double t) const {
    transition(src, dst, nChar, t);  // # nocov
  }

  double root_site_lik(const double* cl) const {
    double s = 0.0;
    for (int a = 0; a < k_; ++a) s += cl[a];
    return s / k_;
  }

 private:
  int k_;
  std::vector<double> Q_;  // row-major
};

struct EdgeList {
  std::vector<int> parent, child;
  std::vector<double> len;
};

// Edges below the root, 1-based, with their lengths
EdgeList Edges(const EvoTempo::TreeView& tr, const NumericVector& lengths) {
  EdgeList el;
  for (int v = 0; v < tr.nAll; ++v) {
    if (tr.parent[v] < 0) continue;
    el.parent.push_back(tr.parent[v] + 1);
    el.child.push_back(v + 1);
    el.len.push_back(lengths[v]);
  }
  return el;
}

// Site likelihoods (averaged over rate categories) of nCol columns given by
// their tip partials, tip-major: nCol blocks of K values per tip
template <class Model>
std::vector<double> SiteLikelihoods(const Model& model, const EdgeList& el, int nTip,
                                    int nCol, const double* tip,
                                    const std::vector<double>& rates) {
  const int K = model.n_states();
  const std::vector<int> states(static_cast<std::size_t>(nTip) * nCol, -1);
  MkPrime::Pruner<Model> pruner(model, nTip, nCol, states.data(), rates);
  pruner.set_tree(el.parent.data(), el.child.data(), el.len.data(),
                  static_cast<int>(el.len.size()));
  const std::size_t block = static_cast<std::size_t>(nCol) * K;
  for (int t = 1; t <= nTip; ++t) {
    pruner.set_tip_partials(t, tip + (t - 1) * block);
  }
  std::vector<double> site(nCol);
  pruner.site_likelihoods(site.data());
  return site;
}

}  // namespace

// Log-likelihood of the partitions of MkCharacters() (R/Characters.R) with
// rate multipliers `rates` (equiprobable categories); with `variable`, the
// Mkv correction per missing-data mask, from each partition's constant
// columns, which follow its nData data columns
// [[Rcpp::export(name = ".MkCharactersCpp")]]
double MkCharactersCpp(List tree, List partitions, NumericVector lengths,
                       NumericVector rates, bool variable) {
  const EvoTempo::TreeView tr(tree);
  const EdgeList el = Edges(tr, lengths);
  const std::vector<double> r(rates.begin(), rates.end());
  double total = 0.0;
  for (R_xlen_t i = 0; i < partitions.size(); ++i) {
    const List g = partitions[i];
    const int K = as<int>(g["nState"]);
    const bool ordered = as<bool>(g["ordered"]);
    const int nData = as<int>(g["nData"]);
    const int nCol = as<int>(g["nCol"]);
    const NumericVector tip = g["tip"];
    if (tip.size() != static_cast<R_xlen_t>(K) * nCol * tr.nTip) {
      stop("Partition tip partials do not match the tree");  // # nocov
    }
    std::vector<double> site;
    try {
      if (ordered) {
        const NumericMatrix Q = g["Q"];
        const UniformizedModel model(K, Q.begin());
        site = SiteLikelihoods(model, el, tr.nTip, nCol, tip.begin(), r);
      } else {
        const MkPrime::JCModel model(K);
        site = SiteLikelihoods(model, el, tr.nTip, nCol, tip.begin(), r);
      }
    } catch (const std::exception& e) {
      stop(e.what());  // # nocov
    }
    for (int c = 0; c < nData; ++c) {
      if (!(site[c] > 0.0)) return -std::numeric_limits<double>::infinity();
      total += std::log(site[c]);
    }
    if (variable) {
      const int nConst = as<int>(g["nConst"]);
      const IntegerVector count = g["count"];
      for (R_xlen_t m = 0; m < count.size(); ++m) {
        double p = 0.0;
        for (int s = 0; s < nConst; ++s) p += site[nData + m * nConst + s];
        if (!ordered) p *= K;
        // A constant pattern certain (to rounding): not a likelihood (R98-2)
        if (!(p < 1.0)) return -std::numeric_limits<double>::infinity();
        total -= count[m] * std::log1p(-p);
      }
    }
  }
  return total;
}

// characters: integer matrix, one row per tip; NA (or negative) if unknown
// [[Rcpp::export(name = ".MkLikelihoodCpp")]]
double MkLikelihoodCpp(List tree, IntegerMatrix characters, NumericVector lengths,
                       int nState) {
  const EvoTempo::TreeView tr(tree);
  if (characters.nrow() != tr.nTip) {
    stop("`characters` needs one row per tip");  // # nocov
  }
  const int nChar = characters.ncol();
  if (nChar == 0) return 0.0;
  std::vector<int> parent, child;
  std::vector<double> len;
  for (int v = 0; v < tr.nAll; ++v) {
    if (tr.parent[v] < 0) continue;
    parent.push_back(tr.parent[v] + 1);
    child.push_back(v + 1);
    len.push_back(lengths[v]);
  }
  std::vector<int> states(characters.begin(), characters.end());
  for (int& s : states) {
    if (s == NA_INTEGER || s < 0) s = -1;
  }
  try {
    const MkPrime::JCModel model(nState);
    return MkPrime::log_likelihood(model, parent.data(), child.data(), len.data(),
                                   static_cast<int>(len.size()), tr.nTip, nChar,
                                   states.data());
  } catch (const std::exception& e) {
    stop(e.what());
  }
}
