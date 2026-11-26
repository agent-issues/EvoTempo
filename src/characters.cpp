// Mk likelihood of discrete characters on a TempoTree, by the MkPrime
// header-only pruning kernels (LinkingTo: MkPrime).  MkLikelihood() and
// .MkCharactersLogLik() in R/Characters.R hold the pure-R references.
//
// For a plain character matrix (.MkCpp), partial likelihoods are cached between
// calls (dev/plans/mk-likelihood-cache.md): one MkPrime::Pruner per character
// matrix, which recomputes only the nodes whose subtree or child edge lengths
// changed since its last evaluation.  The cache never holds the chain state:
// every call passes the whole tree and every length, so a cached result is
// bitwise identical to a one-shot evaluation.  The partitions of MkCharacters()
// (.MkCharactersCpp) are pruned afresh at every call.

#include <Rcpp.h>
#include <MkPrime/MkPrime.h>
#include "tree.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
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

typedef MkPrime::Pruner<MkPrime::JCModel> MkPruner;

// One edge per non-root node, in node order; nodes 1-based
struct MkTree {
  std::vector<int> parent;
  std::vector<int> child;
  std::vector<double> len;
  bool operator==(const MkTree& o) const {
    return parent == o.parent && child == o.child && len == o.len;
  }
};

// Number of edges that differ; INT_MAX if the trees are not comparable
int NDiff(const MkTree& a, const MkTree& b) {
  if (a.child != b.child) return INT_MAX;
  int n = 0;
  for (std::size_t i = 0; i < a.len.size(); ++i) {
    n += a.parent[i] != b.parent[i] || a.len[i] != b.len[i];
  }
  return n;
}

struct MkEntry {
  int nTip;
  int nState;
  std::vector<int> states;
  std::unique_ptr<MkPruner> pruner;
  bool built = false;
  bool open = false;  // a checkpoint awaits accept() or reject()
  MkTree last, prev;
  std::vector<double> lastSite, prevSite;
  unsigned long long tick = 0;
};

const std::size_t kMaxEntries = 32;
std::vector<std::unique_ptr<MkEntry> > registry;
unsigned long long useCount = 0;
double nRecomputed = 0;
int nCheckFailures = 0;

std::unique_ptr<MkPruner> NewPruner(int nState, int nTip, int nChar,
                                    const std::vector<int>& states) {
  return std::unique_ptr<MkPruner>(
    new MkPruner(MkPrime::JCModel(nState), nTip, nChar, states.data()));
}

MkEntry& Lookup(int nTip, int nState, const std::vector<int>& states) {
  for (std::size_t i = 0; i < registry.size(); ++i) {
    MkEntry& e = *registry[i];
    if (e.nTip == nTip && e.nState == nState && e.states == states) {
      e.tick = ++useCount;
      return e;
    }
  }
  std::unique_ptr<MkEntry> e(new MkEntry());
  e->nTip = nTip;
  e->nState = nState;
  e->states = states;
  e->pruner = NewPruner(nState, nTip, static_cast<int>(states.size() / nTip), states);
  e->tick = ++useCount;
  if (registry.size() >= kMaxEntries) {
    std::size_t oldest = 0;
    for (std::size_t i = 1; i < registry.size(); ++i) {
      if (registry[i]->tick < registry[oldest]->tick) oldest = i;
    }
    registry.erase(registry.begin() + oldest);
  }
  registry.push_back(std::move(e));
  return *registry.back();
}

void Forget(const MkEntry* entry) {
  for (std::size_t i = 0; i < registry.size(); ++i) {
    if (registry[i].get() == entry) {
      registry.erase(registry.begin() + i);
      return;
    }
  }
}

void SetTree(MkPruner& pruner, const MkTree& x) {
  pruner.set_tree(x.parent.data(), x.child.data(), x.len.data(),
                  static_cast<int>(x.len.size()));
}

void Score(MkPruner& pruner, std::vector<double>& site) {
  site.resize(pruner.n_char());
  pruner.site_likelihoods(site.data());
  nRecomputed += pruner.n_recomputed();
}

// Brings `e` to tree `x`, inferring from `x` whether the open proposal was
// accepted (x differs less from it) or rejected (x differs less from the state
// before it); either guess gives the same result.
const std::vector<double>& Evaluate(MkEntry& e, const MkTree& x) {
  if (!e.built) {
    SetTree(*e.pruner, x);
    Score(*e.pruner, e.lastSite);
    e.last = x;
    e.built = true;
    return e.lastSite;
  }
  if (x == e.last) return e.lastSite;
  if (e.open) {
    e.open = false;
    if (x == e.prev) {
      e.pruner->reject();
      std::swap(e.last, e.prev);
      e.lastSite.swap(e.prevSite);
      return e.lastSite;
    }
    if (NDiff(x, e.prev) < NDiff(x, e.last)) {
      e.pruner->reject();
      std::swap(e.last, e.prev);
      e.lastSite.swap(e.prevSite);
    } else {
      e.pruner->accept();
    }
  }
  std::swap(e.prev, e.last);
  e.prevSite.swap(e.lastSite);
  e.pruner->checkpoint();
  e.open = true;
  SetTree(*e.pruner, x);
  Score(*e.pruner, e.lastSite);
  e.last = x;
  return e.lastSite;
}

// Sum of log site likelihoods, as MkPrime's Pruner::log_likelihood()
double LogSum(const std::vector<double>& site) {
  double logLik = 0.0;
  for (std::size_t c = 0; c < site.size(); ++c) {
    if (site[c] <= 0.0) return R_NegInf;
    logLik += std::log(site[c]);
  }
  return logLik;
}

bool SameBits(const std::vector<double>& a, const std::vector<double>& b) {
  return a.size() == b.size() &&
    std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

// Bitwise agreement of the cached partials (re-read, not recomputed) and the
// returned site likelihoods with a one-shot evaluation; failures are counted,
// so that a caller that catches the error cannot hide them
void Check(MkEntry& e, const MkTree& x, const std::vector<double>& site) {
  std::vector<double> reread(site.size()), fresh(site.size());
  e.pruner->site_likelihoods(reread.data());
  std::unique_ptr<MkPruner> oneShot = NewPruner(e.nState, e.nTip,
                                                e.pruner->n_char(), e.states);
  SetTree(*oneShot, x);
  oneShot->site_likelihoods(fresh.data());
  if (!SameBits(reread, fresh) || !SameBits(site, fresh)) {
    ++nCheckFailures;
    throw std::logic_error("cached Mk likelihood differs from a fresh evaluation");
  }
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

// characters: integer matrix, one row per tip; NA (or negative) if unknown.
// Returns the log-likelihood, or with `sites` one per character; `cache = false`
// evaluates without the cache.
// [[Rcpp::export(name = ".MkCpp")]]
NumericVector MkCpp(List tree, IntegerMatrix characters, NumericVector lengths,
                    int nState, bool sites, bool cache) {
  const int nTip = as<int>(tree["nTip"]);
  if (characters.nrow() != nTip) {
    stop("`characters` needs one row per tip");
  }
  const int nChar = characters.ncol();
  if (nChar == 0) return sites ? NumericVector(0) : NumericVector::create(0.0);
  if (nTip < 2) stop("Mk likelihood needs at least two tips");
  const IntegerVector par = tree["parent"];
  if (lengths.size() < par.size()) stop("`lengths` needs one value per node");
  MkTree x;
  for (int v = 0; v < par.size(); ++v) {
    if (par[v] == NA_INTEGER || par[v] < 1) continue;  // the root
    x.parent.push_back(par[v]);
    x.child.push_back(v + 1);
    x.len.push_back(lengths[v]);
  }
  std::vector<int> states(characters.begin(), characters.end());
  for (int& s : states) {
    if (s == NA_INTEGER || s < 0) s = -1;
  }
  MkEntry* entry = nullptr;
  std::vector<double> site;
  try {
    if (cache) {
      entry = &Lookup(nTip, nState, states);
      site = Evaluate(*entry, x);
      const char* check = std::getenv("EVOTEMPO_MK_CACHE_CHECK");
      if (check && *check && std::string(check) != "0") Check(*entry, x, site);
    } else {
      std::unique_ptr<MkPruner> pruner = NewPruner(nState, nTip, nChar, states);
      SetTree(*pruner, x);
      site.resize(nChar);
      pruner->site_likelihoods(site.data());
    }
  } catch (const std::exception& e) {
    // A failed evaluation may leave the entry part-way through an update
    if (entry) Forget(entry);
    stop("%s", e.what());
  }
  if (!sites) return NumericVector::create(LogSum(site));
  NumericVector out(nChar);
  for (int c = 0; c < nChar; ++c) out[c] = std::log(site[c]);
  return out;
}

// `corrupt` (for tests of the self-check) alters every entry's stored result
// [[Rcpp::export(name = ".MkCacheCpp")]]
List MkCacheCpp(bool clear, bool corrupt) {
  if (clear) {
    registry.clear();
    nRecomputed = 0;
    nCheckFailures = 0;
  }
  if (corrupt) {
    for (std::size_t i = 0; i < registry.size(); ++i) {
      for (double& v : registry[i]->lastSite) v *= 0.5;
    }
  }
  return List::create(_["entries"] = static_cast<int>(registry.size()),
                      _["recomputed"] = nRecomputed,
                      _["checkFailures"] = nCheckFailures);
}
