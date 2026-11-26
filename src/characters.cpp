// Mk likelihood of discrete characters on a TempoTree, by the MkPrime
// header-only pruning kernels (LinkingTo: MkPrime).  MkLikelihood() and
// .MkCharactersLogLik() in R/Characters.R hold the pure-R references.
//
// For a plain character matrix (.MkCpp), partial likelihoods are cached between
// calls (dev/plans/mk-likelihood-cache.md): one MkPrime::Pruner per character
// matrix, which recomputes only the nodes whose subtree or child edge lengths
// changed since its last evaluation.  The cache never holds the chain state:
// every call passes the whole tree and every length, so a cached result is
// bitwise identical to a one-shot evaluation.  Matrices that are always scored
// together share one diff of the tree through a handle (.MkSetCpp,
// .MkSetLogLikCpp).  The partitions of MkCharacters() (.MkCharactersCpp) are
// pruned afresh at every call.

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

// Pruners of a fixed list of character matrices on one tree, evaluated together
struct MkSet {
  int nTip = 0;
  std::vector<int> nState;
  std::vector<char> sites;  // per matrix: one value per character, not their sum
  std::vector<std::vector<int> > states;
  std::vector<std::unique_ptr<MkPruner> > pruner;  // null for no characters
  std::vector<std::vector<double> > lastSite, prevSite;
  bool built = false;
  bool open = false;  // a checkpoint awaits accept() or reject()
  MkTree last, prev, next;
  std::size_t nOut = 0;
  unsigned long long generation = 0;  // of the cache, when last evaluated

  int NChar(std::size_t i) const {
    return static_cast<int>(states[i].size() / nTip);
  }
  std::unique_ptr<MkPruner> NewPruner(std::size_t i,
                                      const std::vector<int>& tipStates) const {
    return std::unique_ptr<MkPruner>(
      new MkPruner(MkPrime::JCModel(nState[i]), nTip, NChar(i), tipStates.data()));
  }
  void Add(int k, const std::vector<int>& tipStates, bool bySite) {
    nState.push_back(k);
    sites.push_back(bySite);
    states.push_back(tipStates);
    const std::size_t i = states.size() - 1;
    pruner.push_back(NChar(i) ? NewPruner(i, tipStates) : nullptr);
    lastSite.push_back(std::vector<double>());
    prevSite.push_back(std::vector<double>());
    nOut += bySite ? NChar(i) : 1;
  }
  // A failed evaluation may leave a pruner part-way through an update
  void Reset() {
    for (std::size_t i = 0; i < pruner.size(); ++i) {
      if (pruner[i]) pruner[i] = NewPruner(i, states[i]);
    }
    built = false;
    open = false;
  }
};

struct MkEntry {
  MkSet set;  // one matrix
  unsigned long long tick = 0;
};

const std::size_t kMaxEntries = 32;
std::vector<std::unique_ptr<MkEntry> > registry;
unsigned long long useCount = 0;
double nRecomputed = 0;
int nCheckFailures = 0;
unsigned long long generation = 0;  // advanced when the cache is cleared

MkEntry& Lookup(int nTip, int nState, const std::vector<int>& states) {
  for (std::size_t i = 0; i < registry.size(); ++i) {
    MkEntry& e = *registry[i];
    if (e.set.nTip == nTip && e.set.nState[0] == nState && e.set.states[0] == states) {
      e.tick = ++useCount;
      return e;
    }
  }
  std::unique_ptr<MkEntry> e(new MkEntry());
  e->set.nTip = nTip;
  e->set.Add(nState, states, true);
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

// Sets every pruner to tree `x` and reads the site likelihoods into `lastSite`
void Score(MkSet& s, const MkTree& x) {
  for (std::size_t i = 0; i < s.pruner.size(); ++i) {
    if (!s.pruner[i]) continue;
    MkPruner& pruner = *s.pruner[i];
    SetTree(pruner, x);
    s.lastSite[i].resize(pruner.n_char());
    pruner.site_likelihoods(s.lastSite[i].data());
    nRecomputed += pruner.n_recomputed();
  }
  s.last = x;
}

// Brings `s` to tree `x`, inferring from `x` whether the open proposal was
// accepted (x differs less from it) or rejected (x differs less from the state
// before it); either guess gives the same result.
void Evaluate(MkSet& s, const MkTree& x) {
  if (!s.built) {
    Score(s, x);
    s.built = true;
    return;
  }
  if (x == s.last) return;
  if (s.open) {
    s.open = false;
    const bool back = x == s.prev;
    const bool rejected = back || NDiff(x, s.prev) < NDiff(x, s.last);
    for (std::size_t i = 0; i < s.pruner.size(); ++i) {
      if (!s.pruner[i]) continue;
      if (rejected) {
        s.pruner[i]->reject();
      } else {
        s.pruner[i]->accept();
      }
    }
    if (rejected) {
      std::swap(s.last, s.prev);
      s.lastSite.swap(s.prevSite);
      if (back) return;
    }
  }
  std::swap(s.prev, s.last);
  s.prevSite.swap(s.lastSite);
  for (std::size_t i = 0; i < s.pruner.size(); ++i) {
    if (s.pruner[i]) s.pruner[i]->checkpoint();
  }
  s.open = true;
  Score(s, x);
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

bool CheckRequested() {
  const char* check = std::getenv("EVOTEMPO_MK_CACHE_CHECK");
  return check && *check && std::string(check) != "0";
}

// Bitwise agreement of the cached partials (re-read, not recomputed) and the
// returned site likelihoods with a one-shot evaluation of `states` and `nState`,
// the tip states and numbers of states that the caller holds; failures are
// counted, so that a caller that catches the error cannot hide them
void Check(MkSet& s, const std::vector<std::vector<int> >& states,
           const std::vector<int>& nState) {
  if (nState != s.nState || states != s.states) {
    ++nCheckFailures;
    throw std::logic_error("character matrices differ from those registered");
  }
  for (std::size_t i = 0; i < s.pruner.size(); ++i) {
    if (!s.pruner[i]) continue;
    const std::vector<double>& site = s.lastSite[i];
    std::vector<double> reread(site.size()), fresh(site.size());
    s.pruner[i]->site_likelihoods(reread.data());
    std::unique_ptr<MkPruner> oneShot = s.NewPruner(i, states[i]);
    SetTree(*oneShot, s.last);
    oneShot->site_likelihoods(fresh.data());
    if (!SameBits(reread, fresh) || !SameBits(site, fresh)) {
      ++nCheckFailures;
      throw std::logic_error("cached Mk likelihood differs from a fresh evaluation");
    }
  }
}

std::vector<int> TipStates(const IntegerMatrix& characters) {
  std::vector<int> states(characters.begin(), characters.end());
  for (int& s : states) {
    if (s == NA_INTEGER || s < 0) s = -1;
  }
  return states;
}

// One edge per non-root node
void ReadTree(const IntegerVector& par, const NumericVector& lengths, MkTree& x) {
  if (lengths.size() < par.size()) stop("`lengths` needs one value per node");
  x.parent.clear();
  x.child.clear();
  x.len.clear();
  for (int v = 0; v < par.size(); ++v) {
    if (par[v] == NA_INTEGER || par[v] < 1) continue;  // the root
    x.parent.push_back(par[v]);
    x.child.push_back(v + 1);
    x.len.push_back(lengths[v]);
  }
}

SEXP SetTag() {
  static SEXP tag = Rf_install("EvoTempoMkSet");
  return tag;
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
  MkTree x;
  ReadTree(tree["parent"], lengths, x);
  const std::vector<int> states = TipStates(characters);
  MkEntry* entry = nullptr;
  std::vector<double> site;
  try {
    if (cache) {
      entry = &Lookup(nTip, nState, states);
      Evaluate(entry->set, x);
      site = entry->set.lastSite[0];
      if (CheckRequested()) Check(entry->set, entry->set.states, entry->set.nState);
    } else {
      MkSet oneShot;
      oneShot.nTip = nTip;
      oneShot.Add(nState, states, true);
      SetTree(*oneShot.pruner[0], x);
      site.resize(nChar);
      oneShot.pruner[0]->site_likelihoods(site.data());
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

// Registers character matrices (as `characters` of .MkCpp(), each with its own
// number of states) whose likelihoods are always wanted together, and returns
// a handle for .MkSetLogLikCpp()
// [[Rcpp::export(name = ".MkSetCpp")]]
SEXP MkSetCpp(int nTip, List matrices, IntegerVector nState, LogicalVector sites) {
  if (nState.size() != matrices.size() || sites.size() != matrices.size()) {
    stop("`nState` and `sites` need one entry per matrix");
  }
  if (nTip < 2) stop("Mk likelihood needs at least two tips");
  if (is_true(any(is_na(sites)))) stop("`sites` must be TRUE or FALSE");
  std::unique_ptr<MkSet> set(new MkSet());
  set->nTip = nTip;
  set->generation = generation;
  try {
    for (R_xlen_t i = 0; i < matrices.size(); ++i) {
      const IntegerMatrix m = matrices[i];
      if (m.nrow() != nTip) stop("each matrix needs one row per tip");
      set->Add(nState[i], TipStates(m), sites[i]);
    }
  } catch (const std::invalid_argument& e) {
    stop("%s", e.what());
  }
  XPtr<MkSet> handle(set.release(), true, SetTag());
  return handle;
}

// For each matrix of the set, in order: the log-likelihood, or where registered
// with `sites` one value per character.  One diff of the tree serves every
// matrix.  NULL if `handle` is not a live handle (they do not survive
// serialization).  `matrices` and `nState`, the data that the caller means, are
// read only by the EVOTEMPO_MK_CACHE_CHECK check.
// [[Rcpp::export(name = ".MkSetLogLikCpp")]]
SEXP MkSetLogLikCpp(SEXP handle, IntegerVector parent, NumericVector lengths,
                    List matrices, IntegerVector nState) {
  if (TYPEOF(handle) != EXTPTRSXP || R_ExternalPtrTag(handle) != SetTag() ||
      !R_ExternalPtrAddr(handle)) {
    return R_NilValue;
  }
  MkSet& s = *static_cast<MkSet*>(R_ExternalPtrAddr(handle));
  if (s.generation != generation) {
    s.Reset();
    s.generation = generation;
  }
  ReadTree(parent, lengths, s.next);
  try {
    Evaluate(s, s.next);
    if (CheckRequested()) {
      std::vector<std::vector<int> > states;
      for (R_xlen_t i = 0; i < matrices.size(); ++i) {
        states.push_back(TipStates(as<IntegerMatrix>(matrices[i])));
      }
      Check(s, states, std::vector<int>(nState.begin(), nState.end()));
    }
  } catch (const std::exception& e) {
    s.Reset();
    stop("%s", e.what());
  }
  NumericVector out(s.nOut);
  R_xlen_t at = 0;
  for (std::size_t i = 0; i < s.pruner.size(); ++i) {
    const std::vector<double>& site = s.lastSite[i];
    if (s.sites[i]) {
      for (std::size_t c = 0; c < site.size(); ++c) out[at++] = std::log(site[c]);
    } else {
      out[at++] = LogSum(site);
    }
  }
  return out;
}

// `corrupt` (for tests of the self-check) alters every entry's stored result
// [[Rcpp::export(name = ".MkCacheCpp")]]
List MkCacheCpp(bool clear, bool corrupt) {
  if (clear) {
    registry.clear();
    ++generation;
    nRecomputed = 0;
    nCheckFailures = 0;
  }
  if (corrupt) {
    for (std::size_t i = 0; i < registry.size(); ++i) {
      for (double& v : registry[i]->set.lastSite[0]) v *= 0.5;
    }
  }
  return List::create(_["entries"] = static_cast<int>(registry.size()),
                      _["recomputed"] = nRecomputed,
                      _["checkFailures"] = nCheckFailures);
}
