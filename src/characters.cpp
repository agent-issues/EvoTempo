// Mk likelihood of discrete characters on a TempoTree, by the MkPrime
// header-only pruning kernels (LinkingTo: MkPrime).  MkLikelihood() in
// R/Characters.R holds the pure-R reference.
//
// Partial likelihoods are cached between calls (dev/plans/mk-likelihood-cache.md):
// one MkPrime::Pruner per character matrix, which recomputes only the nodes whose
// subtree or child edge lengths changed since its last evaluation.  The cache
// never holds the chain state: every call passes the whole tree and every length,
// so a cached result is bitwise identical to a one-shot evaluation.

#include <Rcpp.h>
#include <MkPrime/MkPrime.h>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace {

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
  unsigned long tick = 0;
};

const std::size_t kMaxEntries = 32;
std::vector<std::unique_ptr<MkEntry> > registry;
unsigned long useCount = 0;
double nRecomputed = 0;

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
  return std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

bool CheckRequested() {
  const char* check = std::getenv("EVOTEMPO_MK_CACHE_CHECK");
  return check && *check && std::string(check) != "0";
}

// Exact agreement of the cached partials (re-read, not recomputed) and the
// returned site likelihoods with a one-shot evaluation of `states`, the tip
// states that the caller holds
void Check(MkSet& s, const std::vector<std::vector<int> >& states) {
  if (states != s.states) {
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
    // # nocov start
    if (!SameBits(reread, fresh) || !SameBits(site, fresh)) {
      throw std::logic_error("cached Mk likelihood differs from a fresh evaluation");
    }
    // # nocov end
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
  x.parent.clear();
  x.child.clear();
  x.len.clear();
  for (int v = 0; v < par.size(); ++v) {
    if (par[v] == NA_INTEGER) continue;
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

// characters: integer matrix, one row per tip; NA (or negative) if unknown.
// Returns the log-likelihood, or with `sites` one per character; `cache = false`
// evaluates without the cache.
// [[Rcpp::export(name = ".MkCpp")]]
NumericVector MkCpp(List tree, IntegerMatrix characters, NumericVector lengths,
                    int nState, bool sites, bool cache) {
  const int nTip = as<int>(tree["nTip"]);
  if (characters.nrow() != nTip) {
    stop("`characters` needs one row per tip");  // # nocov
  }
  const int nChar = characters.ncol();
  if (nChar == 0) return sites ? NumericVector(0) : NumericVector::create(0.0);
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
      if (CheckRequested()) Check(entry->set, entry->set.states);
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
    stop(e.what());
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
  std::unique_ptr<MkSet> set(new MkSet());
  set->nTip = nTip;
  try {
    for (R_xlen_t i = 0; i < matrices.size(); ++i) {
      const IntegerMatrix m = matrices[i];
      if (nTip < 2 || m.nrow() != nTip) stop("each matrix needs one row per tip");
      set->Add(nState[i], TipStates(m), sites[i]);
    }
  } catch (const std::invalid_argument& e) {
    stop(e.what());
  }
  XPtr<MkSet> handle(set.release(), true, SetTag());
  return handle;
}

// For each matrix of the set, in order: the log-likelihood, or where registered
// with `sites` one value per character.  One diff of the tree serves every
// matrix.  NULL if `handle` is not a live handle (they do not survive
// serialization).  `matrices` is read only by the EVOTEMPO_MK_CACHE_CHECK check.
// [[Rcpp::export(name = ".MkSetLogLikCpp")]]
SEXP MkSetLogLikCpp(SEXP handle, IntegerVector parent, NumericVector lengths,
                    List matrices) {
  if (TYPEOF(handle) != EXTPTRSXP || R_ExternalPtrTag(handle) != SetTag() ||
      !R_ExternalPtrAddr(handle)) {
    return R_NilValue;
  }
  MkSet& s = *static_cast<MkSet*>(R_ExternalPtrAddr(handle));
  if (lengths.size() != parent.size()) stop("`lengths` needs one entry per node");
  ReadTree(parent, lengths, s.next);
  try {
    Evaluate(s, s.next);
    if (CheckRequested()) {
      std::vector<std::vector<int> > states;
      for (R_xlen_t i = 0; i < matrices.size(); ++i) {
        states.push_back(TipStates(as<IntegerMatrix>(matrices[i])));
      }
      Check(s, states);
    }
  } catch (const std::exception& e) {
    s.Reset();
    stop(e.what());
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

// [[Rcpp::export(name = ".MkCacheCpp")]]
List MkCacheCpp(bool clear) {
  if (clear) {
    registry.clear();
    nRecomputed = 0;
  }
  return List::create(_["entries"] = static_cast<int>(registry.size()),
                      _["recomputed"] = nRecomputed);
}
