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
#include <cstring>
#include <cmath>
#include <cstdlib>
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
