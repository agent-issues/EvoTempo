// Mk likelihood of discrete characters on a TempoTree, by the MkPrime
// header-only pruning kernels (LinkingTo: MkPrime).  MkLikelihood() in
// R/Characters.R holds the pure-R reference.

#include <Rcpp.h>
#include <MkPrime/MkPrime.h>
#include "paths.h"
#include "tree.h"
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <vector>

using namespace Rcpp;

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

// Persistent character likelihood for TempoMCMC() (#26, item 3; design in
// dev/plans/incremental-pruning.md).  A CharacterCache owns a
// MkPrime::Pruner<JCModel> whose conditional likelihoods persist across
// .Call()s, behind an external pointer, and a memo of each branch's
// tempo-weighted length.
//
// Lengths: for each branch the memo keeps a copy of its path (bins and jump
// times), its end ages and sum_i e^(gamma x_k_i) dwell_i.  A branch whose path
// and ages equal the memo's reuses that sum; others are summed afresh exactly
// as .CharacterLengthsCpp() does, so every length is bitwise the same.  A new
// grid x or gamma clears the memo; `rate` multiplies the sums.
//
// Pruning: each evaluation hands the pruner the parent vector and every
// length.  With the topology unchanged, set_edge_length() marks only the
// changed branches and their ancestors dirty; after a topology change,
// set_tree() diffs the children of each node.  Either way the pruner
// recomputes the dirty nodes only, in the same order and from the same inputs
// as a fresh pruning, so the value is bitwise that of .MkLikelihoodCpp()
// whatever the history of the cache.  The first evaluation after a resolution
// opens a checkpoint(); Resolve() accepts or rejects it, a rejection restoring
// the journaled CLs without recomputation.  Correctness never depends on
// Resolve(): an unresolved proposal is simply diffed against at the next
// evaluation.  All state is owned C++ storage; no R object is held across
// calls.
namespace {

// Element `name` of the list `x`, or R_NilValue
SEXP ListElement(SEXP x, const char* name) {
  const SEXP names = Rf_getAttrib(x, R_NamesSymbol);
  if (names == R_NilValue) return R_NilValue;  // # nocov
  for (R_xlen_t i = 0; i < Rf_xlength(x); ++i) {
    if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) return VECTOR_ELT(x, i);
  }
  return R_NilValue;  // # nocov
}

class CharacterCache {
 public:
  CharacterCache(int nTip, int nChar, const std::vector<int>& states, int nState)
    : pruner_(MkPrime::JCModel(nState), nTip, nChar, states.data()) {}

  double LogLikelihood(const IntegerVector& parent, const double* lengths) {
    if (!built_) {
      Build(parent, lengths);
      built_ = true;
    }
    if (!pending_) {
      pruner_.checkpoint();
      pending_ = true;
    }
    if (SameTopology(parent)) {
      for (std::size_t e = 0; e < child_.size(); ++e) {
        pruner_.set_edge_length(static_cast<int>(e), lengths[child_[e] - 1]);
      }
    } else {
      Build(parent, lengths);
    }
    return pruner_.log_likelihood();
  }

  // Tempo-weighted lengths, as .CharacterLengthsCpp(); fills lengths_
  const std::vector<double>& Lengths(List tree, const NumericVector& x, List paths,
                                     double rate, double gamma) {
    const int nTip = as<int>(tree["nTip"]);
    const int nAll = nTip + as<int>(tree["nNode"]);
    const int root = as<int>(tree["root"]) - 1;
    const double stem = as<double>(tree["stem"]);
    const bool origin = tree.containsElementNamed("origin") && as<bool>(tree["origin"]);
    const IntegerVector parent = tree["parent"];
    const NumericVector age = tree["age"];
    const CharacterVector tipType = tree["tipType"];
    if (gamma != gamma_ || static_cast<int>(x_.size()) != x.size() ||
        !std::equal(x_.begin(), x_.end(), x.begin()) ||
        static_cast<int>(memo_.size()) != nAll) {
      gamma_ = gamma;
      x_.assign(x.begin(), x.end());
      memo_.assign(nAll, Memo());
    }
    lengths_.assign(nAll, 0.0);
    for (int v = 0; v < nAll; ++v) {
      if (parent[v] == NA_INTEGER && !(origin && v == root)) continue;
      if (v < nTip && std::strcmp(CHAR(STRING_ELT(tipType, v)), "ancestor") == 0) continue;
      const double top = v == root ? age[root] + stem : age[parent[v] - 1];
      const double bottom = age[v];
      const SEXP path = VECTOR_ELT(paths, v);
      const IntegerVector k(ListElement(path, "k"));
      const NumericVector t(ListElement(path, "t"));
      Memo& m = memo_[v];
      if (!(m.valid && m.top == top && m.bottom == bottom &&
            static_cast<int>(m.k.size()) == k.size() &&
            static_cast<int>(m.t.size()) == t.size() &&
            std::equal(m.k.begin(), m.k.end(), k.begin()) &&
            std::equal(m.t.begin(), m.t.end(), t.begin()))) {
        const int n = k.size();
        terms_.assign(n, 0.0);
        for (int i = 0; i < n; ++i) {
          const double hi = i == 0 ? top : t[i - 1];
          const double lo = i == n - 1 ? bottom : t[i];
          terms_[i] = std::exp(gamma * x[k[i] - 1]) * (hi - lo);
        }
        m.sum = EvoTempo::LongSum(terms_);
        m.k.assign(k.begin(), k.end());
        m.t.assign(t.begin(), t.end());
        m.top = top;
        m.bottom = bottom;
        m.valid = true;
      }
      lengths_[v] = rate * m.sum;
    }
    return lengths_;
  }

  void Resolve(bool accept) {
    if (!pending_) return;
    if (accept) {
      pruner_.accept();
    } else {
      pruner_.reject();
    }
    pending_ = false;
  }

  int Recomputed() const { return pruner_.n_recomputed(); }

 private:
  struct Memo {
    bool valid = false;
    double top = 0.0, bottom = 0.0, sum = 0.0;
    std::vector<int> k;
    std::vector<double> t;
  };

  // Edges: every non-root node, in increasing order (as .MkLikelihoodCpp())
  void Build(const IntegerVector& parent, const double* lengths) {
    const int nAll = parent.size();
    parentOf_.clear();
    child_.clear();
    len_.clear();
    for (int v = 0; v < nAll; ++v) {
      if (parent[v] == NA_INTEGER) continue;
      parentOf_.push_back(parent[v]);
      child_.push_back(v + 1);
      len_.push_back(lengths[v]);
    }
    pruner_.set_tree(parentOf_.data(), child_.data(), len_.data(),
                     static_cast<int>(len_.size()));
  }

  // The pruner's own topology is the reference, so that reject() needs no
  // book-keeping here
  bool SameTopology(const IntegerVector& parent) const {
    const int nAll = parent.size();
    if (nAll != pruner_.n_node()) return false;  // # nocov
    for (int v = 0; v < nAll; ++v) {
      const int p = parent[v] == NA_INTEGER ? 0 : parent[v];
      if (pruner_.parent_of(v + 1) != p) return false;
    }
    return true;
  }

  MkPrime::Pruner<MkPrime::JCModel> pruner_;
  bool built_ = false;
  bool pending_ = false;
  std::vector<int> parentOf_, child_;
  std::vector<double> len_;
  double gamma_ = std::nan("");
  std::vector<double> x_, lengths_, terms_;
  std::vector<Memo> memo_;
};

CharacterCache* CacheFrom(SEXP cache) {
  XPtr<CharacterCache> ptr(cache);
  if (!ptr.get()) {
    stop("Character cache is no longer valid (was it saved and reloaded?)");
  }
  return ptr.get();
}

}  // namespace

// characters: integer matrix, one row per tip, at least one column
// [[Rcpp::export(name = ".CharacterCacheCpp", rng = false)]]
SEXP CharacterCacheCpp(IntegerMatrix characters, int nState) {
  std::vector<int> states(characters.begin(), characters.end());
  for (int& s : states) {
    if (s == NA_INTEGER || s < 0) s = -1;
  }
  try {
    XPtr<CharacterCache> ptr(
      new CharacterCache(characters.nrow(), characters.ncol(), states, nState), true);
    return ptr;
  } catch (const std::exception& e) {
    stop(e.what());
  }
}

// Log-likelihood on given lengths (one per node)
// [[Rcpp::export(name = ".CachedMkLikelihoodCpp", rng = false)]]
double CachedMkLikelihoodCpp(SEXP cache, IntegerVector parent, NumericVector lengths) {
  CharacterCache* cc = CacheFrom(cache);
  if (lengths.size() != parent.size()) {
    stop("`lengths` needs one entry per node");
  }
  try {
    return cc->LogLikelihood(parent, lengths.begin());
  } catch (const std::exception& e) {
    stop(e.what());  // # nocov
  }
}

// Log-likelihood on the tempo-weighted lengths of `paths` (.CharacterLengths())
// [[Rcpp::export(name = ".CachedCharacterLikelihoodCpp", rng = false)]]
double CachedCharacterLikelihoodCpp(SEXP cache, List tree, NumericVector x, List paths,
                                    double rate, double gamma) {
  CharacterCache* cc = CacheFrom(cache);
  const std::vector<double>& len = cc->Lengths(tree, x, paths, rate, gamma);
  const IntegerVector parent = tree["parent"];
  try {
    return cc->LogLikelihood(parent, len.data());
  } catch (const std::exception& e) {
    stop(e.what());  // # nocov
  }
}

// The memoized lengths, for tests
// [[Rcpp::export(name = ".CachedCharacterLengthsCpp", rng = false)]]
NumericVector CachedCharacterLengthsCpp(SEXP cache, List tree, NumericVector x,
                                        List paths, double rate, double gamma) {
  const std::vector<double>& len = CacheFrom(cache)->Lengths(tree, x, paths, rate, gamma);
  return NumericVector(len.begin(), len.end());
}

// [[Rcpp::export(name = ".ResolveCharacterCacheCpp", rng = false)]]
void ResolveCharacterCacheCpp(SEXP cache, bool accept) {
  CacheFrom(cache)->Resolve(accept);
}

// Number of nodes recomputed by the last evaluation
// [[Rcpp::export(name = ".CharacterCacheRecomputedCpp", rng = false)]]
int CharacterCacheRecomputedCpp(SEXP cache) {
  return CacheFrom(cache)->Recomputed();
}
