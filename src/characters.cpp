// Mk likelihood of discrete characters on a TempoTree, by the MkPrime
// header-only pruning kernels (LinkingTo: MkPrime).  MkLikelihood() in
// R/Characters.R holds the pure-R reference.

#include <Rcpp.h>
#include <MkPrime/MkPrime.h>
#include "tree.h"
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
// .Call()s, behind an external pointer.  Each evaluation hands the pruner the
// tree's parent vector and every branch length: with the topology unchanged,
// set_edge_length() marks only the changed branches and their ancestors
// dirty; after a topology change, set_tree() diffs the children of each node.
// Either way the pruner recomputes the dirty nodes only, in the same order and
// from the same inputs as a fresh pruning, so the value is bitwise that of
// .MkLikelihoodCpp() whatever the history of the cache.  The first
// evaluation after a resolution opens a checkpoint(); Resolve() then accepts
// or rejects it, a rejection restoring the journaled CLs without recomputing.
// Correctness never depends on Resolve(): an unresolved proposal is simply
// diffed against at the next evaluation.  All state is owned C++ storage; no
// R vector is held across calls.
namespace {

class CharacterCache {
 public:
  CharacterCache(int nTip, int nChar, const std::vector<int>& states, int nState)
    : pruner_(MkPrime::JCModel(nState), nTip, nChar, states.data()) {}

  double LogLikelihood(const IntegerVector& parent, const NumericVector& lengths) {
    const int nAll = parent.size();
    if (lengths.size() != nAll) {
      stop("`lengths` needs one entry per node");  // # nocov
    }
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
  // Edges: every non-root node, in increasing order (as .MkLikelihoodCpp())
  void Build(const IntegerVector& parent, const NumericVector& lengths) {
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

// [[Rcpp::export(name = ".CachedMkLikelihoodCpp", rng = false)]]
double CachedMkLikelihoodCpp(SEXP cache, IntegerVector parent, NumericVector lengths) {
  CharacterCache* cc = CacheFrom(cache);
  try {
    return cc->LogLikelihood(parent, lengths);
  } catch (const std::exception& e) {
    stop(e.what());  // # nocov
  }
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
