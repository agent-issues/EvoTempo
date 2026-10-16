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
  if (characters.nrow() != tr.nTip) stop("`characters` needs one row per tip");
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
