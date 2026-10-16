#ifndef EVOTEMPO_TREE_H
#define EVOTEMPO_TREE_H

// Zero-based view of a TempoTree (see TempoTree() in R/TreeLikelihood.R)

#include <Rcpp.h>
#include <string>
#include <vector>

namespace EvoTempo {

enum TipType { kExtant = 0, kFossil = 1, kAncestor = 2 };

struct TreeView {
  int nTip;
  int nNode;
  int nAll;
  int root;                 // zero-based
  double stem;
  std::vector<int> parent;  // -1 at the root
  std::vector<int> kids;    // 2 per internal node, row-major, zero-based
  std::vector<double> age;
  std::vector<int> tipType;
  std::vector<int> sampled;   // per internal node
  std::vector<int> postorder; // zero-based internal nodes

  explicit TreeView(Rcpp::List tree) {
    nTip = Rcpp::as<int>(tree["nTip"]);
    nNode = Rcpp::as<int>(tree["nNode"]);
    nAll = nTip + nNode;
    root = Rcpp::as<int>(tree["root"]) - 1;
    stem = Rcpp::as<double>(tree["stem"]);
    const Rcpp::IntegerVector par = tree["parent"];
    parent.resize(nAll);
    for (int i = 0; i < nAll; ++i) {
      parent[i] = par[i] == NA_INTEGER ? -1 : par[i] - 1;
    }
    const Rcpp::IntegerMatrix k = tree["kids"];
    kids.resize(2 * nNode);
    for (int i = 0; i < nNode; ++i) {
      kids[2 * i] = k(i, 0) - 1;
      kids[2 * i + 1] = k(i, 1) - 1;
    }
    const Rcpp::NumericVector a = tree["age"];
    age.assign(a.begin(), a.end());
    const Rcpp::CharacterVector tt = tree["tipType"];
    tipType.resize(nTip);
    for (int i = 0; i < nTip; ++i) {
      const std::string s = Rcpp::as<std::string>(tt[i]);
      tipType[i] = s == "extant" ? kExtant : s == "fossil" ? kFossil : kAncestor;
    }
    const Rcpp::LogicalVector sa = tree["sampled"];
    sampled.assign(sa.begin(), sa.end());
    const Rcpp::IntegerVector po = tree["postorder"];
    postorder.resize(po.size());
    for (int i = 0; i < po.size(); ++i) postorder[i] = po[i] - 1;
  }

  bool IsAncestorTip(int v) const { return v < nTip && tipType[v] == kAncestor; }

  // Branches carrying paths: every non-root node except sampled-ancestor tips
  // (.PathNodes), in increasing order
  std::vector<int> PathNodes() const {
    std::vector<int> out;
    for (int v = 0; v < nAll; ++v) {
      if (parent[v] >= 0 && !IsAncestorTip(v)) out.push_back(v);
    }
    return out;
  }

  // .Compatible(tree, condition = "extant", origin = FALSE)
  bool CrownCompatible() const {
    std::vector<int> covered(nAll, 0);
    for (int v = 0; v < nTip; ++v) {
      if (tipType[v] != kExtant) continue;
      int u = v;
      while (parent[u] != root) u = parent[u];
      covered[u] = 1;
    }
    const int r = root - nTip;
    return covered[kids[2 * r]] && covered[kids[2 * r + 1]];
  }
};

}  // namespace EvoTempo

#endif  // EVOTEMPO_TREE_H
