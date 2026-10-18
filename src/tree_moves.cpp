// Time-tree moves: C++ engine for the kernels of R/TreeMoves.R.  Each function
// reproduces its namesake there (without the "Cpp" suffix) and draws random
// numbers from R's generator in the same order.  The composite kernels
// (.DrawBridgesCpp, .LogBridgesCpp, .NodeLogWeightsCpp, .DrawTipCpp,
// .LogTipQCpp) read the frozen model once per call.

#include <Rcpp.h>
#include "frozen.h"
#include "paths.h"
#include "tree.h"
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace EvoTempo;

namespace {

// max() as R computes it: NaN if any element is NaN
double RMax(const std::vector<double>& x) {
  double top = R_NegInf;
  for (double v : x) {
    if (std::isnan(v)) return v;
    if (v > top) top = v;
  }
  return top;
}

double LogSumExp(const std::vector<double>& x) {
  const double top = RMax(x);
  std::vector<double> e(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) e[i] = std::exp(x[i] - top);
  return top + std::log(LongSum(e));
}

// .LogTransferRow(pc, k) (row = true) or .LogTransferCol(pc, k); k one-based
std::vector<double> LogTransfer(const Pieces& pc, int k, bool row, int K) {
  if (k < 1 || k > K) Rcpp::stop("Bin %d out of range", k);
  std::vector<double> r(K, 0.0), x(K);
  r[k - 1] = 1.0;
  double logScale = 0.0;
  const int m = pc.M.size();
  for (int step = 0; step < m; ++step) {
    const std::vector<double>& M = pc.M[row ? step : m - 1 - step];
    if (row) {
      for (int j = 0; j < K; ++j) {
        double s = 0.0;
        for (int i = 0; i < K; ++i) s += r[i] * M[i + static_cast<std::size_t>(j) * K];
        x[j] = s;
      }
    } else {
      MatVec(M, r, x, K);
    }
    const double top = RMax(x);
    for (int i = 0; i < K; ++i) r[i] = x[i] / top;
    logScale += std::log(top);
  }
  for (int i = 0; i < K; ++i) r[i] = std::log(std::max(r[i], 0.0)) + logScale;
  return r;
}

int PiecesK(const Pieces& pc) {
  return static_cast<int>(std::lround(std::sqrt(static_cast<double>(pc.M[0].size()))));
}

double LogPathExact(const TreeView& tr, const FrozenSet& fs, int v, const Path& p) {
  const int n = p.k.size();
  std::vector<double> terms(n);
  for (int i = 0; i < n; ++i) {
    const double hi = i == 0 ? tr.age[tr.parent[v]] : p.t[i - 1];
    const double lo = i == n - 1 ? tr.age[v] : p.t[i];
    const int k = p.k[i] - 1;
    terms[i] = (fs.q(k, k) - fs.sigma[k]) * (hi - lo) +
      2 * fs.lambda[k] * (fs.Hermite(hi, k, true) - fs.Hermite(lo, k, true));
  }
  long double jumps = 0.0;
  for (int i = 0; i + 1 < n; ++i) jumps += std::log(fs.q(p.k[i] - 1, p.k[i + 1] - 1));
  return LongSum(terms) + static_cast<double>(jumps);
}

double LogBridge(const TreeView& tr, const FrozenSet& fs, int v, const Path& p,
                 const Pieces& pc) {
  return LogPathExact(tr, fs, v, p) - PathWeight(fs, pc, p) -
    LogTransfer(pc, p.k.front(), true, fs.K)[p.k.back() - 1];
}

std::vector<double> TipLogWeights(const TreeView& tr, const FrozenSet& fs, int v,
                                  int kTop, const Pieces& pc) {
  std::vector<double> lw = LogTransfer(pc, kTop, true, fs.K);
  for (int k = 0; k < fs.K; ++k) {
    lw[k] += std::log(fs.TipValue(tr.tipType[v], tr.age[v], k));
  }
  return lw;
}

// logQ: log density of the tip bin and the bridge, as .DrawTip / .LogTipQ
double TipLogQ(const TreeView& tr, const FrozenSet& fs, int v, const Path& p,
               const Pieces& pc, const std::vector<double>& lw) {
  return LogBridge(tr, fs, v, p, pc) + lw[p.k.back() - 1] - LogSumExp(lw);
}

}  // namespace

// [[Rcpp::export(name = ".LogTransferRowCpp")]]
NumericVector LogTransferRowCpp(List pc, int k) {
  const Pieces pieces = PiecesFromR(pc);
  const std::vector<double> r = LogTransfer(pieces, k, true, PiecesK(pieces));
  return NumericVector(r.begin(), r.end());
}

// [[Rcpp::export(name = ".LogTransferColCpp")]]
NumericVector LogTransferColCpp(List pc, int j) {
  const Pieces pieces = PiecesFromR(pc);
  const std::vector<double> r = LogTransfer(pieces, j, false, PiecesK(pieces));
  return NumericVector(r.begin(), r.end());
}

// [[Rcpp::export(name = ".LogPathExactCpp")]]
double LogPathExactCpp(List tree, List set, int v, List path) {
  return LogPathExact(TreeView(tree), FrozenSet(set, false), v - 1, PathFromR(path));
}

// [[Rcpp::export(name = ".TipLogWeightsCpp")]]
NumericVector TipLogWeightsCpp(List tree, List set, int v, int kTop, List pc) {
  const FrozenSet fs(set, false);
  const std::vector<double> lw = TipLogWeights(TreeView(tree), fs, v - 1, kTop,
                                               PiecesFromR(pc));
  return NumericVector(lw.begin(), lw.end());
}

// .NodeLogWeights(); ends NA for a fossil tip whose bin is marginalized
// [[Rcpp::export(name = ".NodeLogWeightsCpp")]]
NumericVector NodeLogWeightsCpp(List tree, List set, int v, NumericVector factor,
                                int kTop, IntegerVector kids, NumericVector ends) {
  const TreeView tr(tree);
  const FrozenSet fs(set);
  const int K = fs.K;
  std::vector<double> lw = LogTransfer(MakePieces(tr, fs, v - 1), kTop, true, K);
  for (int k = 0; k < K; ++k) lw[k] = std::log(factor[k]) + lw[k];
  std::vector<double> tipD(K), x(K);
  for (int i = 0; i < kids.size(); ++i) {
    const int kid = kids[i] - 1;
    const Pieces pc = MakePieces(tr, fs, kid);
    if (ISNAN(ends[i])) {
      for (int l = 0; l < K; ++l) {
        tipD[l] = std::log(fs.TipValue(tr.tipType[kid], tr.age[kid], l));
      }
      for (int l = 0; l < K; ++l) {
        const std::vector<double> row = LogTransfer(pc, l + 1, true, K);
        for (int k = 0; k < K; ++k) x[k] = row[k] + tipD[k];
        lw[l] += LogSumExp(x);
      }
    } else {
      const int end = static_cast<int>(ends[i]);
      const std::vector<double> col = LogTransfer(pc, end, false, K);
      for (int k = 0; k < K; ++k) lw[k] += col[k];
    }
  }
  return NumericVector(lw.begin(), lw.end());
}

// Bridges on branches `nodes` between the fixed bins `bin`; returns the
// updated paths and their log proposal density
// [[Rcpp::export(name = ".DrawBridgesCpp")]]
List DrawBridgesCpp(List tree, List set, IntegerVector nodes, IntegerVector bin,
                    List paths) {
  const TreeView tr(tree);
  const FrozenSet fs(set);
  int len = paths.size();
  for (int u : nodes) len = std::max(len, u);
  List out(len);
  for (int u = 0; u < paths.size(); ++u) out[u] = paths[u];
  double logQ = 0.0;
  try {
    for (int u1 : nodes) {
      const int u = u1 - 1;
      const Pieces pc = MakePieces(tr, fs, u);
      const Path p = SampleBranch(fs, pc, nullptr, bin[tr.parent[u]], bin[u]);
      out[u] = PathToR(p);
      logQ = logQ + LogBridge(tr, fs, u, p, pc);
    }
  } catch (const InfeasibleProposal& e) {
    return Infeasible(e);
  }
  return List::create(_["paths"] = out, _["logQ"] = logQ);
}

// Sum of the log bridge densities of the paths on branches `nodes`
// [[Rcpp::export(name = ".LogBridgesCpp")]]
double LogBridgesCpp(List tree, List set, IntegerVector nodes, List paths) {
  const TreeView tr(tree);
  const FrozenSet fs(set);
  std::vector<double> terms(nodes.size());
  for (int i = 0; i < nodes.size(); ++i) {
    const int u = nodes[i] - 1;
    terms[i] = LogBridge(tr, fs, u, PathFromR(paths[u]), MakePieces(tr, fs, u));
  }
  return LongSum(terms);
}

// [[Rcpp::export(name = ".DrawTipCpp")]]
List DrawTipCpp(List tree, List set, int v, int kTop) {
  const TreeView tr(tree);
  const FrozenSet fs(set);
  const Pieces pc = MakePieces(tr, fs, v - 1);
  const std::vector<double> lw = TipLogWeights(tr, fs, v - 1, kTop, pc);
  const double top = RMax(lw);
  std::vector<double> w(lw.size());
  for (std::size_t k = 0; k < lw.size(); ++k) w[k] = std::exp(lw[k] - top);
  try {
    const int kTip = Sample1(w);
    const Path p = SampleBranch(fs, pc, nullptr, kTop, kTip);
    return List::create(_["path"] = PathToR(p), _["bin"] = kTip,
                        _["logQ"] = TipLogQ(tr, fs, v - 1, p, pc, lw));
  } catch (const InfeasibleProposal& e) {
    return Infeasible(e);
  }
}

// [[Rcpp::export(name = ".LogTipQCpp")]]
double LogTipQCpp(List tree, List set, int v, List path) {
  const TreeView tr(tree);
  const FrozenSet fs(set);
  const Path p = PathFromR(path);
  const Pieces pc = MakePieces(tr, fs, v - 1);
  return TipLogQ(tr, fs, v - 1, p, pc, TipLogWeights(tr, fs, v - 1, p.k.front(), pc));
}

// .Rebuild(): children, sampled ancestors and their ages, and postorder
// recomputed from `parent`, `age` and `tipType`
// [[Rcpp::export(name = ".RebuildCpp")]]
List RebuildCpp(List tree) {
  const int nTip = as<int>(tree["nTip"]), nNode = as<int>(tree["nNode"]);
  const int nAll = nTip + nNode, root = as<int>(tree["root"]) - 1;
  const IntegerVector par = tree["parent"];
  const CharacterVector tipType = tree["tipType"];
  std::vector<std::vector<int> > children(nAll);
  for (int v = 0; v < nAll; ++v) {
    if (par[v] != NA_INTEGER) children[par[v] - 1].push_back(v);
  }
  IntegerMatrix kids(nNode, 2);
  std::vector<int> isAncestor(nAll, 0);
  for (int v = 0; v < nTip; ++v) isAncestor[v] = tipType[v] == "ancestor";
  LogicalVector sampled(nNode);
  for (int i = 0; i < nNode; ++i) {
    const std::vector<int>& ch = children[nTip + i];
    if (ch.size() != 2) {
      Rcpp::stop("Node %d is not a bifurcation", nTip + i + 1);  // # nocov
    }
    kids(i, 0) = ch[0] + 1;
    kids(i, 1) = ch[1] + 1;
    sampled[i] = isAncestor[ch[0]] || isAncestor[ch[1]];
  }
  NumericVector age = clone(as<NumericVector>(tree["age"]));
  for (int v = 0; v < nTip; ++v) {
    if (isAncestor[v]) age[v] = age[par[v] - 1];
  }
  // .Preorder(): depth first, children in increasing order
  std::vector<int> order, stack(1, root);
  while (!stack.empty()) {
    const int v = stack.back();
    stack.pop_back();
    order.push_back(v);
    stack.insert(stack.end(), children[v].rbegin(), children[v].rend());
  }
  std::vector<int> post;
  for (auto it = order.rbegin(); it != order.rend(); ++it) {
    if (*it >= nTip) post.push_back(*it + 1);
  }
  List out = Rf_shallow_duplicate(tree);
  out["kids"] = kids;
  out["sampled"] = sampled;
  out["age"] = age;
  out["postorder"] = IntegerVector(post.begin(), post.end());
  return out;
}
