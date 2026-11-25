// [[Rcpp::depends(MkPrime)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <MkPrime/MkPrime.h>
#include <chrono>
#include <memory>
#include <vector>
using namespace Rcpp;
typedef MkPrime::Pruner<MkPrime::JCModel> P;
typedef std::chrono::steady_clock Clock;

// Per-phase time (ms per evaluation of all matrices) of the proposal cycle
// accept/reject, checkpoint, set_tree, site_likelihoods
// mode 0: one tip edge changes; 1: all edges; 2: one tip edge via set_edge_length
// [[Rcpp::export]]
NumericVector Phases(IntegerVector parent, List matrices, IntegerVector nState,
                     NumericVector lengths, int reps, int mode, bool reject) {
  std::vector<int> par, child;
  std::vector<double> len;
  for (int v = 0; v < parent.size(); ++v) {
    if (parent[v] == NA_INTEGER) continue;
    par.push_back(parent[v]); child.push_back(v + 1); len.push_back(lengths[v]);
  }
  const int nEdge = len.size();
  std::vector<std::unique_ptr<P> > pr;
  std::vector<std::vector<double> > site;
  for (int i = 0; i < matrices.size(); ++i) {
    IntegerMatrix m = matrices[i];
    std::vector<int> st(m.begin(), m.end());
    for (int& s : st) if (s == NA_INTEGER || s < 0) s = -1;
    pr.emplace_back(new P(MkPrime::JCModel(nState[i]), m.nrow(), m.ncol(), st.data()));
    site.emplace_back(m.ncol());
    pr.back()->set_tree(par.data(), child.data(), len.data(), nEdge);
    pr.back()->site_likelihoods(site.back().data());
  }
  double tEnd = 0, tCheck = 0, tSet = 0, tLik = 0, nRec = 0;
  bool open = false;
  std::vector<double> l2 = len;
  for (int r = 0; r < reps; ++r) {
    const int e = r % 50;  // tips come first
    if (mode == 1) {
      for (int j = 0; j < nEdge; ++j) l2[j] = len[j] * (1 + 1e-6 * (r + 1));
    } else {
      l2 = len;
      l2[e] = len[e] * (1 + 1e-3 * (1 + r % 7));
    }
    Clock::time_point t0 = Clock::now();
    if (open) for (auto& p : pr) { if (reject) p->reject(); else p->accept(); }
    Clock::time_point t1 = Clock::now();
    for (auto& p : pr) p->checkpoint();
    Clock::time_point t2 = Clock::now();
    if (mode == 2) {
      for (auto& p : pr) p->set_edge_length(e, l2[e]);
    } else {
      for (auto& p : pr) p->set_tree(par.data(), child.data(), l2.data(), nEdge);
    }
    Clock::time_point t3 = Clock::now();
    for (std::size_t i = 0; i < pr.size(); ++i) {
      pr[i]->site_likelihoods(site[i].data());
      nRec += pr[i]->n_recomputed();
    }
    Clock::time_point t4 = Clock::now();
    open = true;
    tEnd += std::chrono::duration<double, std::milli>(t1 - t0).count();
    tCheck += std::chrono::duration<double, std::milli>(t2 - t1).count();
    tSet += std::chrono::duration<double, std::milli>(t3 - t2).count();
    tLik += std::chrono::duration<double, std::milli>(t4 - t3).count();
  }
  return NumericVector::create(_["acceptOrReject"] = tEnd / reps, _["checkpoint"] = tCheck / reps,
                               _["setTree"] = tSet / reps, _["siteLikelihoods"] = tLik / reps,
                               _["nodesRecomputed"] = nRec / reps);
}
