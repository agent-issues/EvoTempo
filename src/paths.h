#ifndef EVOTEMPO_PATHS_H
#define EVOTEMPO_PATHS_H

// Path helpers defined in paths.cpp and shared with tree_moves.cpp

#include <Rcpp.h>
#include "frozen.h"
#include "tree.h"
#include <vector>

namespace EvoTempo {

struct Path {
  std::vector<double> t;
  std::vector<int> k;  // one-based bins
};

typedef std::vector<std::vector<double> > Beta;

int Sample1(const std::vector<double>& w);
double LongSum(const std::vector<double>& x);
Rcpp::List PathToR(const Path& p);
Path PathFromR(SEXP x);
Pieces MakePieces(const TreeView& tr, const FrozenSet& set, int v);
Pieces PiecesFromR(SEXP x);
void MatVec(const std::vector<double>& M, const std::vector<double>& x,
            std::vector<double>& out, int K);
Path SampleBranch(const FrozenSet& set, const Pieces& pc, const Beta* betaIn,
                  int kTop, int end);
double PathWeight(const FrozenSet& set, const Pieces& pc, const Path& path);
Rcpp::List Infeasible(const InfeasibleProposal& e);

}  // namespace EvoTempo

#endif  // EVOTEMPO_PATHS_H
