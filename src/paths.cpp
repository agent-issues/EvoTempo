// Tempo paths on a fixed time tree: C++ engine for the frozen messages,
// forward-filtering backward-sampling with endpoint-conditioned bridges, and
// the exact augmented density.  Each function reproduces its namesake in
// R/Paths.R (without the "Cpp" suffix), draws random numbers from R's
// generator in the same order, and reads and returns the same R structures,
// so the two engines are interchangeable within one chain.

#include <Rcpp.h>
#include "frozen.h"
#include "paths.h"
#include "tree.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;
using EvoTempo::FindInterval;
using EvoTempo::FrozenSet;
using EvoTempo::InfeasibleProposal;
using EvoTempo::Pieces;
using EvoTempo::TreeView;

namespace EvoTempo {

double Unif() { return R::runif(0.0, 1.0); }

// .Sample1(): one-based draw from weights w
int Sample1(const std::vector<double>& w) {
  std::vector<double> cw(w.size());
  long double total = 0.0;
  for (std::size_t i = 0; i < w.size(); ++i) {
    total += std::max(w[i], 0.0);
    cw[i] = static_cast<double>(total);
  }
  const double tot = cw.empty() ? 0.0 : cw.back();
  if (!std::isfinite(tot) || !(tot > 0)) {
    throw InfeasibleProposal("No positive finite weight to sample from");
  }
  return FindInterval(Unif() * tot, cw) + 1;
}

double LongSum(const std::vector<double>& x) {
  long double s = 0.0;
  for (double v : x) s += v;
  return static_cast<double>(s);
}

List PathToR(const Path& p) {
  return List::create(_["t"] = NumericVector(p.t.begin(), p.t.end()),
                      _["k"] = IntegerVector(p.k.begin(), p.k.end()));
}

Path PathFromR(SEXP x) {
  List l(x);
  const NumericVector t = l["t"];
  const IntegerVector k = l["k"];
  Path p;
  p.t.assign(t.begin(), t.end());
  p.k.assign(k.begin(), k.end());
  return p;
}

// Jump time of a bridge: the root of F(tau) = sum_k cf_k J_k(tau) - target on
// (0, remain), where J_k(tau) integrates exp(gaa u + lam_k (remain - u)) over
// [0, tau], rd_k = gaa - lam_k, eLam_k = exp(lam_k remain), and
// F(0) = -target < 0 < F(remain) = total - target.  F' is the jump-time density,
// sum_k cf_k exp(gaa tau + lam_k (remain - tau)), and F'' and F''' multiply its
// terms by rd_k and rd_k^2, so the derivatives come with F from the same
// transcendentals.  exp(lam_k remain) is not factored out of the rd_k > 0 terms:
// it underflows where expm1(rd_k tau) overflows.
//
// The iteration takes Halley steps where their correction to the Newton step is
// modest, Newton steps elsewhere, and bisects when a step leaves the bracket or
// fails to halve the last move.  It stops when the predicted error of the new
// iterate is below tol / 2: |step| for a Newton step, the cubic term of Halley's
// error for a Halley step within 1e-4 * span of the last iterate.  R's .Bridge()
// uses uniroot() at the same tolerance, so the engines agree to ~tol, not bit for
// bit.  Draws no random numbers.
double JumpTime(const std::vector<double>& cf, const std::vector<double>& lam,
                const std::vector<double>& rateDiff,
                const std::vector<double>& eLam, double gaa, double eGaa,
                double remain, double target, double total, double span) {
  const int K = cf.size();
  const double tol = 1e-12 * span;
  // Start with a Newton step from the nearer end, using F'(0) = sum cf eLam and
  // F'(remain) = eGaa sum cf, eGaa = exp(gaa remain); regula falsi if it leaves
  // (0, remain)
  const double u = target / total;
  double slope0 = 0.0, sumCf = 0.0;
  for (int k = 0; k < K; ++k) {
    slope0 += cf[k] * eLam[k];
    sumCf += cf[k];
  }
  double x = u < 0.5 ? target / slope0 :
    remain - (total - target) / (eGaa * sumCf);
  if (!(x > 0 && x < remain)) x = remain * u;
  double lo = 0.0, hi = remain, prev = remain;
  for (int it = 0; it < 100; ++it) {
    long double sumJ = 0.0;
    double d1 = 0.0, d2 = 0.0, d3 = 0.0;
    for (int k = 0; k < K; ++k) {
      const double rd = rateDiff[k];
      double e, jk;
      if (rd > 0) {
        e = std::exp(gaa * x + lam[k] * (remain - x));
        jk = e * -std::expm1(-rd * x) / rd;
      } else {
        const double m = std::expm1(rd * x);
        jk = rd < 0 ? eLam[k] * -m / -rd : x * eLam[k];
        e = eLam[k] * (1 + m);
      }
      sumJ += cf[k] * jk;
      const double ce = cf[k] * e;
      d1 += ce;
      d2 += ce * rd;
      d3 += ce * rd * rd;
    }
    const double f = static_cast<double>(sumJ) - target;
    if (f < 0) lo = x; else hi = x;
    double step = f / d1;
    const double corr = step * d2 / (2 * d1);
    const bool halley = corr > -1 && corr < 0.5;
    if (halley) step /= 1 - corr;
    double next = x - step;
    const bool ok = d1 > 0 && next >= lo && next <= hi &&
      std::fabs(step) < 0.5 * prev;
    if (!ok) next = (lo + hi) / 2;
    prev = std::fabs(next - x);
    const double c3 = (d2 / (2 * d1)) * (d2 / (2 * d1)) - d3 / (6 * d1);
    const double err = !ok ? R_PosInf :
      halley && std::fabs(step) < 1e-4 * span ?
        std::fabs(c3) * std::fabs(step) * step * step : std::fabs(step);
    if (hi - lo < tol || err < tol / 2) return next;
    x = next;
  }
  return (lo + hi) / 2;  // # nocov
}

// .Bridge(): endpoint-conditioned path in slice s (zero-based) from bin i to
// bin j (one-based) over `span`; jump times measured from the old end
Path Bridge(const FrozenSet& set, int s, double span, int i, int j) {
  const int K = set.K;
  const EvoTempo::Slice& sl = set.slices[s];
  const std::vector<double>& lam = sl.values;
  const std::vector<double>& V = sl.vectors;
  const std::vector<double>& hp = set.halfPi;
  Path out;
  out.k.push_back(i);
  int a = i;
  double now = 0.0;
  std::vector<double> rateDiff(K), jAll(K), eLam(K), w;
  std::vector<int> nbr;
  std::vector<std::vector<double> > coef;
  while (true) {
    const double remain = span - now;
    const int a0 = a - 1;
    const double gaa = set.q(a0, a0) +
      (2 * set.lambda[a0] * sl.eBar[a0] - set.sigma[a0]);
    // J_k(remain); for rd > 0 its exponential is exp(gaa remain), and
    // exp(lam_k remain) = exp(gaa remain) (1 + expm1(-rd remain)) is kept for
    // the root finder's start only
    const double eGaa = std::exp(gaa * remain);
    for (int k = 0; k < K; ++k) {
      const double rd = rateDiff[k] = gaa - lam[k];
      if (rd > 0) {
        const double m = std::expm1(-rd * remain);
        jAll[k] = eGaa * -m / rd;
        eLam[k] = eGaa * (1 + m);
      } else {
        eLam[k] = std::exp(lam[k] * remain);
        jAll[k] = rd < 0 ? eLam[k] * -std::expm1(rd * remain) / -rd :
          remain * eLam[k];
      }
    }
    nbr.clear();
    coef.clear();
    for (int c = 0; c < K; ++c) {
      if (set.q(a0, c) > 0) {
        nbr.push_back(c + 1);
        std::vector<double> cf(K);
        for (int k = 0; k < K; ++k) {
          cf[k] = set.q(a0, c) * hp[j - 1] / hp[c] * V[c + k * K] *
            V[(j - 1) + k * K];
        }
        coef.push_back(cf);
      }
    }
    w.assign(1, a == j ? eGaa : 0.0);
    std::vector<double> prod(K);
    for (const std::vector<double>& cf : coef) {
      for (int k = 0; k < K; ++k) prod[k] = cf[k] * jAll[k];
      w.push_back(LongSum(prod));
    }
    long double pos = 0.0;
    for (double x : w) pos += std::max(x, 0.0);
    if (!(static_cast<double>(pos) > 0)) {
      throw InfeasibleProposal("Infeasible bridge");
    }
    const int pick = Sample1(w);
    if (pick == 1) break;
    const double target = Unif() * w[pick - 1];
    now += JumpTime(coef[pick - 2], lam, rateDiff, eLam, gaa, eGaa, remain,
                    target, w[pick - 1], span);
    a = nbr[pick - 2];
    out.t.push_back(now);
    out.k.push_back(a);
  }
  return out;
}

// .Pieces(tree, set, v), v zero-based
Pieces MakePieces(const TreeView& tr, const FrozenSet& set, int v) {
  const double top = tr.age[tr.parent[v]], bottom = tr.age[v];
  Pieces pc;
  pc.edges.push_back(top);
  for (int i = static_cast<int>(set.bounds.size()) - 1; i >= 0; --i) {
    const double b = set.bounds[i];
    if (b > bottom + 1e-12 && b < top - 1e-12) pc.edges.push_back(b);
  }
  pc.edges.push_back(bottom);
  const int m = pc.edges.size() - 1;
  pc.slice.resize(m);
  pc.M.resize(m);
  for (int j = 0; j < m; ++j) {
    const double mid = (pc.edges[j + 1] + pc.edges[j]) / 2;
    pc.slice[j] = std::min(FindInterval(mid, set.bounds), set.nSlice);
    set.SliceExp(pc.slice[j] - 1, pc.edges[j] - pc.edges[j + 1], pc.M[j]);
  }
  return pc;
}

List PiecesToR(const Pieces& pc, int K) {
  List M(pc.M.size());
  for (std::size_t j = 0; j < pc.M.size(); ++j) {
    NumericMatrix m(K, K);
    std::copy(pc.M[j].begin(), pc.M[j].end(), m.begin());
    M[j] = m;
  }
  return List::create(_["edges"] = NumericVector(pc.edges.begin(), pc.edges.end()),
                      _["slice"] = IntegerVector(pc.slice.begin(), pc.slice.end()),
                      _["M"] = M);
}

Pieces PiecesFromR(SEXP x) {
  List l(x);
  Pieces pc;
  const NumericVector e = l["edges"];
  pc.edges.assign(e.begin(), e.end());
  const IntegerVector s = l["slice"];
  pc.slice.assign(s.begin(), s.end());
  List M = l["M"];
  pc.M.resize(M.size());
  for (int j = 0; j < M.size(); ++j) {
    NumericMatrix m = M[j];
    pc.M[j].assign(m.begin(), m.end());
  }
  return pc;
}

Beta BetaFromR(SEXP x) {
  List l(x);
  Beta b(l.size());
  for (int j = 0; j < l.size(); ++j) {
    const NumericVector v = l[j];
    b[j].assign(v.begin(), v.end());
  }
  return b;
}

// out = M x for K x K column-major M
void MatVec(const std::vector<double>& M, const std::vector<double>& x,
            std::vector<double>& out, int K) {
  out.assign(K, 0.0);
  for (int j = 0; j < K; ++j) {
    const double xj = x[j];
    for (int i = 0; i < K; ++i) out[i] += M[i + static_cast<std::size_t>(j) * K] * xj;
  }
}

// .SampleBranch(); end (one-based) <= 0 to draw the bottom bin from beta
Path SampleBranch(const FrozenSet& set, const Pieces& pc, const Beta* betaIn,
                  int kTop, int end) {
  const int K = set.K;
  const int m = pc.M.size();
  Beta local;
  const Beta* beta = betaIn;
  if (end > 0) {
    local.assign(m + 1, std::vector<double>(K, 0.0));
    local[m][end - 1] = 1.0;
    std::vector<double> x;
    for (int j = m - 1; j >= 0; --j) {
      MatVec(pc.M[j], local[j + 1], x, K);
      const double mx = *std::max_element(x.begin(), x.end());
      for (int k = 0; k < K; ++k) local[j][k] = x[k] / mx;
    }
    beta = &local;
  }
  std::vector<int> b(m + 1);
  b[0] = kTop;
  std::vector<double> w(K);
  for (int j = 0; j < m; ++j) {
    const std::vector<double>& M = pc.M[j];
    for (int k = 0; k < K; ++k) {
      w[k] = M[(b[j] - 1) + static_cast<std::size_t>(k) * K] * (*beta)[j + 1][k];
    }
    b[j + 1] = Sample1(w);
  }
  Path out;
  out.k.push_back(kTop);
  for (int j = 0; j < m; ++j) {
    const double span = pc.edges[j] - pc.edges[j + 1];
    const Path br = Bridge(set, pc.slice[j] - 1, span, b[j], b[j + 1]);
    for (double t : br.t) out.t.push_back(pc.edges[j] - t);
    out.k.insert(out.k.end(), br.k.begin() + 1, br.k.end());
  }
  return out;
}

// Messages in C++ form
struct Msg {
  std::vector<double> D;  // nAll x K column-major
  std::vector<Pieces> pieces;
  std::vector<Beta> beta;
  std::vector<int> has;
  double logScale = 0.0;
};

Msg MsgFromR(List msg, int nAll) {
  Msg out;
  const NumericMatrix D = msg["D"];
  out.D.assign(D.begin(), D.end());
  List pcs = msg["pieces"], bet = msg["beta"];
  out.pieces.resize(nAll);
  out.beta.resize(nAll);
  out.has.assign(nAll, 0);
  for (int v = 0; v < nAll && v < pcs.size(); ++v) {
    if (Rf_isNull(pcs[v])) continue;
    out.pieces[v] = PiecesFromR(pcs[v]);
    out.beta[v] = BetaFromR(bet[v]);
    out.has[v] = 1;
  }
  out.logScale = as<double>(msg["logScale"]);
  return out;
}

struct Draw {
  std::vector<Path> paths;
  std::vector<int> hasPath;
  std::vector<int> bin;  // one-based
};

void SampleClade(const TreeView& tr, const FrozenSet& set, const Msg& msg,
                 Draw& d, int v, int kTop) {
  if (!msg.has[v]) Rcpp::stop("No pieces for branch %d", v + 1);
  d.paths[v] = SampleBranch(set, msg.pieces[v], &msg.beta[v], kTop, 0);
  d.hasPath[v] = 1;
  d.bin[v] = d.paths[v].k.back();
  if (v >= tr.nTip) {
    for (int c = 0; c < 2; ++c) {
      const int kid = tr.kids[2 * (v - tr.nTip) + c];
      if (tr.IsAncestorTip(kid)) {
        d.bin[kid] = d.bin[v];
      } else {
        SampleClade(tr, set, msg, d, kid, d.bin[v]);
      }
    }
  }
}

List DrawToR(const Draw& d, const List* pathsIn) {
  const int nAll = d.bin.size();
  List paths = pathsIn ? clone(*pathsIn) : List(nAll);
  if (paths.size() < nAll) {
    List grown(nAll);
    for (int v = 0; v < paths.size(); ++v) grown[v] = paths[v];
    paths = grown;
  }
  for (int v = 0; v < nAll; ++v) {
    if (d.hasPath[v]) paths[v] = PathToR(d.paths[v]);
  }
  return List::create(_["paths"] = paths,
                      _["bin"] = IntegerVector(d.bin.begin(), d.bin.end()));
}

List Infeasible(const InfeasibleProposal& e) {
  return List::create(_["infeasible"] = std::string(e.what()));
}

double PathWeight(const FrozenSet& set, const Pieces& pc, const Path& path) {
  std::vector<double> cuts(pc.edges);
  cuts.insert(cuts.end(), path.t.begin(), path.t.end());
  std::sort(cuts.begin(), cuts.end(), std::greater<double>());
  cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
  // -c(edges[1], path$t), increasing
  std::vector<double> neg(1, -pc.edges[0]);
  for (double t : path.t) neg.push_back(-t);
  const int n = cuts.size() - 1;
  std::vector<double> terms(n);
  for (int i = 0; i < n; ++i) {
    const double hi = cuts[i], lo = cuts[i + 1], mid = (hi + lo) / 2;
    const int k = path.k[FindInterval(-mid, neg) - 1] - 1;
    const int s = std::min(FindInterval(mid, set.bounds), set.nSlice) - 1;
    const double eBar = set.slices[s].eBar[k];
    terms[i] = 2 * set.lambda[k] * (set.Hermite(hi, k, true) -
                                     set.Hermite(lo, k, true) - eBar * (hi - lo));
  }
  return LongSum(terms);
}

}  // namespace EvoTempo

using namespace EvoTempo;

// [[Rcpp::export(name = ".PiecesCpp")]]
List PiecesCpp(List tree, List set, int v) {
  const TreeView tr(tree);
  const FrozenSet fs(set);
  return PiecesToR(MakePieces(tr, fs, v - 1), fs.K);
}

// [[Rcpp::export(name = ".SliceExpCpp")]]
NumericMatrix SliceExpCpp(List set, int s, double span) {
  const FrozenSet fs(set);
  std::vector<double> out;
  fs.SliceExp(s - 1, span, out);
  NumericMatrix m(fs.K, fs.K);
  std::copy(out.begin(), out.end(), m.begin());
  return m;
}

// [[Rcpp::export(name = ".MessagesCpp")]]
List MessagesCpp(List tree, List set) {
  const TreeView tr(tree);
  const FrozenSet fs(set);
  const int K = fs.K, nTip = tr.nTip, nAll = tr.nAll;
  std::vector<double> D(static_cast<std::size_t>(nAll) * K, 0.0);
  for (int i = 0; i < nTip; ++i) {
    for (int k = 0; k < K; ++k) {
      D[i + static_cast<std::size_t>(k) * nAll] = fs.TipValue(tr.tipType[i], tr.age[i], k);
    }
  }
  List pieces(nAll), beta(nAll);
  double logScale = 0.0;
  std::vector<double> x;
  auto Up = [&](int v) {
    const Pieces pc = MakePieces(tr, fs, v);
    const int m = pc.M.size();
    std::vector<std::vector<double> > b(m + 1, std::vector<double>(K));
    for (int k = 0; k < K; ++k) b[m][k] = D[v + static_cast<std::size_t>(k) * nAll];
    for (int j = m - 1; j >= 0; --j) {
      MatVec(pc.M[j], b[j + 1], x, K);
      const double mx = *std::max_element(x.begin(), x.end());
      for (int k = 0; k < K; ++k) b[j][k] = x[k] / mx;
      logScale += std::log(mx);
    }
    pieces[v] = PiecesToR(pc, K);
    List bl(m + 1);
    for (int j = 0; j <= m; ++j) bl[j] = NumericVector(b[j].begin(), b[j].end());
    beta[v] = bl;
    return b[0];
  };
  std::vector<double> value(K);
  for (int v : tr.postorder) {
    const int k1 = tr.kids[2 * (v - nTip)], k2 = tr.kids[2 * (v - nTip) + 1];
    if (tr.sampled[v - nTip]) {
      const int cont = tr.IsAncestorTip(k1) ? k2 : k1;
      const std::vector<double> u = Up(cont);
      for (int k = 0; k < K; ++k) value[k] = fs.psi[k] * (1 - fs.removal) * u[k];
    } else {
      const std::vector<double> u1 = Up(k1);
      const std::vector<double> u2 = Up(k2);
      for (int k = 0; k < K; ++k) {
        value[k] = u1[k] * u2[k] * (v == tr.root ? 1.0 : fs.lambda[k]);
      }
    }
    const double mx = *std::max_element(value.begin(), value.end());
    for (int k = 0; k < K; ++k) D[v + static_cast<std::size_t>(k) * nAll] = value[k] / mx;
    logScale += std::log(mx);
  }
  NumericMatrix Dm(nAll, K);
  std::copy(D.begin(), D.end(), Dm.begin());
  return List::create(_["D"] = Dm, _["pieces"] = pieces, _["beta"] = beta,
                      _["logScale"] = logScale);
}

// [[Rcpp::export(name = ".BridgeCpp")]]
List BridgeCpp(List set, int s, double span, int i, int j) {
  const FrozenSet fs(set);
  try {
    return PathToR(Bridge(fs, s - 1, span, i, j));
  } catch (const InfeasibleProposal& e) {
    return Infeasible(e);
  }
}

// [[Rcpp::export(name = ".SampleBranchCpp")]]
List SampleBranchCpp(List set, List pc, SEXP beta, int kTop, SEXP end) {
  const FrozenSet fs(set);
  const Pieces pieces = PiecesFromR(pc);
  const int e = Rf_isNull(end) ? 0 : as<int>(end);
  Beta b;
  if (!e) b = BetaFromR(beta);
  try {
    return PathToR(SampleBranch(fs, pieces, &b, kTop, e));
  } catch (const InfeasibleProposal& ex) {
    return Infeasible(ex);
  }
}

// [[Rcpp::export(name = ".SampleTreeCpp")]]
List SampleTreeCpp(List tree, List set, List msg) {
  const TreeView tr(tree);
  const FrozenSet fs(set);
  const Msg m = MsgFromR(msg, tr.nAll);
  const int K = fs.K;
  Draw d;
  d.paths.resize(tr.nAll);
  d.hasPath.assign(tr.nAll, 0);
  d.bin.assign(tr.nAll, 0);
  try {
    std::vector<double> w(K);
    for (int k = 0; k < K; ++k) {
      w[k] = fs.prior[k] * m.D[tr.root + static_cast<std::size_t>(k) * tr.nAll];
      if (fs.joint) w[k] = w[k] / fs.survive[k];
    }
    d.bin[tr.root] = Sample1(w);
    const int r = tr.root - tr.nTip;
    for (int c = 0; c < 2; ++c) {
      SampleClade(tr, fs, m, d, tr.kids[2 * r + c], d.bin[tr.root]);
    }
  } catch (const InfeasibleProposal& e) {
    return Infeasible(e);
  }
  return DrawToR(d, nullptr);
}

// [[Rcpp::export(name = ".SampleCladeCpp")]]
List SampleCladeCpp(List tree, List set, List msg, List paths, IntegerVector bin,
                    int v, int kTop) {
  const TreeView tr(tree);
  const FrozenSet fs(set);
  const Msg m = MsgFromR(msg, tr.nAll);
  Draw d;
  d.paths.resize(tr.nAll);
  d.hasPath.assign(tr.nAll, 0);
  d.bin.assign(bin.begin(), bin.end());
  d.bin.resize(tr.nAll, 0);
  try {
    SampleClade(tr, fs, m, d, v - 1, kTop);
  } catch (const InfeasibleProposal& e) {
    return Infeasible(e);
  }
  return DrawToR(d, &paths);
}

// [[Rcpp::export(name = ".PathWeightCpp")]]
double PathWeightCpp(List set, List pc, List path) {
  const FrozenSet fs(set);
  return PathWeight(fs, PiecesFromR(pc), PathFromR(path));
}

// [[Rcpp::export(name = ".LogWeightCpp")]]
double LogWeightCpp(List tree, List set, List msg, List paths, IntegerVector nodes) {
  const FrozenSet fs(set);
  List pcs = msg["pieces"];
  std::vector<double> terms(nodes.size());
  for (int i = 0; i < nodes.size(); ++i) {
    const int v = nodes[i] - 1;
    terms[i] = PathWeight(fs, PiecesFromR(pcs[v]), PathFromR(paths[v]));
  }
  (void) tree;
  return LongSum(terms);
}

// .LogAugmented(): exact log augmented density, start at the root
// [[Rcpp::export(name = ".LogAugmentedCpp")]]
double LogAugmentedCpp(List tree, List set, List paths, IntegerVector bin) {
  const TreeView tr(tree);
  const FrozenSet fs(set, false);
  if (fs.condition == "extant" && !tr.CrownCompatible()) {
    return -std::numeric_limits<double>::infinity();
  }
  const int nTip = tr.nTip;
  std::vector<double> part;
  long double jumpSum = 0.0;
  bool anyJump = false;
  for (int v : tr.PathNodes()) {
    const Path p = PathFromR(paths[v]);
    const int n = p.k.size();
    for (int i = 0; i < n; ++i) {
      const double hi = i == 0 ? tr.age[tr.parent[v]] : p.t[i - 1];
      const double lo = i == n - 1 ? tr.age[v] : p.t[i];
      const int k = p.k[i] - 1;
      part.push_back((fs.q(k, k) - fs.sigma[k]) * (hi - lo) +
                     2 * fs.lambda[k] * (fs.Hermite(hi, k, true) -
                                         fs.Hermite(lo, k, true)));
    }
    for (int i = 0; i + 1 < n; ++i) {
      jumpSum += std::log(fs.q(p.k[i] - 1, p.k[i + 1] - 1));
      anyJump = true;
    }
  }
  const double pathPart = LongSum(part) + (anyJump ? static_cast<double>(jumpSum) : 0.0);
  std::vector<double> lamTerms, psiTerms, tipTerms;
  for (int v : tr.postorder) {
    const int b = bin[v] - 1;
    if (tr.sampled[v - nTip]) {
      psiTerms.push_back(std::log(fs.psi[b] * (1 - fs.removal)));
    } else if (v != tr.root) {
      lamTerms.push_back(std::log(fs.lambda[b]));
    }
  }
  for (int i = 0; i < nTip; ++i) {
    if (tr.tipType[i] == EvoTempo::kAncestor) continue;
    tipTerms.push_back(std::log(fs.TipValue(tr.tipType[i], tr.age[i], bin[i] - 1)));
  }
  const int rootBin = bin[tr.root] - 1;
  double condition;
  if (fs.joint) {
    condition = std::log(fs.survive[rootBin]);
  } else {
    std::vector<double> ps(fs.K);
    for (int k = 0; k < fs.K; ++k) ps[k] = fs.prior[k] * fs.survive[k];
    condition = std::log(LongSum(ps));
  }
  return std::log(fs.prior[rootBin]) + pathPart + LongSum(lamTerms) +
    LongSum(psiTerms) + LongSum(tipTerms) - condition;
}

// Tempo-weighted branch lengths (.CharacterLengths)
// [[Rcpp::export(name = ".CharacterLengthsCpp")]]
NumericVector CharacterLengthsCpp(List tree, NumericVector x, List paths,
                                  double rate, double gamma) {
  const TreeView tr(tree);
  NumericVector out(tr.nAll);
  std::vector<double> terms;
  for (int v : tr.PathNodes()) {
    const Path p = PathFromR(paths[v]);
    const int n = p.k.size();
    terms.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
      const double hi = i == 0 ? tr.age[tr.parent[v]] : p.t[i - 1];
      const double lo = i == n - 1 ? tr.age[v] : p.t[i];
      terms[i] = std::exp(gamma * x[p.k[i] - 1]) * (hi - lo);
    }
    out[v] = rate * LongSum(terms);
  }
  return out;
}
