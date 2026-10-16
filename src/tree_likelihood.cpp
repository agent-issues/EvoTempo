// Marginal dCET tree likelihood and extinction probabilities: C++ engine
// for TreeLikelihood(), ExtinctionProbability() and .EGrid().
// R/TreeLikelihood.R, R/Extinction.R and R/Paths.R hold the reference
// implementations that these reproduce.

#include <Rcpp.h>
#include "dcet_ode.h"
#include "tree.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

using namespace Rcpp;
using EvoTempo::EDProblem;
using EvoTempo::EGProblem;
using EvoTempo::EProblem;
using EvoTempo::BandGenerator;
using EvoTempo::Radau5;

namespace {

// E at each of `times` (non-negative, non-decreasing), from E(0) = 1 - rho;
// out is row-major, times.size() x K
void ExtinctionAt(const BandGenerator& Q, const std::vector<double>& lambda,
                  const std::vector<double>& mu, const std::vector<double>& sigma,
                  double rho, const double* times, int nTime, double rtol,
                  double atol, std::vector<double>& out) {
  const int K = Q.n;
  EProblem prob(Q, lambda.data(), mu.data(), sigma.data());
  Radau5 solver(rtol, atol);
  std::vector<double> e(K, 1.0 - rho);
  out.assign(static_cast<std::size_t>(nTime) * K, 0.0);
  double now = 0.0;
  for (int i = 0; i < nTime; ++i) {
    if (times[i] > now) {
      if (!solver.Integrate(prob, e.data(), times[i] - now)) {
        Rcpp::stop("ODE solver failed");
      }
      now = times[i];
    }
    std::copy(e.begin(), e.end(), out.begin() + static_cast<std::size_t>(i) * K);
  }
}

std::vector<double> AsVector(const NumericVector& x, int n) {
  std::vector<double> out(n);
  for (int i = 0; i < n; ++i) out[i] = x[x.size() == 1 ? 0 : i];
  return out;
}

double LongSum(const std::vector<double>& x) {
  long double s = 0.0;
  for (double v : x) s += v;
  return static_cast<double>(s);
}

}  // namespace

// [[Rcpp::export(name = ".ExtinctionCpp")]]
NumericMatrix ExtinctionCpp(NumericVector times, NumericMatrix Q,
                            NumericVector lambda, NumericVector mu,
                            NumericVector psi, double rho, double rtol,
                            double atol) {
  const int K = Q.nrow();
  const BandGenerator gen(Q.begin(), K);
  const std::vector<double> lam = AsVector(lambda, K), m = AsVector(mu, K),
    ps = AsVector(psi, K);
  std::vector<double> sigma(K);
  for (int k = 0; k < K; ++k) sigma[k] = lam[k] + m[k] + ps[k];
  std::vector<double> out;
  ExtinctionAt(gen, lam, m, sigma, rho, times.begin(), times.size(), rtol, atol,
               out);
  NumericMatrix ret(times.size(), K);
  for (int i = 0; i < times.size(); ++i) {
    for (int k = 0; k < K; ++k) ret(i, k) = out[static_cast<std::size_t>(i) * K + k];
  }
  return ret;
}

// E, G = integral of E, and dE/dt at sorted ages starting at 0 (.EGrid)
// [[Rcpp::export(name = ".EGridCpp")]]
List EGridCpp(NumericVector ages, NumericMatrix Q, NumericVector lambda,
              NumericVector mu, NumericVector psi, double rho, double rtol,
              double atol) {
  const int K = Q.nrow();
  const BandGenerator gen(Q.begin(), K);
  const std::vector<double> lam = AsVector(lambda, K), m = AsVector(mu, K),
    ps = AsVector(psi, K);
  std::vector<double> sigma(K);
  for (int k = 0; k < K; ++k) sigma[k] = lam[k] + m[k] + ps[k];
  EGProblem prob(gen, lam.data(), m.data(), sigma.data());
  Radau5 solver(rtol, atol);
  std::vector<double> y(2 * K, 0.0), f(2 * K);
  for (int k = 0; k < K; ++k) y[k] = 1.0 - rho;
  const int nAge = ages.size();
  NumericMatrix E(nAge, K), G(nAge, K), dE(nAge, K);
  double now = 0.0;
  for (int i = 0; i < nAge; ++i) {
    if (ages[i] > now) {
      if (!solver.Integrate(prob, y.data(), ages[i] - now)) {
        Rcpp::stop("ODE solver failed");
      }
      now = ages[i];
    }
    prob.RhsE(y.data(), f.data());
    for (int k = 0; k < K; ++k) {
      E(i, k) = y[k];
      G(i, k) = y[K + k];
      dE(i, k) = f[k];
    }
  }
  return List::create(_["E"] = E, _["G"] = G, _["dE"] = dE);
}

// .PruneTree() numerics, after its argument checks and .Compatible()
// [[Rcpp::export(name = ".PruneTreeCpp")]]
double PruneTreeCpp(List tree, NumericMatrix Q, NumericVector lambda,
                    NumericVector mu, NumericVector psi, double rho,
                    double removal, NumericVector prior, std::string condition,
                    bool joint, bool origin, double rtol, double atol,
                    Nullable<NumericMatrix> tipMask) {
  const EvoTempo::TreeView tr(tree);
  const int K = Q.nrow();
  const BandGenerator gen(Q.begin(), K);
  const std::vector<double> lam = AsVector(lambda, K), m = AsVector(mu, K),
    ps = AsVector(psi, K);
  std::vector<double> sigma(K);
  for (int k = 0; k < K; ++k) sigma[k] = lam[k] + m[k] + ps[k];
  const bool masked = tipMask.isNotNull();
  NumericMatrix mask = masked ? NumericMatrix(tipMask.get()) : NumericMatrix(1, 1);
  const int nTip = tr.nTip, nAll = tr.nAll;
  const double topAge = tr.age[tr.root] + (origin ? tr.stem : 0.0);

  // E at every node age and the top
  std::vector<double> sorted(tr.age.begin(), tr.age.end());
  sorted.push_back(topAge);
  std::sort(sorted.begin(), sorted.end());
  sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
  std::vector<double> eAt;
  ExtinctionAt(gen, lam, m, sigma, rho, sorted.data(), sorted.size(), rtol, atol,
               eAt);
  auto E = [&](double a) -> const double* {
    const std::size_t i = std::lower_bound(sorted.begin(), sorted.end(), a) -
      sorted.begin();
    return eAt.data() + i * K;
  };

  std::vector<double> D(static_cast<std::size_t>(nAll) * K, 0.0);
  auto Dv = [&](int v) { return D.data() + static_cast<std::size_t>(v) * K; };
  for (int i = 0; i < nTip; ++i) {
    double* d = Dv(i);
    if (tr.tipType[i] == EvoTempo::kExtant) {
      for (int k = 0; k < K; ++k) d[k] = rho;
    } else if (tr.tipType[i] == EvoTempo::kFossil) {
      const double* e = E(tr.age[i]);
      for (int k = 0; k < K; ++k) d[k] = ps[k] * (removal + (1.0 - removal) * e[k]);
    }
    if (masked) {
      for (int k = 0; k < K; ++k) d[k] *= mask(i, k);
    }
  }

  EDProblem prob(gen, lam.data(), m.data(), sigma.data());
  std::vector<double> y(2 * K);
  // D after integrating node v's D up to toAge
  auto Branch = [&](int v, double toAge, std::vector<double>& out) {
    const double* d0 = Dv(v);
    out.assign(d0, d0 + K);
    const double span = toAge - tr.age[v];
    const double scale = *std::max_element(d0, d0 + K);
    if (span <= 0 || !(scale > 0)) return;
    const double* e0 = E(tr.age[v]);
    for (int k = 0; k < K; ++k) {
      y[k] = e0[k];
      y[K + k] = d0[k] / scale;
    }
    Radau5 solver(rtol, atol);
    if (!solver.Integrate(prob, y.data(), span)) Rcpp::stop("ODE solver failed");
    for (int k = 0; k < K; ++k) out[k] = std::max(y[K + k], 0.0) * scale;
  };

  double logScale = 0.0;
  std::vector<double> b1, b2, value(K);
  for (int v : tr.postorder) {
    const int k1 = tr.kids[2 * (v - nTip)], k2 = tr.kids[2 * (v - nTip) + 1];
    if (tr.sampled[v - nTip]) {
      const bool sa1 = tr.IsAncestorTip(k1);
      const int cont = sa1 ? k2 : k1, sa = sa1 ? k1 : k2;
      Branch(cont, tr.age[v], b1);
      for (int k = 0; k < K; ++k) {
        value[k] = ps[k] * (1.0 - removal) * b1[k] * (masked ? mask(sa, k) : 1.0);
      }
    } else {
      Branch(k1, tr.age[v], b1);
      Branch(k2, tr.age[v], b2);
      const bool noLambda = v == tr.root && !origin;
      for (int k = 0; k < K; ++k) {
        value[k] = b1[k] * b2[k] * (noLambda ? 1.0 : lam[k]);
      }
    }
    const double scale = *std::max_element(value.begin(), value.end());
    if (!(scale > 0)) return -std::numeric_limits<double>::infinity();
    double* d = Dv(v);
    for (int k = 0; k < K; ++k) d[k] = value[k] / scale;
    logScale += std::log(scale);
  }
  std::vector<double> top(Dv(tr.root), Dv(tr.root) + K);
  if (origin) Branch(tr.root, topAge, top);

  std::vector<double> survive(K, 1.0);
  if (condition != "none") {
    if (condition == "extant") {
      std::vector<double> mu2(K), sigma2(K), ex;
      for (int k = 0; k < K; ++k) {
        mu2[k] = m[k] + removal * ps[k];
        sigma2[k] = lam[k] + mu2[k];
      }
      ExtinctionAt(gen, lam, mu2, sigma2, rho, &topAge, 1, rtol, atol, ex);
      for (int k = 0; k < K; ++k) survive[k] = 1.0 - ex[k];
    } else {
      const double* e = E(topAge);
      for (int k = 0; k < K; ++k) survive[k] = 1.0 - e[k];
    }
    if (!origin) {
      for (int k = 0; k < K; ++k) survive[k] *= survive[k];
    }
  }
  std::vector<double> terms(K);
  if (joint) {
    for (int k = 0; k < K; ++k) terms[k] = prior[k] * top[k] / survive[k];
    return std::log(LongSum(terms)) + logScale;
  }
  for (int k = 0; k < K; ++k) terms[k] = prior[k] * top[k];
  const double num = LongSum(terms);
  for (int k = 0; k < K; ++k) terms[k] = prior[k] * survive[k];
  return std::log(num) + logScale - std::log(LongSum(terms));
}

// Solver counters since the last call, then reset (profiling)
// [[Rcpp::export(name = ".SolverCountsCpp")]]
IntegerVector SolverCountsCpp() {
  EvoTempo::SolverCounter& n = EvoTempo::Counter();
  IntegerVector out = IntegerVector::create(
    _["steps"] = n.steps, _["rejects"] = n.rejects, _["newton"] = n.newton,
    _["rhs"] = n.rhs);
  n = EvoTempo::SolverCounter();
  return out;
}
