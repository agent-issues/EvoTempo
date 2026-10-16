#ifndef EVOTEMPO_DCET_ODE_H
#define EVOTEMPO_DCET_ODE_H

// Backward equations of dCET (spec equations 2.1 and 3.1) for the Radau5
// integrator: E alone, E with D, and E with G = integral of E.
// Generators are dense column-major matrices with known bandwidths, so a
// tridiagonal tempo generator costs O(K) per solve and a dense one O(K^3).

#include "band.h"
#include "radau5.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

namespace EvoTempo {

// Generator with lower and upper bandwidths
struct Generator {
  int n;
  int kl;
  int ku;
  std::vector<double> q;  // column-major n x n

  Generator() : n(0), kl(0), ku(0) {}
  Generator(const double* qq, int nn) : n(nn), kl(0), ku(0), q(qq, qq + nn * nn) {
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        if (i != j && q[i + j * n] != 0.0) {
          if (i > j) kl = std::max(kl, i - j);
          if (j > i) ku = std::max(ku, j - i);
        }
      }
    }
  }

  double operator()(int i, int j) const { return q[i + j * n]; }

  // out = Q x
  void Multiply(const double* x, double* out) const {
    for (int i = 0; i < n; ++i) {
      const int lo = std::max(0, i - kl), hi = std::min(n - 1, i + ku);
      double s = 0.0;
      for (int j = lo; j <= hi; ++j) s += q[i + j * n] * x[j];
      out[i] = s;
    }
  }
};

// Shared part: A = Q + diag(g), g = 2 lambda E - sigma; LUs of s I - A
class ExtinctionCore {
 public:
  ExtinctionCore(const Generator& Q, const double* lambda, const double* mu,
                 const double* sigma)
    : Q_(Q), lambda_(lambda), mu_(mu), sigma_(sigma), g_(Q.n), qe_(Q.n) {
    real_.Resize(Q.n, Q.kl, Q.ku);
    cplx_.Resize(Q.n, Q.kl, Q.ku);
  }

  int K() const { return Q_.n; }

  void RhsE(const double* e, double* f) {
    Q_.Multiply(e, qe_.data());
    for (int k = 0; k < Q_.n; ++k) {
      f[k] = mu_[k] - sigma_[k] * e[k] + lambda_[k] * e[k] * e[k] + qe_[k];
    }
  }

  void JacobianE(const double* e) {
    for (int k = 0; k < Q_.n; ++k) g_[k] = 2.0 * lambda_[k] * e[k] - sigma_[k];
  }

  const double* g() const { return g_.data(); }

  bool Factor(double sr, std::complex<double> sc) {
    sr_ = sr;
    sc_ = sc;
    real_.Zero();
    cplx_.Zero();
    const int n = Q_.n;
    for (int j = 0; j < n; ++j) {
      const int lo = std::max(0, j - Q_.ku), hi = std::min(n - 1, j + Q_.kl);
      for (int i = lo; i <= hi; ++i) {
        const double a = Q_(i, j) + (i == j ? g_[i] : 0.0);
        real_(i, j) = (i == j ? sr : 0.0) - a;
        cplx_(i, j) = (i == j ? sc : std::complex<double>(0.0, 0.0)) - a;
      }
    }
    return real_.Factor() && cplx_.Factor();
  }

  void SolveRealK(double* b) const { real_.Solve(b); }
  void SolveComplexK(std::complex<double>* b) const { cplx_.Solve(b); }
  double sr() const { return sr_; }
  std::complex<double> sc() const { return sc_; }

 protected:
  const Generator& Q_;
  const double* lambda_;
  const double* mu_;
  const double* sigma_;
  std::vector<double> g_, qe_;
  BandLU<double> real_;
  BandLU<std::complex<double> > cplx_;
  double sr_ = 0.0;
  std::complex<double> sc_;
};

// dE/dt = mu - sigma E + lambda E^2 + Q E
class EProblem : public ExtinctionCore {
 public:
  using ExtinctionCore::ExtinctionCore;
  int Size() const { return K(); }
  void Rhs(const double* y, double* f) { RhsE(y, f); }
  void Jacobian(const double* y) { JacobianE(y); }
  void SolveReal(double* b) const { SolveRealK(b); }
  void SolveComplex(std::complex<double>* b) const { SolveComplexK(b); }
};

// (E, D): dD/dt = (2 lambda E - sigma) D + Q D
class EDProblem : public ExtinctionCore {
 public:
  EDProblem(const Generator& Q, const double* lambda, const double* mu,
            const double* sigma)
    : ExtinctionCore(Q, lambda, mu, sigma), c_(Q.n), qd_(Q.n),
      tmpC_(Q.n) {}
  int Size() const { return 2 * K(); }
  void Rhs(const double* y, double* f) {
    const int n = K();
    RhsE(y, f);
    Q_.Multiply(y + n, qd_.data());
    for (int k = 0; k < n; ++k) {
      f[n + k] = (2.0 * lambda_[k] * y[k] - sigma_[k]) * y[n + k] + qd_[k];
    }
  }
  void Jacobian(const double* y) {
    const int n = K();
    JacobianE(y);
    for (int k = 0; k < n; ++k) c_[k] = 2.0 * lambda_[k] * y[n + k];
  }
  // (s I - J) x = b with J = [[A, 0], [diag(c), A]]
  void SolveReal(double* b) const {
    const int n = K();
    SolveRealK(b);
    for (int k = 0; k < n; ++k) b[n + k] += c_[k] * b[k];
    SolveRealK(b + n);
  }
  void SolveComplex(std::complex<double>* b) const {
    const int n = K();
    SolveComplexK(b);
    for (int k = 0; k < n; ++k) b[n + k] += c_[k] * b[k];
    SolveComplexK(b + n);
  }

 private:
  std::vector<double> c_, qd_, tmpC_;
};

// (E, G): dG/dt = E
class EGProblem : public ExtinctionCore {
 public:
  using ExtinctionCore::ExtinctionCore;
  int Size() const { return 2 * K(); }
  void Rhs(const double* y, double* f) {
    const int n = K();
    RhsE(y, f);
    for (int k = 0; k < n; ++k) f[n + k] = y[k];
  }
  void Jacobian(const double* y) { JacobianE(y); }
  // (s I - J) x = b with J = [[A, 0], [I, 0]]
  void SolveReal(double* b) const {
    const int n = K();
    SolveRealK(b);
    for (int k = 0; k < n; ++k) b[n + k] = (b[n + k] + b[k]) / sr();
  }
  void SolveComplex(std::complex<double>* b) const {
    const int n = K();
    SolveComplexK(b);
    for (int k = 0; k < n; ++k) b[n + k] = (b[n + k] + b[k]) / sc();
  }
};

}  // namespace EvoTempo

#endif  // EVOTEMPO_DCET_ODE_H
