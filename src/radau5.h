#ifndef EVOTEMPO_RADAU5_H
#define EVOTEMPO_RADAU5_H

// Three-stage Radau IIA (order 5) integrator for stiff systems, after
// RADAU5 of Hairer & Wanner (1996, Solving Ordinary Differential Equations
// II, section IV.8): simplified Newton iterations on the transformed stage
// equations (one real and one complex linear system per step), the embedded
// error estimate filtered by (I - h gamma0 J)^-1, and the same step-size
// control and tolerance transformation.
//
// A problem supplies
//   int Size() const;
//   void Rhs(const double* y, double* f);
//   void Jacobian(const double* y);          // J at y, kept by the problem
//   bool Factor(double sr, std::complex<double> sc);
//                                            // LUs of sr I - J and sc I - J
//   void SolveReal(double* b) const;         // b <- (sr I - J)^-1 b
//   void SolveComplex(std::complex<double>* b) const;

#include <algorithm>
#include <cmath>
#include <complex>
#include <initializer_list>
#include <limits>
#include <vector>

namespace EvoTempo {

// Steps and rejections over all Radau5 integrations, for profiling
struct SolverCounter {
  long steps = 0, rejects = 0, newton = 0, rhs = 0;
};
inline SolverCounter& Counter() {
  static SolverCounter counter;
  return counter;
}

namespace radau {
// Collocation constants of Radau IIA, s = 3, computed to 40 digits
const double kGamma = 3.6378342527444957322;  // real eigenvalue of A^-1
const double kAlpha = 2.6810828736277521339;  // complex pair alpha +- i beta
const double kBeta = 3.0504301992474105694;
const double kDD1 = -10.048809399827415562;  // -(13 + 7 sqrt 6) / 3
const double kDD2 = 1.3821427331607488958;   // (-13 + 7 sqrt 6) / 3
const double kDD3 = -0.33333333333333333333;
// T and T^-1 with T^-1 A^-1 T = [[gamma, 0, 0], [0, alpha, beta],
// [0, -beta, alpha]]
const double kT[3][3] = {
  {0.099580415535919304075, 0.024647066019156465398, -0.15697287755978402247},
  {0.26383580323107993402, 0.31334462687440969388, 0.36027433855335156149},
  {1.0544443077337475481, -0.39296916392505438235, 1.0277301834059091024}};
const double kTI[3][3] = {
  {3.9629580821892511305, 0.31076351624992975469, 0.49635285776667563311},
  {0.92948941545931617518, 2.289684736920919586, -0.66068899877361969815},
  {-3.7105633119033522723, 0.55665648941319216267, 0.21113824867815223251}};
}  // namespace radau

class Radau5 {
 public:
  // rtol, atol as given to an ODE solver; transformed as in RADAU5
  Radau5(double rtol, double atol, double h0 = 1e-6, long maxSteps = 1000000)
    : h_(h0), faccon_(1.0), maxSteps_(maxSteps), nStep_(0), nReject_(0) {
    const double ratio = atol / rtol;
    rtol_ = 0.1 * std::pow(rtol, 2.0 / 3.0);
    atol_ = rtol_ * ratio;
    const double uround = std::numeric_limits<double>::epsilon();
    fnewt_ = std::max(10.0 * uround / rtol_, std::min(0.03, std::sqrt(rtol_)));
  }

  long steps() const { return nStep_; }
  long rejects() const { return nReject_; }

  // Advance y over [0, span]; returns false on failure.
  template <class P>
  bool Integrate(P& prob, double* y, double span) {
    if (!(span > 0)) return true;
    const int n = prob.Size();
    Allocate(n);
    const double uround = std::numeric_limits<double>::epsilon();
    const int nitMax = 7;
    double t = 0.0;
    double h = std::min(h_, span);
    bool first = true, reject = false, last = false;
    long nLocal = 0;
    while (true) {
      if (t + 1.0001 * h >= span) {
        h = span - t;
        last = true;
      }
      if (!(h > 0)) break;
      prob.Jacobian(y);
      prob.Rhs(y, f0_.data());
      for (int i = 0; i < n; ++i) scal_[i] = atol_ + rtol_ * std::fabs(y[i]);
      while (true) {
        if (++nLocal > maxSteps_ || h < 1e-14 * span) return false;
        const double sr = radau::kGamma / h;
        const std::complex<double> sc(radau::kAlpha / h, -radau::kBeta / h);
        if (!prob.Factor(sr, sc)) {
          h *= 0.5;
          reject = true;
          last = false;
          continue;
        }
        std::fill(z1_.begin(), z1_.end(), 0.0);
        std::fill(z2_.begin(), z2_.end(), 0.0);
        std::fill(z3_.begin(), z3_.end(), 0.0);
        std::fill(w1_.begin(), w1_.end(), 0.0);
        std::fill(w2_.begin(), w2_.end(), 0.0);
        std::fill(w3_.begin(), w3_.end(), 0.0);
        faccon_ = std::pow(std::max(faccon_, uround), 0.8);
        double dynOld = 0.0, thqOld = 0.0, theta = 0.0;
        bool converged = false;
        int it = 0;
        for (; it < nitMax; ++it) {
          ++Counter().newton;
          Counter().rhs += 3;
          for (int i = 0; i < n; ++i) tmp_[i] = y[i] + z1_[i];
          prob.Rhs(tmp_.data(), f1_.data());
          for (int i = 0; i < n; ++i) tmp_[i] = y[i] + z2_[i];
          prob.Rhs(tmp_.data(), f2_.data());
          for (int i = 0; i < n; ++i) tmp_[i] = y[i] + z3_[i];
          prob.Rhs(tmp_.data(), f3_.data());
          for (int i = 0; i < n; ++i) {
            const double a1 = radau::kTI[0][0] * f1_[i] + radau::kTI[0][1] * f2_[i] +
              radau::kTI[0][2] * f3_[i];
            const double a2 = radau::kTI[1][0] * f1_[i] + radau::kTI[1][1] * f2_[i] +
              radau::kTI[1][2] * f3_[i];
            const double a3 = radau::kTI[2][0] * f1_[i] + radau::kTI[2][1] * f2_[i] +
              radau::kTI[2][2] * f3_[i];
            r1_[i] = a1 - sr * w1_[i];
            rc_[i] = std::complex<double>(a2, a3) -
              sc * std::complex<double>(w2_[i], w3_[i]);
          }
          prob.SolveReal(r1_.data());
          prob.SolveComplex(rc_.data());
          double sum = 0.0;
          for (int i = 0; i < n; ++i) {
            const double d1 = r1_[i] / scal_[i];
            const double d2 = rc_[i].real() / scal_[i];
            const double d3 = rc_[i].imag() / scal_[i];
            sum += d1 * d1 + d2 * d2 + d3 * d3;
            w1_[i] += r1_[i];
            w2_[i] += rc_[i].real();
            w3_[i] += rc_[i].imag();
          }
          const double dyno = std::sqrt(sum / (3.0 * n));
          for (int i = 0; i < n; ++i) {
            z1_[i] = radau::kT[0][0] * w1_[i] + radau::kT[0][1] * w2_[i] +
              radau::kT[0][2] * w3_[i];
            z2_[i] = radau::kT[1][0] * w1_[i] + radau::kT[1][1] * w2_[i] +
              radau::kT[1][2] * w3_[i];
            z3_[i] = radau::kT[2][0] * w1_[i] + radau::kT[2][1] * w2_[i] +
              radau::kT[2][2] * w3_[i];
          }
          if (!std::isfinite(dyno)) break;
          if (it > 0) {
            const double thq = dyno / dynOld;
            theta = it == 1 ? thq : std::sqrt(thq * thqOld);
            thqOld = thq;
            if (theta >= 0.99) break;
            faccon_ = theta / (1.0 - theta);
          }
          dynOld = std::max(dyno, uround);
          if (faccon_ * dyno <= fnewt_) {
            converged = true;
            break;
          }
        }
        if (!converged) {
          h *= 0.5;
          reject = true;
          last = false;
          faccon_ = 1.0;
          continue;
        }
        // Error estimate
        for (int i = 0; i < n; ++i) {
          ynew_[i] = y[i] + z3_[i];
          fz_[i] = (radau::kDD1 * z1_[i] + radau::kDD2 * z2_[i] +
                    radau::kDD3 * z3_[i]) / h;
          err_[i] = f0_[i] + fz_[i];
          scalE_[i] = atol_ + rtol_ * std::max(std::fabs(y[i]), std::fabs(ynew_[i]));
        }
        prob.SolveReal(err_.data());
        double err = Norm(n);
        if (err >= 1.0 && (first || reject)) {
          for (int i = 0; i < n; ++i) tmp_[i] = y[i] + err_[i];
          prob.Rhs(tmp_.data(), f1_.data());
          for (int i = 0; i < n; ++i) err_[i] = f1_[i] + fz_[i];
          prob.SolveReal(err_.data());
          err = Norm(n);
        }
        err = std::max(err, 1e-10);
        if (!std::isfinite(err)) err = 1e10;
        const double fac = std::min(0.9, 0.9 * (2.0 * nitMax + 1.0) /
                                       (2.0 * nitMax + it + 1.0));
        const double quot = std::max(0.125, std::min(5.0, std::pow(err, 0.25) / fac));
        double hNew = h / quot;
        if (err < 1.0) {
          ++nStep_;
          ++Counter().steps;
          t += h;
          for (int i = 0; i < n; ++i) y[i] = ynew_[i];
          first = false;
          if (reject) hNew = std::min(hNew, h);
          reject = false;
          if (last) {
            h_ = std::max(hNew, h_ * 1e-3);
            return true;
          }
          h = hNew;
          break;
        }
        ++nReject_;
        ++Counter().rejects;
        reject = true;
        last = false;
        h = first ? 0.1 * h : hNew;
      }
    }
    return true;
  }

 private:
  void Allocate(int n) {
    if (static_cast<int>(f0_.size()) == n) return;
    for (std::vector<double>* v : {&f0_, &f1_, &f2_, &f3_, &z1_, &z2_, &z3_, &w1_,
                                   &w2_, &w3_, &r1_, &tmp_, &ynew_, &err_, &fz_,
                                   &scal_, &scalE_}) {
      v->assign(n, 0.0);
    }
    rc_.assign(n, std::complex<double>(0.0, 0.0));
  }

  double Norm(int n) const {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
      const double d = err_[i] / scalE_[i];
      sum += d * d;
    }
    return std::sqrt(sum / n);
  }

  double rtol_, atol_, fnewt_, h_, faccon_;
  long maxSteps_, nStep_, nReject_;
  std::vector<double> f0_, f1_, f2_, f3_, z1_, z2_, z3_, w1_, w2_, w3_, r1_, tmp_,
    ynew_, err_, fz_, scal_, scalE_;
  std::vector<std::complex<double> > rc_;
};

}  // namespace EvoTempo

#endif  // EVOTEMPO_RADAU5_H
