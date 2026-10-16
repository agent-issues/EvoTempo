#ifndef EVOTEMPO_BAND_H
#define EVOTEMPO_BAND_H

// Banded LU factorization with partial pivoting (the algorithm of LAPACK's
// dgbtf2 / dgbtrs), for real or complex scalars.  A dense matrix is the
// special case kl = ku = n - 1.

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

namespace EvoTempo {

inline double Magnitude(double x) { return std::fabs(x); }
inline double Magnitude(const std::complex<double>& x) {
  return std::fabs(x.real()) + std::fabs(x.imag());
}

template <class T>
class BandLU {
 public:
  BandLU() : n_(0), kl_(0), ku_(0), ld_(1), singular_(false) {}

  void Resize(int n, int kl, int ku) {
    n_ = n;
    kl_ = std::max(0, std::min(kl, n - 1));
    ku_ = std::max(0, std::min(ku, n - 1));
    ld_ = 2 * kl_ + ku_ + 1;
    ab_.assign(static_cast<std::size_t>(ld_) * n_, T(0));
    piv_.assign(n_, 0);
  }

  int n() const { return n_; }
  int kl() const { return kl_; }
  int ku() const { return ku_; }

  void Zero() { std::fill(ab_.begin(), ab_.end(), T(0)); }

  // Element (i, j) of the matrix before factorization; needs
  // j - ku <= i <= j + kl.
  T& operator()(int i, int j) { return ab_[Index(i, j)]; }

  // Factorize in place; returns false if a pivot is exactly zero.
  bool Factor() {
    singular_ = false;
    int ju = 0;
    for (int j = 0; j < n_; ++j) {
      const int km = std::min(kl_, n_ - 1 - j);
      int p = j;
      double best = Magnitude(ab_[Index(j, j)]);
      for (int i = j + 1; i <= j + km; ++i) {
        const double m = Magnitude(ab_[Index(i, j)]);
        if (m > best) {
          best = m;
          p = i;
        }
      }
      piv_[j] = p;
      if (best == 0.0) {
        singular_ = true;
        continue;
      }
      ju = std::max(ju, std::min(j + ku_ + p - j, n_ - 1));
      if (p != j) {
        for (int c = j; c <= ju; ++c) {
          std::swap(ab_[Index(p, c)], ab_[Index(j, c)]);
        }
      }
      const T inv = T(1) / ab_[Index(j, j)];
      for (int i = j + 1; i <= j + km; ++i) ab_[Index(i, j)] *= inv;
      for (int c = j + 1; c <= ju; ++c) {
        const T u = ab_[Index(j, c)];
        if (u == T(0)) continue;
        for (int i = j + 1; i <= j + km; ++i) {
          ab_[Index(i, c)] -= ab_[Index(i, j)] * u;
        }
      }
    }
    return !singular_;
  }

  bool singular() const { return singular_; }

  // Solve A x = b in place, after Factor().
  void Solve(T* b) const {
    const int kv = kl_ + ku_;
    for (int j = 0; j < n_ - 1; ++j) {
      const int lm = std::min(kl_, n_ - 1 - j);
      const int p = piv_[j];
      if (p != j) std::swap(b[p], b[j]);
      const T bj = b[j];
      if (bj == T(0)) continue;
      for (int i = j + 1; i <= j + lm; ++i) b[i] -= ab_[Index(i, j)] * bj;
    }
    for (int j = n_ - 1; j >= 0; --j) {
      b[j] /= ab_[Index(j, j)];
      const T bj = b[j];
      if (bj == T(0)) continue;
      for (int i = std::max(0, j - kv); i < j; ++i) b[i] -= ab_[Index(i, j)] * bj;
    }
  }

 private:
  // Row kl + ku + i - j of column j in LAPACK band storage
  std::size_t Index(int i, int j) const {
    return static_cast<std::size_t>(kl_ + ku_ + i - j) +
      static_cast<std::size_t>(j) * ld_;
  }

  int n_, kl_, ku_, ld_;
  bool singular_;
  std::vector<T> ab_;
  std::vector<int> piv_;
};

}  // namespace EvoTempo

#endif  // EVOTEMPO_BAND_H
