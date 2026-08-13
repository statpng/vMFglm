#include <RcppArmadillo.h>
#include <cmath>


// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Optimized BesselI function
Rcpp::NumericVector besselI_scaled_opt(Rcpp::NumericVector x, double nu, int k_max = 5) {
  int n = x.size();
  Rcpp::NumericVector result(n);
  
  double pi = 3.141592653589793238462643383280;
  
  for(int i = 0; i < n; i++) {
    double z = x[i] / nu;
    double sz = std::sqrt(1 + z*z);
    double t = 1 / sz;
    double eta = 1 / (sz + std::abs(z)) + std::log(z / (1 + sz));
    
    double d = 0;
    if (k_max >= 1) {
      double t2 = t * t;
      double u1_t = (t * (3 - 5 * t2)) / 24;
      
      if (k_max >= 2) {
        double u2_t = t2 * (81 + t2 * (-462 + t2 * 385)) / 1152;
        
        if (k_max >= 3) {
          double u3_t = t * t2 * (30375 + t2 * (-369603 + t2 * (765765 - t2 * 425425))) / 414720;
          
          if (k_max >= 4) {
            double t4 = t2 * t2;
            double u4_t = t4 * (4465125 + t2 * (-94121676 + t2 * (349922430 + t2 * (-446185740 + t2 * 185910725)))) / 39813120;
            
            if (k_max == 5) {
              double u5_t = t * t4 * (1519035525 + t2 * (-49286948607 + t2 * (284499769554 + t2 * (-614135872350 + t2 * (566098157625 - t2 * 188699385875))))) / 6688604160;
              d = (u1_t + (u2_t + (u3_t + (u4_t + u5_t/nu)/nu)/nu)/nu)/nu;
            } else {
              d = (u1_t + (u2_t + (u3_t + u4_t/nu)/nu)/nu)/nu;
            }
          } else {
            d = (u1_t + (u2_t + u3_t/nu)/nu)/nu;
          }
        } else {
          d = (u1_t + u2_t/nu)/nu;
        }
      } else {
        d = u1_t/nu;
      }
    }
    
    result[i] = (1 + d) * std::exp(nu * eta) / std::sqrt(2 * pi * nu * sz);
  }
  
  return result;
}




// old version
// double Cq_cpp(arma::vec theta, bool logarithm = true) {
//   int q = theta.n_elem;
//   double NORM = norm(theta);
//   
//   double result;
//   if (logarithm) {
//     Rcpp::NumericVector bessel_result = besselI_scaled_opt(Rcpp::NumericVector::create(NORM), q/2.0 - 1);
//     result = (q/2.0 - 1.0) * log(NORM) - ((q/2.0) * log(2*M_PI) +
//       log(bessel_result[0]) +
//       NORM);
//   } else {
//     Rcpp::NumericVector bessel_result = besselI_scaled_opt(Rcpp::NumericVector::create(NORM), q/2.0 - 1);
//     result = exp((q/2.0 - 1.0) * log(NORM) - ((q/2.0) * log(2*M_PI) +
//       log(bessel_result[0]) +
//       NORM));
//   }
//   
//   return result;
// }













// [[Rcpp::export]]
double Bq_cpp(arma::vec theta) {
  int q = theta.n_elem;
  double NORM = norm(theta);
  
  Rcpp::NumericVector bessel_result1 = besselI_scaled_opt(Rcpp::NumericVector::create(NORM), q/2.0);
  Rcpp::NumericVector bessel_result2 = besselI_scaled_opt(Rcpp::NumericVector::create(NORM), q/2.0 - 1);
  
  return bessel_result1[0] / bessel_result2[0];
}

// [[Rcpp::export]]
double Hq_cpp(arma::vec theta) {
  int q = theta.n_elem;
  double NORM = norm(theta);
  
  Rcpp::NumericVector bessel_result1 = besselI_scaled_opt(Rcpp::NumericVector::create(NORM), q/2.0);
  Rcpp::NumericVector bessel_result2 = besselI_scaled_opt(Rcpp::NumericVector::create(NORM), q/2.0 - 1);
  
  double I0 = bessel_result2[0];
  double I1 = bessel_result1[0];
  
  return 1 - pow(I1/I0, 2) - (q-1)/NORM * I1/I0;
}










// bessel_stable_general.cpp
// Numerically stable vMF helpers for *general* q >= 2:
//   Cq_cpp  — log normalizing constant
//   b1_vMF  — gradient of log partition (mean direction E[Y])
//   b2_vMF  — Hessian of log partition (Var[Y])
//
// Strategy:
//   • κ < 100         : direct boost::math::cyl_bessel_i (no overflow)
//   • κ ≥ 100         : asymptotic expansion in log scale, optimal truncation
//
// Drop-in replacement for the corresponding functions in src/bessel.cpp.

#include <RcppArmadillo.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <cmath>
#include <limits>
// [[Rcpp::depends(RcppArmadillo, BH)]]

namespace {

// -----------------------------------------------------------------
// Asymptotic series F(z;ν) := I_ν(z) · √(2πz) · e^{−z}
//   F(z;ν) = Σ_{k≥0} (-1)^k a_k(ν) / z^k
//   a_0 = 1,    a_k = a_{k-1} · (4ν² − (2k−1)²) / (8 k)
//
// Recurrence on the term itself:
//   t_0 = 1
//   t_k / t_{k-1} = ((2k-1)² − 4ν²) / (8 k z)
//
// Asymptotic series — diverges eventually; we stop at optimal truncation
// (when |t_k| stops decreasing).
// -----------------------------------------------------------------
double bessel_asymp_series(double nu, double z) {
  const double mu = 4.0 * nu * nu;
  double sum  = 1.0;
  double term = 1.0;
  double prev_abs = 1.0;
  const int max_terms = 200;
  for (int k = 1; k <= max_terms; ++k) {
    const double odd   = 2.0 * k - 1.0;
    const double ratio = (odd * odd - mu) / (8.0 * k * z);
    term *= ratio;
    const double aterm = std::abs(term);
    if (aterm > prev_abs) break;                       // optimal truncation
    sum += term;
    if (aterm < 1e-17 * std::abs(sum)) break;          // converged
    prev_abs = aterm;
  }
  return sum;
}

// -----------------------------------------------------------------
// log I_ν(κ), stable for ν ≥ 0 and κ ≥ 0.
// -----------------------------------------------------------------
double log_Iv(double nu, double kappa) {
  if (kappa < 0.0) Rcpp::stop("log_Iv: negative kappa");
  if (kappa == 0.0) {
    if (nu == 0.0) return 0.0;                          // I_0(0) = 1
    return -std::numeric_limits<double>::infinity();    // I_ν(0) = 0, ν > 0
  }

  // Threshold: boost overflows around κ ≈ 700; asymptotic needs κ ≫ ν².
  // Choose direct boost when κ is small enough not to overflow AND
  // not large enough that asymptotic is clearly superior.
  const double asymp_threshold = std::max(100.0, 4.0 * nu * nu + 50.0);
  if (kappa < asymp_threshold) {
    try {
      double Iv = boost::math::cyl_bessel_i(nu, kappa);
      if (std::isfinite(Iv) && Iv > 0.0) return std::log(Iv);
    } catch (const std::exception&) { /* fall through to asymptotic */ }
  }

  // Asymptotic: log I_ν(κ) = κ − ½ log(2πκ) + log F(κ; ν)
  const double F = bessel_asymp_series(nu, kappa);
  if (F <= 0.0 || !std::isfinite(F))
    Rcpp::stop("log_Iv: asymptotic failed at nu=%g, kappa=%g", nu, kappa);
  return kappa - 0.5 * std::log(2.0 * M_PI * kappa) + std::log(F);
}

// -----------------------------------------------------------------
// A_q(κ) = I_{q/2}(κ) / I_{q/2-1}(κ), stable for any q ≥ 2.
// Computed as exp(log_I_top − log_I_bot); the κ and ½log(2πκ) terms
// cancel exactly when both use the asymptotic branch.
// -----------------------------------------------------------------
double Aq_general(double kappa, int q) {
  if (kappa < 1e-12) return 0.0;
  const double nu_top = 0.5 * q;
  const double nu_bot = nu_top - 1.0;
  return std::exp(log_Iv(nu_top, kappa) - log_Iv(nu_bot, kappa));
}

// -----------------------------------------------------------------
// dA_q/dκ. Uses identity  A_q'(κ) = 1 − A_q² − (q−1)·A_q/κ.
// Accepts mild cancellation at very large κ (result O(1/κ²));
// finite for any κ.
// -----------------------------------------------------------------
double Aq_prime_general(double kappa, int q) {
  if (kappa < 1e-12) return 1.0 / static_cast<double>(q);
  const double A = Aq_general(kappa, q);
  return 1.0 - A * A - (q - 1.0) * A / kappa;
}

} // anonymous namespace


// =================================================================
// Public, exported. Replace originals in src/bessel.cpp.
// =================================================================

// [[Rcpp::export]]
double Cq_cpp(arma::vec theta, bool logarithm = true) {
  const int    q     = theta.n_elem;
  const double kappa = arma::norm(theta, 2);
  if (q < 2) Rcpp::stop("Cq_cpp: q must be >= 2 (got %d)", q);

  // log c_q(κ) = (q/2 − 1) log κ − (q/2) log(2π) − log I_{q/2-1}(κ)
  const double half_q = 0.5 * q;
  double logC;

  if (kappa < 1e-12) {
    // Limit: log c_q(0) = log Γ(q/2) − log 2 − (q/2) log π
    logC = boost::math::lgamma(half_q) - std::log(2.0) - half_q * std::log(M_PI);
  } else {
    logC = (half_q - 1.0) * std::log(kappa)
         - half_q * std::log(2.0 * M_PI)
         - log_Iv(half_q - 1.0, kappa);
  }

  if (!std::isfinite(logC))
    Rcpp::stop("Cq_cpp: non-finite at kappa=%g, q=%d", kappa, q);
  return logarithm ? logC : std::exp(logC);
}

// [[Rcpp::export]]
arma::vec b1_vMF(arma::vec theta) {
  const int    q     = theta.n_elem;
  const double kappa = arma::norm(theta, 2);
  if (q < 2) Rcpp::stop("b1_vMF: q must be >= 2 (got %d)", q);
  if (kappa < 1e-12) return arma::zeros<arma::vec>(q);

  const double A = Aq_general(kappa, q);
  if (!std::isfinite(A))
    Rcpp::stop("b1_vMF: non-finite Aq at kappa=%g, q=%d", kappa, q);
  return (A / kappa) * theta;          // = A_q(κ) · μ,  μ = θ/κ
}

// [[Rcpp::export]]
arma::mat b2_vMF(arma::vec theta) {
  const int    q     = theta.n_elem;
  const double kappa = arma::norm(theta, 2);
  if (q < 2) Rcpp::stop("b2_vMF: q must be >= 2 (got %d)", q);

  arma::mat Iq = arma::eye(q, q);
  if (kappa < 1e-12) return (1.0 / q) * Iq;            // uniform vMF limit

  const arma::vec mu  = theta / kappa;
  const arma::mat MMt = mu * mu.t();
  const double A  = Aq_general(kappa, q);
  const double Ap = Aq_prime_general(kappa, q);

  if (!std::isfinite(A) || !std::isfinite(Ap))
    Rcpp::stop("b2_vMF: non-finite at kappa=%g, q=%d", kappa, q);

  // Var(Y|θ) = A_q'(κ)·μμ' + (A_q(κ)/κ)·(I − μμ')
  return Ap * MMt + (A / kappa) * (Iq - MMt);
}



// [[Rcpp::export]]
arma::mat diag_matrix(const arma::mat& X, double lambda) {
  return lambda * arma::eye(X.n_rows, X.n_cols);
}






// [[Rcpp::export]]
arma::mat calculate_b1_vMF(const arma::mat& Offset, const Rcpp::List& Xt_list, const arma::vec& beta, int n) {
  int q = Offset.n_cols;
  arma::mat eta(n, q);
  
  for (int i = 0; i < n; ++i) {
    arma::mat Xt_i = Rcpp::as<arma::mat>(Xt_list[i]);
    arma::vec theta = Offset.row(i).t() + Xt_i.t() * beta;
    
    // b1_vMF 계산
    eta.row(i) = b1_vMF(theta).t();
  }
  
  return eta;
}



// [[Rcpp::export]]
arma::mat calculate_b2_vMF(const arma::mat& Offset, const Rcpp::List& Xt_list, const arma::vec& beta, int n) {
  int q = Offset.n_cols;
  arma::mat W(n * q, n * q, arma::fill::zeros);
  
  for (int i = 0; i < n; ++i) {
    arma::mat Xt_i = Rcpp::as<arma::mat>(Xt_list[i]);
    arma::vec theta = Offset.row(i).t() + Xt_i.t() * beta;
    
    // b2_vMF 계산 및 W 행렬에 할당
    W.submat(i*q, i*q, (i+1)*q-1, (i+1)*q-1) = b2_vMF(theta);
  }
  
  return W;
}