// vMFglm_iteration_FixedGamma.cpp
// Numerically robust IRLS for sphere GLM (vMF), suitable for rank-1 S^q models.
// Drop-in replacement: same signature as the original, plus optional `verbose`.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// ---- forward decls (assumed defined elsewhere in your package) -------------
arma::mat diag_matrix(const arma::mat& X, double lambda);
arma::mat b2_vMF(arma::vec theta);
arma::vec b1_vMF(arma::vec theta);
arma::vec subgrad(arma::vec theta);
double Hq_cpp(arma::vec theta);
double Bq_cpp(arma::vec theta);
double Cq_cpp(arma::vec theta, bool logarithm = true);
arma::mat calculate_b1_vMF(const arma::mat& Offset, const Rcpp::List& Xt_list,
                           const arma::vec& beta, int n);
arma::mat calculate_b2_vMF(const arma::mat& Offset, const Rcpp::List& Xt_list,
                           const arma::vec& beta, int n);
// ---------------------------------------------------------------------------

namespace {

inline bool is_finite_mat(const arma::mat& M) { return M.is_finite(); }
inline bool is_finite_vec(const arma::vec& v) { return v.is_finite(); }

// Robust SPD inverse: try inv_sympd, fall back to inv, then to ridge-bumped inv.
// Returns false if everything failed.
bool robust_inv_spd(const arma::mat& A, arma::mat& Ainv, double base_lambda) {
  if (arma::inv_sympd(Ainv, A))            return is_finite_mat(Ainv);
  if (arma::inv(Ainv, A))                  return is_finite_mat(Ainv);
  // Last-ditch: bump the diagonal until invertible (or give up after 6 tries).
  double bump = std::max(1e-8, base_lambda);
  arma::mat I = arma::eye(A.n_rows, A.n_cols);
  for (int k = 0; k < 6; ++k) {
    if (arma::inv_sympd(Ainv, A + bump * I) && is_finite_mat(Ainv)) return true;
    if (arma::inv(Ainv, A + bump * I)       && is_finite_mat(Ainv)) return true;
    bump *= 10.0;
  }
  return false;
}

// Robust linear solve. Tries plain solve, then ridge bump, then pinv.
bool robust_solve(const arma::mat& A, const arma::vec& b, arma::vec& x,
                  double base_lambda) {
  // No solve_opts::fast — we want LAPACK to flag rank deficiency.
  if (arma::solve(x, A, b) && is_finite_vec(x)) return true;

  arma::mat I = arma::eye(A.n_rows, A.n_cols);
  double bump = std::max(1e-8, base_lambda);
  for (int k = 0; k < 6; ++k) {
    if (arma::solve(x, A + bump * I, b) && is_finite_vec(x)) return true;
    bump *= 10.0;
  }
  // Pseudo-inverse fallback: stable for genuine rank deficiency.
  arma::mat Apinv;
  if (arma::pinv(Apinv, A)) {
    x = Apinv * b;
    return is_finite_vec(x);
  }
  return false;
}

// Compute log-likelihood for a candidate beta. Returns -inf on any NaN/Inf.
double compute_loglik(const arma::mat& Y, const arma::mat& Offset,
                      const Rcpp::List& Xt_list, const arma::vec& beta, int n) {
  double ll = 0.0;
  for (int i = 0; i < n; ++i) {
    arma::mat Xt_i = Rcpp::as<arma::mat>(Xt_list[i]);
    arma::vec theta = Offset.row(i).t() + Xt_i.t() * beta;
    if (!theta.is_finite()) return -arma::datum::inf;
    double C = Cq_cpp(theta, true);
    if (!std::isfinite(C))  return -arma::datum::inf;
    ll += arma::dot(theta, Y.row(i).t()) + C;
  }
  return ll;
}

} // anonymous namespace


// [[Rcpp::export]]
List vMFglm_iteration_FixedGamma(arma::mat X, arma::mat Y, arma::mat Offset,
                                    arma::vec beta,
                                    arma::mat Xt, List Xt_list,
                                    double eps, int maxit,
                                    double lambda, bool orthogonal,
                                    arma::vec gamma,
                                    Rcpp::Nullable<Rcpp::IntegerVector> zero_beta = R_NilValue,
                                    bool verbose = false) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  const int q = Y.n_cols;

  arma::vec beta_old = beta;
  arma::vec beta_new = beta;
  List beta_list, Fn_list;
  std::vector<double> loglik_list, crit_list;
  double crit = arma::datum::inf;
  std::string status = "ok";

  double loglik_curr = compute_loglik(Y, Offset, Xt_list, beta, n);
  if (!std::isfinite(loglik_curr)) {
    Rcpp::warning("Initial beta yields non-finite log-likelihood.");
  }

  int l1 = 1;
  while (crit > eps && l1 <= maxit) {
    beta_list.push_back(beta);
    beta_old = beta;

    // ----- Compute working response and weights -----
    arma::mat eta = calculate_b1_vMF(Offset, Xt_list, beta, n);
    arma::mat W   = calculate_b2_vMF(Offset, Xt_list, beta, n);

    if (!eta.is_finite() || !W.is_finite()) {
      status = "non-finite eta/W";
      if (verbose) Rcpp::Rcout << "[iter " << l1 << "] " << status << std::endl;
      break;
    }

    // Symmetrize W (b2 should be symmetric; tiny asymmetry from FP).
    W = 0.5 * (W + W.t());

    arma::mat W_reg = W + diag_matrix(W, lambda);
    arma::mat W_inv;
    if (!robust_inv_spd(W_reg, W_inv, lambda)) {
      status = "W inversion failed";
      if (verbose) Rcpp::Rcout << "[iter " << l1 << "] " << status << std::endl;
      break;
    }

    arma::vec yt = arma::vectorise(Y.t());
    arma::vec Z  = Xt.t() * beta + W_inv * (yt - arma::vectorise(eta.t()));
    if (!Z.is_finite()) {
      status = "non-finite Z";
      if (verbose) Rcpp::Rcout << "[iter " << l1 << "] " << status << std::endl;
      break;
    }

    // ----- Build penalty -----
    arma::vec mu = beta.subvec(0, q - 1);
    arma::mat XWX     = Xt * W * Xt.t();
    XWX = 0.5 * (XWX + XWX.t());
    arma::mat penalty = arma::zeros(XWX.n_rows, XWX.n_cols);

    if (orthogonal) {
      // Identifiability for rank-1: mu acts as a unit direction.
      // Normalize mu when building the orthogonal penalty so rank-1 doesn't
      // smuggle scale into gamma.
      double mu_nrm = arma::norm(mu, 2);
      arma::vec mu_dir = (mu_nrm > 1e-12) ? (mu / mu_nrm) : mu;
      arma::mat MMt = mu_dir * mu_dir.t();
      for (int i = 1; i <= p; ++i) {
        penalty.submat(i * q, i * q, (i + 1) * q - 1, (i + 1) * q - 1)
            = gamma(i - 1) * MMt;
      }
    }
    if (zero_beta.isNotNull()) {
      Rcpp::IntegerVector zero_indices = Rcpp::as<Rcpp::IntegerVector>(zero_beta);
      for (int idx : zero_indices) {
        if (idx > 0 && idx <= p) {
          penalty.submat(idx * q, idx * q, (idx + 1) * q - 1, (idx + 1) * q - 1)
              += 1e+12 * arma::eye(q, q);
        }
      }
    }

    // ----- Solve for beta candidate -----
    arma::mat A = XWX + penalty;
    A = 0.5 * (A + A.t());
    arma::vec rhs = Xt * W * Z;

    if (verbose) {
      double rc = arma::rcond(A);
      Rcpp::Rcout << "[iter " << l1 << "] rcond(XWX+pen)=" << rc
                  << "  ||beta||=" << arma::norm(beta, 2) << std::endl;
    }

    arma::vec beta_cand;
    if (!robust_solve(A, rhs, beta_cand, lambda)) {
      status = "solve failed";
      if (verbose) Rcpp::Rcout << "[iter " << l1 << "] " << status << std::endl;
      break;
    }

    // ----- Step-halving: accept only if loglik does not get worse -----
    arma::vec direction = beta_cand - beta_old;
    double step = 1.0;
    arma::vec beta_try = beta_cand;
    double loglik_try = compute_loglik(Y, Offset, Xt_list, beta_try, n);

    int halvings = 0;
    const int max_halvings = 20;
    while ((!std::isfinite(loglik_try) || loglik_try < loglik_curr - 1e-10)
           && halvings < max_halvings) {
      step *= 0.5;
      beta_try   = beta_old + step * direction;
      loglik_try = compute_loglik(Y, Offset, Xt_list, beta_try, n);
      ++halvings;
    }

    if (!std::isfinite(loglik_try)) {
      status = "step-halving failed";
      if (verbose) Rcpp::Rcout << "[iter " << l1 << "] " << status << std::endl;
      // Keep last good beta and stop.
      break;
    }

    if (verbose && halvings > 0) {
      Rcpp::Rcout << "[iter " << l1 << "] step-halved " << halvings
                  << " times, step=" << step << std::endl;
    }

    beta        = beta_try;
    beta_new    = beta;
    loglik_curr = loglik_try;

    crit = arma::norm(beta_old - beta_new, 2);
    loglik_list.push_back(loglik_curr);
    crit_list.push_back(crit);
    Fn_list.push_back(XWX);
    ++l1;
  }

  return List::create(
    Named("beta")        = beta,
    Named("beta.list")   = beta_list,
    Named("loglik.list") = loglik_list,
    Named("crit.list")   = crit_list,
    Named("Fn.list")     = Fn_list,
    Named("iterations")  = l1 - 1,
    Named("gamma")       = gamma,
    Named("status")      = status
  );
}