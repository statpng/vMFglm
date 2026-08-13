// vMFglm_fast.cpp
// ============================================================
// 원본 vMFglm_iteration_FixedGamma의 drop-in 대체 (동일 시그니처 + 추가 옵션).
//
// 핵심 개선
//  (1) [메모리/속도] calculate_b2_vMF의 (nq x nq) 블록대각 W와
//      robust_inv_spd(W) 를 완전히 제거.
//      W_i = c1_i I_q + c2_i s_i s_i'  구조를 이용:
//        XWX             = (X1' diag(c1) X1) (x) I_q  +  M' diag(c2) M
//        W(W+la I)^{-1}  = d1 I + (dpar - d1) s s'    (Sherman-Morrison)
//      -> 반복당 O(n (pq)^2), 원본은 O((nq)^3).
//  (2) [정확성] method="newton" (기본): 벌점 Fisher scoring + step-halving.
//      고집중(kappa 큰) 데이터에서 IRLS의 s-방향 감쇠 h/(h+la)로 인한
//      집중도 효과 소멸 문제를 해결. method="irls"는 원본 알고리즘을
//      반복 단위까지 재현(정확 재현 검증용).
//  (3) [안정성] 업로드된 bessel.cpp의 점근 급수 log_Iv를 사용하되,
//      중간 kappa 구간은 boost 대신 R 내장 bessel_i(scaled)로 대체
//      (BH 의존성 제거, 임의 q 지원).
//
// R 래퍼 호환: Xt, Xt_list 인자는 받되 사용하지 않음 (기존 호출부 무수정).
// ============================================================

#include <RcppArmadillo.h>
#include <Rmath.h>
#include <cmath>
#include <limits>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

namespace {

// ---------- 업로드된 bessel.cpp의 점근 급수 (최적 절단) ----------
double bessel_asymp_series(double nu, double z) {
  const double mu4 = 4.0 * nu * nu;
  double sum = 1.0, term = 1.0, prev_abs = 1.0;
  for (int k = 1; k <= 200; ++k) {
    const double odd = 2.0 * k - 1.0;
    term *= (odd * odd - mu4) / (8.0 * k * z);
    const double a = std::abs(term);
    if (a > prev_abs) break;                       // optimal truncation
    sum += term;
    if (a < 1e-17 * std::abs(sum)) break;
    prev_abs = a;
  }
  return sum;
}

// log I_nu(kappa): 중간 구간은 R 내장(스케일드), 큰 kappa는 점근 급수
double log_Iv(double nu, double kappa) {
  if (kappa <= 0.0) {
    if (kappa == 0.0 && nu == 0.0) return 0.0;
    return -std::numeric_limits<double>::infinity();
  }
  const double thr = std::max(100.0, 4.0 * nu * nu + 50.0);
  if (kappa < thr) {
    double s = R::bessel_i(kappa, nu, 2.0);        // e^{-kappa} I_nu(kappa)
    if (std::isfinite(s) && s > 0.0) return std::log(s) + kappa;
  }
  const double F = bessel_asymp_series(nu, kappa);
  return kappa - 0.5 * std::log(2.0 * M_PI * kappa) + std::log(F);
}

// A_q(kappa) = I_{q/2} / I_{q/2-1}
double Aq_general(double kappa, int q) {
  if (kappa < 1e-12) return 0.0;
  return std::exp(log_Iv(0.5 * q, kappa) - log_Iv(0.5 * q - 1.0, kappa));
}

// 관측별 공통량: r, S(단위방향), A, c1 = A/r, c2 = A' - c1
struct Moments { vec r, A, c1, c2; mat S; };

Moments vmf_moments_cpp(const mat& Theta, int q) {
  const uword n = Theta.n_rows;
  Moments mo;
  mo.r  = sqrt(sum(square(Theta), 1));
  mo.S  = Theta;
  mo.A.set_size(n); mo.c1.set_size(n); mo.c2.set_size(n);
  const double tol = 1e-12;
  for (uword i = 0; i < n; ++i) {
    double ri = mo.r(i);
    if (ri < tol) {
      mo.S.row(i).zeros();
      mo.A(i) = 0.0; mo.c1(i) = 1.0 / q; mo.c2(i) = 0.0;
    } else {
      mo.S.row(i) /= ri;
      double A = Aq_general(ri, q);
      double a = A / ri;
      mo.A(i) = A; mo.c1(i) = a;
      mo.c2(i) = (1.0 - A * A - (q - 1.0) * a) - a;   // A'(r) - A/r
    }
  }
  return mo;
}

// XWX = (X1' diag(c1) X1) (x) I_q + M' diag(c2) M,  M[i,] = x1_i (x) s_i
mat build_XWX(const mat& X1, const Moments& mo, int q) {
  const uword n = X1.n_rows, p1 = X1.n_cols;
  mat G1 = X1.t() * (X1.each_col() % mo.c1);
  mat M(n, p1 * q);
  for (uword j = 0; j < p1; ++j)
    M.cols(j * q, j * q + q - 1) = mo.S.each_col() % X1.col(j);
  return kron(G1, eye(q, q)) + M.t() * (M.each_col() % mo.c2);
}

// 로그가능도 (스케일드 Bessel, 벡터화)
double loglik_cpp(const mat& Theta, const mat& Y, int q) {
  const uword n = Theta.n_rows;
  const double hq = 0.5 * q;
  double ll = accu(Theta % Y);
  for (uword i = 0; i < n; ++i) {
    double r = std::sqrt(dot(Theta.row(i), Theta.row(i)));
    double C = (r < 1e-12)
      ? (R::lgammafn(hq) - std::log(2.0) - hq * std::log(M_PI))
      : ((hq - 1.0) * std::log(r) - hq * std::log(2.0 * M_PI) - log_Iv(hq - 1.0, r));
    if (!std::isfinite(C)) return -datum::inf;
    ll += C;
  }
  return ll;
}

// 벌점 행렬 (orthogonal MMt + zero_beta), 원본과 동일
mat build_penalty(const vec& beta, int p, int q, bool orthogonal,
                  const vec& gamma, Nullable<IntegerVector> zero_beta) {
  const int d = (p + 1) * q;
  mat pen(d, d, fill::zeros);
  if (orthogonal) {
    vec mu = beta.subvec(0, q - 1);
    double nrm = norm(mu, 2);
    vec mu_dir = (nrm > 1e-12) ? vec(mu / nrm) : mu;
    mat MMt = mu_dir * mu_dir.t();
    for (int j = 1; j <= p; ++j)
      pen.submat(j * q, j * q, (j + 1) * q - 1, (j + 1) * q - 1) = gamma(j - 1) * MMt;
  }
  if (zero_beta.isNotNull()) {
    IntegerVector zb = as<IntegerVector>(zero_beta);
    for (int idx : zb)
      if (idx > 0 && idx <= p)
        pen.submat(idx * q, idx * q, (idx + 1) * q - 1, (idx + 1) * q - 1)
          += 1e+12 * eye(q, q);
  }
  return pen;
}

bool robust_solve(const mat& A, const vec& b, vec& x, double base_lambda) {
  if (solve(x, A, b) && x.is_finite()) return true;
  mat I = eye(A.n_rows, A.n_cols);
  double bump = std::max(1e-8, base_lambda);
  for (int k = 0; k < 6; ++k) {
    if (solve(x, A + bump * I, b) && x.is_finite()) return true;
    bump *= 10.0;
  }
  mat Ap;
  if (pinv(Ap, A)) { x = Ap * b; return x.is_finite(); }
  return false;
}

} // namespace


// [[Rcpp::export]]
List vMFglm_iteration_FixedGamma_fast(
    arma::mat X, arma::mat Y, arma::mat Offset, arma::vec beta,
    arma::mat Xt, List Xt_list,                       // 미사용 (호환용)
    double eps, int maxit, double lambda, bool orthogonal, arma::vec gamma,
    Rcpp::Nullable<Rcpp::IntegerVector> zero_beta = R_NilValue,
    bool verbose = false,
    std::string method = "newton",
    double kappa_max = 1e6) {

  const int n = X.n_rows, p = X.n_cols, q = Y.n_cols;
  const int d = (p + 1) * q;
  mat X1 = join_horiz(ones(n, 1), X);

  // ridge (newton 전용, 절편(mu) 제외): 원본 R IRLS의 bdiag(0, la I)와 동일 역할
  vec ridge_diag(d, fill::zeros);
  ridge_diag.subvec(q, d - 1).fill(lambda);

  auto Bmat = [&](const vec& b) { return mat(reshape(b, q, p + 1).t()); };

  List beta_list, Fn_list;
  std::vector<double> loglik_list, crit_list;
  std::string status = "ok";
  double crit = datum::inf;

  mat pen = build_penalty(beta, p, q, orthogonal, gamma, zero_beta);
  auto penll = [&](const vec& b) {
    double ll = loglik_cpp(X1 * Bmat(b) + Offset, Y, q);
    if (!std::isfinite(ll)) return -datum::inf;
    return ll - 0.5 * dot(b, pen * b)
              - 0.5 * lambda * dot(b.subvec(q, d - 1), b.subvec(q, d - 1));
  };

  double loglik_curr = loglik_cpp(X1 * Bmat(beta) + Offset, Y, q);
  double penll_curr  = penll(beta);

  int l1 = 1;
  while (crit > eps && l1 <= maxit) {
    beta_list.push_back(beta);
    vec beta_old = beta;

    mat Theta = X1 * Bmat(beta) + Offset;
    Moments mo = vmf_moments_cpp(Theta, q);
    mat XWX = build_XWX(X1, mo, q);
    XWX = 0.5 * (XWX + XWX.t());
    mat Resid = Y - (mo.S.each_col() % mo.A);         // y_i - b1(theta_i)
    pen = build_penalty(beta, p, q, orthogonal, gamma, zero_beta);

    vec beta_new;

    if (method == "irls") {
      // ── 원본 알고리즘 재현: W(W+la)^{-1} 필터 (Sherman-Morrison) ──
      vec d1   = mo.c1 / (mo.c1 + lambda);
      vec dpar = (mo.c1 + mo.c2) / (mo.c1 + mo.c2 + lambda);
      mat V = Resid.each_col() % d1
            + mo.S.each_col() % ((dpar - d1) % sum(mo.S % Resid, 1));
      vec rhs = XWX * beta + vectorise((X1.t() * V).t());
      mat Asys = 0.5 * ((XWX + pen) + (XWX + pen).t());
      vec beta_cand;
      if (!robust_solve(Asys, rhs, beta_cand, lambda)) { status = "solve failed"; break; }

      // 원본과 동일한 step-halving (비벌점 loglik 기준)
      vec dir = beta_cand - beta_old;
      double step = 1.0; vec beta_try = beta_cand;
      double ll_try = loglik_cpp(X1 * Bmat(beta_try) + Offset, Y, q);
      int h = 0;
      while ((!std::isfinite(ll_try) || ll_try < loglik_curr - 1e-10) && h < 20) {
        step *= 0.5; beta_try = beta_old + step * dir;
        ll_try = loglik_cpp(X1 * Bmat(beta_try) + Offset, Y, q); ++h;
      }
      if (!std::isfinite(ll_try)) { status = "step-halving failed"; break; }
      beta_new = beta_try; loglik_curr = ll_try;
      crit = norm(beta_old - beta_new, 2);            // 원본과 동일(절대 기준)

    } else {
      // ── Newton (벌점 Fisher scoring + step-halving) ──
      vec score = vectorise((X1.t() * Resid).t())
                - pen * beta - ridge_diag % beta;
      mat H = XWX + pen; H.diag() += ridge_diag;
      H = 0.5 * (H + H.t());
      vec step_v;
      if (!robust_solve(H, score, step_v, lambda)) { status = "solve failed"; break; }

      double fac = 1.0; vec beta_try; double pl_try;
      for (int h = 0; h < 40; ++h) {
        beta_try = beta_old + fac * step_v;
        pl_try = penll(beta_try);
        if (std::isfinite(pl_try) && pl_try >= penll_curr - 1e-10) break;
        fac *= 0.5;
      }
      if (!std::isfinite(pl_try)) { status = "step-halving failed"; break; }
      beta_new = beta_try; penll_curr = pl_try;
      loglik_curr = loglik_cpp(X1 * Bmat(beta_new) + Offset, Y, q);
      crit = norm(beta_old - beta_new, 2) / (1.0 + norm(beta_old, 2));  // 상대 기준

      if (norm(beta_new, 2) > kappa_max) {
        beta = beta_new;
        loglik_list.push_back(loglik_curr); crit_list.push_back(crit);
        Fn_list.push_back(XWX);
        status = "kappa_max reached (MLE may diverge)";
        Rcpp::warning("||beta|| > kappa_max: 데이터가 거의 무잡음이어서 MLE가 "
                      "발산할 수 있습니다. method='irls' 사용을 고려하세요.");
        ++l1; break;
      }
    }

    beta = beta_new;
    loglik_list.push_back(loglik_curr);
    crit_list.push_back(crit);
    Fn_list.push_back(XWX);
    if (verbose)
      Rcout << "[iter " << l1 << "] ll=" << loglik_curr
            << " crit=" << crit << " ||beta||=" << norm(beta, 2) << std::endl;
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