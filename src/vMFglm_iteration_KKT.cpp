#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::mat diag_matrix(const arma::mat& X, double lambda);
arma::mat b2_vMF(arma::vec theta);
arma::vec b1_vMF(arma::vec theta);
arma::vec subgrad(arma::vec theta);
double Hq_cpp(arma::vec theta);
double Bq_cpp(arma::vec theta);
double Cq_cpp(arma::vec theta, bool logarithm = true);
arma::mat calculate_b1_vMF(const arma::mat& Offset, const Rcpp::List& Xt_list, const arma::vec& beta, int n);
arma::mat calculate_b2_vMF(const arma::mat& Offset, const Rcpp::List& Xt_list, const arma::vec& beta, int n);

// [[Rcpp::export]]
List sphereGLM_iteration_KKT(arma::mat X, arma::mat Y, arma::mat Offset, arma::vec mu, arma::vec beta,
                             arma::mat Xt, List Xt_list, double eps, int maxit,
                             double lambda) {
  bool inv_success;
  
  int n = X.n_rows;
  int p = X.n_cols;
  int q = Y.n_cols;
  
  int l1 = 1;
  double crit = 1.0;
  
  // initialization
  vec bmu = mu;
  vec bbeta = beta;
  vec gamma = ones<vec>(p); // gamma 초기화 추가
  
  // container for saving results
  List beta_list;
  List Fn_list;
  std::vector<double> loglik_list;
  std::vector<double> crit_list;
  
  // double loglik_new;
  
  while (crit > eps && l1 <= maxit) {
    beta_list.push_back( join_cols(bmu, bbeta) );
    
    // b1, b2 계산
    List b1_list(n);
    List b2_list(n);
    
    
    for(int i = 0; i < n; i++) {
      mat Xi = as<mat>(Xt_list[i]);
      vec theta = Offset.row(i).t() + bmu + Xi.t() * bbeta;
      b1_list[i] = b1_vMF(theta);
      b2_list[i] = b2_vMF(theta);
    }
    
    
    // s1, s2, s3 calculation
    vec s1 = zeros<vec>(q);
    vec s2 = zeros<vec>(p * q);
    vec s3 = zeros<vec>(p);
    
    for(int i = 0; i < n; i++) {
      vec b1i = as<vec>(b1_list[i]);
      s1 += Y.row(i).t() - b1i;
      mat Xi = as<mat>(Xt_list[i]);
      s2 += Xi * (Y.row(i).t() - b1i);
    }

    
    mat beta_mat = reshape(bbeta, q, p).t();
    for(int j = 0; j < p; j++) {
      s1 -= gamma(j) * beta_mat.row(j).t();
      s3(j) = -dot(beta_mat.row(j), bmu);
    }
    
    s2 -= vectorise(kron(gamma, bmu));  // bmu로 수정
    
    
    vec sn = join_cols(s1, s2, s3);
    
    // Hessian matrix configuration
    mat H11 = zeros<mat>(q, q);
    mat H12 = zeros<mat>(q, p * q);
    mat H13 = -beta_mat.t();
    mat H22 = zeros<mat>(p * q, p * q);
    mat H23 = -kron(eye<mat>(p,p), bmu);
    mat H33 = zeros<mat>(p, p);
    
    
    for(int i = 0; i < n; i++) {
      mat b2i = as<mat>(b2_list[i]);
      mat Xi = as<mat>(Xt_list[i]);
      H11 -= b2i;
      H12 -= b2i * Xi.t();
      H22 -= Xi * b2i * Xi.t();
    }
    
    
    // gamma
    mat temp = zeros<mat>(q, p * q);
    for(int j = 0; j < p; j++) {
      temp.cols(j*q, (j+1)*q-1) = gamma(j) * eye<mat>(q, q);
    }
    
    H12 = H12 - temp;
    

    
    // Fn matrix
    int total_dim = (p + 1) * q + p;
    mat Fn = zeros<mat>(total_dim, total_dim);
    
    
      
    Fn.submat(0, 0, q-1, q-1) = -H11;
    Fn.submat(0, q, q-1, (p+1)*q-1) = -H12;
    Fn.submat(0, (p+1)*q, q-1, total_dim-1) = -H13;
    Fn.submat(q, q, (p+1)*q-1, (p+1)*q-1) = -H22;
    Fn.submat(q, (p+1)*q, (p+1)*q-1, total_dim-1) = -H23;
    Fn.submat((p+1)*q, (p+1)*q, total_dim-1, total_dim-1) = -H33;

    for(int i = 0; i < Fn.n_rows; i++) {
      for(int j = 0; j < i; j++) {
        Fn(i,j) = Fn(j,i);
      }
    }
    
    
    Fn_list.push_back(Fn);
    
    
    // update
    vec solutions;
    inv_success = true;
    try {
      solutions = join_cols(bmu, bbeta, gamma) + solve(Fn, sn);
    } catch(...) {
      solutions = join_cols(bmu, bbeta, gamma) + pinv(Fn) * sn;
      inv_success = false;
    }
    
    
    
    // updated parameter
    vec bmu_new = solutions.subvec(0, q-1);
    vec bbeta_new = solutions.subvec(q, (p+1)*q-1);
    vec gamma_new = solutions.subvec((p+1)*q, total_dim-1);
    
    // loglik calculation
    double loglik = 0;
    for (int i = 0; i < n; ++i) {
      mat Xi = as<mat>(Xt_list[i]);
      vec theta = Offset.row(i).t() + bmu_new + Xi.t() * bbeta_new;
      loglik += dot(theta, Y.row(i).t()) + Cq_cpp(theta, true);
    }
    
    // loglik_new = loglik;
    
    // check convergence
    crit = norm(join_cols(bmu_new - bmu, bbeta_new - bbeta), 2);
    // crit = abs(loglik_new - loglik_old);
    
    // save result
    loglik_list.push_back(loglik);
    crit_list.push_back(crit);
    
    // update parameter
    bmu = bmu_new;
    bbeta = bbeta_new;
    gamma = gamma_new;
    
    l1++;
  }
  
  return List::create(
    Named("beta") = join_cols(bmu, bbeta),
    Named("gamma") = gamma,
    Named("beta.list") = beta_list,
    Named("loglik.list") = loglik_list,
    Named("crit.list") = crit_list,
    Named("Fn.list") = Fn_list,
    Named("iterations") = l1 - 1,
    Named("inv.success") = inv_success
  );
}