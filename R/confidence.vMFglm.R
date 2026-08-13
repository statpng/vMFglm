# ============================================================
#  confidence_vMFglm.R
#
#  변경점:
#  - confidence.vMFglm() 내에서 cov.mu / cov.ortho 사용 방식을
#    새 Jacobian 기반 공분산으로 업데이트.
#  - inference_beta_j(): 단독 호출 가능한 j번째 공변량 추론 함수 추가.
#    → 논문 Proposition 1, Corollary 1, Theorem (신뢰 영역) 구현.
# ============================================================


# ------------------------------------------------------------------
#  inference_beta_j
#
#  단독 호출 가능한 함수 — 모형 적합 후 j번째 공변량에 대한
#  위치/집중도 효과 분리 추론 전체를 수행.
#
#  Inputs:
#    fit   : vMFglm 객체
#    j     : 공변량 인덱스 (1-based)
#    alpha : 유의 수준 (default 0.05)
#
#  Returns (list):
#    mu_hat, beta_j_hat, s0
#    beta_loc   = P_{mu_perp} beta_j_hat  (위치 효과)
#    beta_conc  = P_{mu}      beta_j_hat  (집중도 효과)
#    Omega_loc_n, Omega_conc_n            (점근 공분산 / n)
#    W_loc, p_loc   (위치 효과 검정, chi^2_{q-1})
#    W_conc, p_conc (집중도 효과 검정, chi^2_1)
#    ci_conc        (집중도 스칼라 신뢰구간)
#    C_loc_centre, C_loc_shape, C_loc_cv  (위치 신뢰 타원 정보)
# ------------------------------------------------------------------
#' @export inference_beta_j
inference_beta_j <- function(fit, j, alpha = 0.05) {

  X  <- fit$X;  Y <- fit$Y
  p  <- ncol(X); q <- ncol(Y); n <- nrow(Y)
  mu        <- fit$mu
  beta_j    <- fit$beta[j + 1, ]
  gamma     <- fit$gamma
  orthogonal <- fit$params$orthogonal

  # Fisher info and covariance
  fit.Fisher <- FisherMatrix(X = X, beta_new = fit$beta,
                             orthogonal = orthogonal, gamma = gamma)
  cov.mat    <- MASS::ginv(fit.Fisher$Fn)   # = Sigma / n

  # Index helpers
  idx_mu <- seq_len(q)
  idx_j  <- q + (j - 1) * q + seq_len(q)

  # Joint 2q x 2q covariance block
  cov_joint <- cov.mat[c(idx_mu, idx_j), c(idx_mu, idx_j)]

  # Jacobians + projected covariances
  jac <- compute_proj_jacobians(mu, beta_j, cov_joint)

  P_mu       <- jac$P_mu
  P_perp     <- jac$P_perp
  s0         <- jac$s0
  Omega_loc_n  <- jac$Omega_loc_n    # Omega_loc  / n
  Omega_conc_n <- jac$Omega_conc_n   # Omega_conc / n

  # Projected estimates
  beta_loc  <- drop(P_perp %*% beta_j)
  beta_conc <- drop(P_mu   %*% beta_j)

  # Test statistics  W = proj' * ginv(Omega_n) * proj
  #   = n * proj' * ginv(Omega) * proj  ~ chi^2
  W_loc  <- as.numeric(t(beta_loc)  %*% MASS::ginv(Omega_loc_n)  %*% beta_loc)
  W_conc <- as.numeric(t(beta_conc) %*% MASS::ginv(Omega_conc_n) %*% beta_conc)

  p_loc  <- pchisq(W_loc,  df = q - 1, lower.tail = FALSE)
  p_conc <- pchisq(W_conc, df = 1,     lower.tail = FALSE)

  # 95% CI for scalar concentration effect  s = mu' beta_j / ||mu||
  mu_unit <- mu / sqrt(sum(mu^2))
  s_hat   <- sum(mu_unit * beta_j)
  var_s   <- as.numeric(t(mu_unit) %*% Omega_conc_n %*% mu_unit)
  se_s    <- sqrt(max(var_s, 0))
  z_a2    <- qnorm(1 - alpha / 2)
  ci_conc <- c(lower = s_hat - z_a2 * se_s,
               upper = s_hat + z_a2 * se_s)

  list(
    j          = j,
    mu_hat     = mu,
    beta_j_hat = beta_j,
    s0         = s0,
    # Projected estimates
    beta_loc   = beta_loc,
    beta_conc  = beta_conc,
    # Asymptotic covariances (= Omega / n)
    Omega_loc_n  = Omega_loc_n,
    Omega_conc_n = Omega_conc_n,
    # Test statistics
    W_loc  = W_loc,  df_loc  = q - 1, p_loc  = p_loc,
    W_conc = W_conc, df_conc = 1,     p_conc = p_conc,
    # Confidence intervals / regions
    ci_conc        = ci_conc,
    C_loc_centre   = beta_loc,
    C_loc_shape    = Omega_loc_n,
    C_loc_cv       = qchisq(1 - alpha, df = q - 1),
    alpha          = alpha
  )
}


# ------------------------------------------------------------------
#  print helper for inference_beta_j output
# ------------------------------------------------------------------
#' @export print_inference_beta_j
print_inference_beta_j <- function(res) {
  sep <- strrep("-", 58)
  cat(sep, "\n")
  cat(sprintf("  Projection-based inference  (covariate j = %d)\n", res$j))
  cat(sep, "\n")
  cat("  mu_hat     :", round(res$mu_hat,     4), "\n")
  cat("  beta_j_hat :", round(res$beta_j_hat, 4), "\n")
  cat(sprintf("  s0 = mu' beta_j / ||mu||^2 = %.4f\n\n", res$s0))
  cat("  Projected components\n")
  cat("    [Location]      P_{mu_perp} beta_j :", round(res$beta_loc,  4), "\n")
  cat("    [Concentration] P_{mu}      beta_j :", round(res$beta_conc, 4), "\n\n")
  cat("  Hypothesis tests\n")
  cat(sprintf("    H0_loc  (P_{mu_perp} beta_j = 0) : W = %8.3f, df = %d, p = %.4e  %s\n",
              res$W_loc,  res$df_loc,  res$p_loc,
              ifelse(res$p_loc  < res$alpha, "<-- REJECT", "")))
  cat(sprintf("    H0_conc (P_{mu}      beta_j = 0) : W = %8.3f, df = %d, p = %.4e  %s\n",
              res$W_conc, res$df_conc, res$p_conc,
              ifelse(res$p_conc < res$alpha, "<-- REJECT", "")))
  cat(sprintf("  %.0f%% CI for concentration (s = mu' beta_j / ||mu||): [%.4f, %.4f]\n",
              (1 - res$alpha) * 100, res$ci_conc["lower"], res$ci_conc["upper"]))
  cat(sep, "\n")
}


# ------------------------------------------------------------------
#  lrt_orthogonality
#
#  Theorem: 직교 제약 H0: beta_j perp mu (모든 j)에 대한 LRT
#    Lambda_n = 2 { ell(gamma_hat) - ell(gamma_hat_c) } ~ chi^2_p
#
#  Inputs:
#    fit_uncon : 비제약 vMFglm 객체
#    fit_con   : 제약(orthogonal=TRUE) vMFglm 객체
#    alpha     : 유의 수준
# ------------------------------------------------------------------
#' @export lrt_orthogonality
lrt_orthogonality <- function(fit_uncon, fit_con, alpha = 0.05) {

  ll_uncon <- tail(fit_uncon$loglik.list, 1)
  ll_con   <- tail(fit_con$loglik.list,   1)
  Lambda   <- 2 * (ll_uncon - ll_con)
  df       <- ncol(fit_uncon$X)                # p (number of constraints)
  p_val    <- pchisq(Lambda, df = df, lower.tail = FALSE)
  cv       <- qchisq(1 - alpha, df = df)

  cat(sprintf(
    "LRT (H0: beta_j perp mu for all j):\n  Lambda_n = %.4f, df = %d, p-value = %.4e\n",
    Lambda, df, p_val))
  cat(sprintf("  chi^2_{%d, %.2f} critical value = %.4f\n", df, alpha, cv))
  cat(if (p_val < alpha)
    "  => REJECT: use unconstrained model\n"
  else
    "  => FAIL TO REJECT: constrained model preferred (efficiency gain)\n")

  invisible(list(Lambda = Lambda, df = df, p_value = p_val))
}


# ------------------------------------------------------------------
#  confidence.vMFglm
#  (의존: summary.vMFglm의 새 cov.mu / cov.ortho)
# ------------------------------------------------------------------
#' @export confidence.vMFglm
confidence.vMFglm <- function(fit, conf.level = 0.95, j = NULL) {

  sm  <- summary(fit, conf.level = conf.level)
  p   <- ncol(fit$X)
  q   <- ncol(fit$Y)
  mu  <- fit$mu

  if (is.null(j)) j <- seq_len(p)

  lapply(j, function(jj) {
    beta_j    <- fit$beta[jj + 1, ]
    Omega_c_n <- sm$cov.mu[[jj]]      # rank-1 concentration
    Omega_l_n <- sm$cov.ortho[[jj]]   # rank-(q-1) location

    mu_unit <- mu / sqrt(sum(mu^2))
    s_hat   <- sum(mu_unit * beta_j)
    var_s   <- as.numeric(t(mu_unit) %*% Omega_c_n %*% mu_unit)
    se_s    <- sqrt(max(var_s, 0))
    z_a2    <- qnorm(1 - (1 - conf.level) / 2)
    ci_conc <- c(lower = s_hat - z_a2 * se_s,
                 upper = s_hat + z_a2 * se_s)

    list(
      j              = jj,
      # Location: ellipsoid centred at P_{mu_perp} beta_j_hat
      C_loc_centre   = drop(sm$projmat.ortho %*% beta_j),
      C_loc_shape    = Omega_l_n,
      C_loc_cv       = qchisq(conf.level, df = q - 1),
      # Concentration: interval for s = mu' beta_j / ||mu||
      s_hat          = s_hat,
      se_s           = se_s,
      ci_conc        = ci_conc
    )
  })
}
