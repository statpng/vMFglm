#  summary_vMFglm.R
#
#  주요 변경점:
#  ─────────────────────────────────────────────────────────────
#  [핵심] reduce_jacobian 파라미터 추가
#
#  reduce_jacobian = TRUE  (기본값, 권장):
#    축소 야코비안 사용
#    U_perp : SVD로 구한 μ⊥ 공간의 ONB  (q × (q-1))
#    J_perp_red = U_perp' J_perp_full  →  (q-1) × 2q,  full rank
#    Omega_loc_n : (q-1)×(q-1), full rank → solve() 사용
#    W_loc = proj_red' solve(Omega_loc_n) proj_red  ~  chi²(q-1) 정확
#    → ||β_j|| / ||μ|| 에 무관하게 type I error 제어
#
#  reduce_jacobian = FALSE (비교용):
#    원래 야코비안 사용
#    J_perp : q × 2q  →  Omega_loc_n : q×q, rank = q-1 (이론적)
#    W_loc = proj' ginv(Omega_loc_n) proj  ~  chi²(q-1) 점근적
#    → ||β_j|| / ||μ|| 이 클 때 유한 표본에서 type I error 팽창
#       (μ 방향 오염: ~ ||β_j||² / ||μ||² 크기)
#  ─────────────────────────────────────────────────────────────


# Helper: delta method Jacobians
#
#  reduce = TRUE:
#    J_perp_red  : (q-1) × 2q  (full rank)
#    Omega_loc_n : (q-1)×(q-1) (full rank, solve() 사용)
#    Omega_conc_n_scalar : 스칼라
#
#  reduce = FALSE:
#    J_perp_full : q × 2q
#    Omega_loc_n : q×q (rank q-1, ginv() 사용)
#    Omega_conc_n : q×q (rank 1, ginv() 사용)

compute_proj_jacobians <- function(mu, beta_j, cov_joint, reduce = TRUE) {

  q       <- length(mu)
  mu2     <- sum(mu^2)
  mu_norm <- sqrt(mu2)
  mu_unit <- mu / mu_norm
  s0      <- sum(mu * beta_j) / mu2
  P_mu    <- tcrossprod(mu_unit)
  P_perp  <- diag(q) - P_mu

  # ── Full J_perp  (q × 2q) ────────────────────────────────
  J_perp_full <- cbind(
    -tcrossprod(mu, beta_j) / mu2 + s0 * P_mu - s0 * P_perp,
    P_perp
  )

  # ── Full J_par  (q × 2q) ─────────────────────────────────
  J_par_full <- cbind(
    tcrossprod(mu, beta_j) / mu2 - s0 * P_mu + s0 * P_perp,
    P_mu
  )

  {

    # ── μ⊥ 공간의 ONB: U_perp (q × (q-1)) ─────────────────
    sv     <- svd(P_perp, nu = q, nv = 0)
    U_perp <- sv$u[, seq_len(q - 1L), drop = FALSE]

    # ── Reduced J_perp  (q-1) × 2q ─────────────────────────
    J_perp_red <- t(U_perp) %*% J_perp_full

    # ── 스칼라 Jacobian for g_par = μ_unit' β_j ────────────
    # dg_par/dμ   = P_perp β_j / ||μ||
    # dg_par/dβ_j = μ_unit
    J_par_red <- c(
      drop(P_perp %*% beta_j) / mu_norm,
      mu_unit
    )

    # ── 공분산 ─────────────────────────────────────────────
    Omega_loc_n         <- J_perp_red %*% cov_joint %*% t(J_perp_red)  # (q-1)×(q-1)
    Omega_loc_full_n    <- J_perp_full %*% cov_joint %*% t(J_perp_full)  # q×q
    Omega_conc_n <- as.numeric(t(J_par_red) %*% cov_joint %*% J_par_red)
    # Omega_conc_n_mat    <- Omega_conc_n_scalar * tcrossprod(mu_unit)    # q×q rank-1
    Omega_conc_full_n    <- J_par_full %*% cov_joint %*% t(J_par_full)  # q×q
    
    
    list(
      reduced             = TRUE,
      J_perp              = J_perp_red,
      J_perp_full         = J_perp_full,
      J_par_red           = J_par_red,
      J_par_full         = J_par_full,
      U_perp              = U_perp,
      P_mu                = P_mu,
      P_perp              = P_perp,
      s0                  = s0,
      Omega_loc_n         = Omega_loc_n,          # (q-1)×(q-1), full rank
      Omega_loc_full_n    = Omega_loc_full_n,          # (q-1)×(q-1), full rank
      Omega_conc_n        = Omega_conc_n,
      Omega_conc_full_n   = Omega_conc_full_n     # q×q rank-1, 하위 호환
      
    )

  }
}




# loglik.vmf, score, FisherMatrix, matsum  (변경 없음)

#' @export loglik.vmf
loglik.vmf <- function(X, Y, beta, gamma = 0) {

  if (length(gamma) == 1) gamma <- rep(gamma, ncol(X))
  q  <- ncol(Y)
  X1 <- cbind(1, X)

  Cq0 <- function(NORM, q) {
    (q/2 - 1) * log(NORM) -
      (q/2 * log(2 * pi) + log(besselI(NORM, q/2 - 1, expon.scaled = TRUE)) + NORM)
  }

  B     <- matrix(beta, ncol = q)
  Theta <- X1 %*% B
  NORM  <- sqrt(Rfast::rowsums(Theta^2))

  sum(Cq0(NORM, q)) + sum(Theta * Y) +
    sum(gamma * apply(beta[-1, , drop = FALSE], 1,
                      function(bj) crossprod(bj, beta[1, ])))
}


#' @export score
score <- function(X, Y, beta_new, gamma, orthogonal) {

  if (orthogonal) gamma <- rep(0, ncol(X))
  Offset <- matrix(0, nrow(Y), ncol(Y))
  n <- nrow(X); p <- ncol(X); q <- ncol(Y)
  bmu   <- beta_new[1, ]
  bB    <- beta_new[-1, , drop = FALSE]
  gamma <- as.vector(gamma)

  b1.list <- lapply(seq_len(n), function(i)
    b1.vMF(Offset[i, ] + bmu + t(bB) %*% X[i, ]))

  score_ll <- matsum(seq_len(n), function(i) {
    kronecker(c(1, X[i, ]), diag(1, q, q)) %*% (Y[i, ] - b1.list[[i]])
  })

  s1 <- -matsum(seq_len(p), function(j) gamma[j] * bB[j, ])
  s2 <- -kronecker(gamma, bmu)
  score_ll + c(s1, s2)
}


#' @export matsum
matsum <- function(seq, FUN) Reduce("+", lapply(seq, FUN))


#' @export FisherMatrix
FisherMatrix <- function(X, beta_new, orthogonal = FALSE, gamma = 0, OffSet = NULL) {

  n <- nrow(X); p <- ncol(X); q <- ncol(beta_new)
  if (length(gamma) == 1) gamma <- rep(gamma, p)

  bbeta   <- as.vector(t(beta_new[-1, , drop = FALSE]))
  bmu     <- beta_new[1, ]
  Xt.list <- apply(X, 1, function(Xi) kronecker(Xi, diag(1, q)), simplify = FALSE)
  b2.list <- lapply(seq_len(n), function(i)
    b2.vMF(bmu + crossprod(Xt.list[[i]], bbeta)))

  if (!orthogonal) {
    H11 <- -matsum(seq_len(n), function(i) b2.list[[i]])
    H22 <- -matsum(seq_len(n), function(i)
      Xt.list[[i]] %*% b2.list[[i]] %*% t(Xt.list[[i]]))
    H21 <- -matsum(seq_len(n), function(i) Xt.list[[i]] %*% b2.list[[i]])

    Fn <- matrix(NA, (p+1)*q, (p+1)*q)
    Fn[1:q, 1:q] <- -H11
    Fn[(q+1):((p+1)*q), (q+1):((p+1)*q)] <- -H22
    Fn[(q+1):((p+1)*q), 1:q] <- -H21
    Fn[1:q, (q+1):((p+1)*q)] <- t(-H21)

  } else {
    H11 <- -Reduce(`+`, b2.list)
    H12 <- -Reduce(`+`, lapply(seq_len(n),
                               function(i) b2.list[[i]] %*% t(Xt.list[[i]]))) -
      do.call("cbind", lapply(gamma, function(gj) gj * diag(1,q,q)))
    H13 <- -t(matrix(bbeta, p, q, byrow = TRUE))
    H22 <- -Reduce(`+`, lapply(seq_len(n),
                               function(i) Xt.list[[i]] %*% b2.list[[i]] %*% t(Xt.list[[i]])))
    H23 <- -kronecker(diag(p), bmu)
    H33 <- matrix(0, p, p)

    Fn <- matrix(0, (p+1)*q+p, (p+1)*q+p)
    Fn[1:q, 1:q] <- -H11
    Fn[1:q, (q+1):((p+1)*q)] <- -H12
    Fn[1:q, ((p+1)*q+1):((p+1)*q+p)] <- -H13
    Fn[(q+1):((p+1)*q), (q+1):((p+1)*q)] <- -H22
    Fn[(q+1):((p+1)*q), ((p+1)*q+1):((p+1)*q+p)] <- -H23
    Fn[((p+1)*q+1):((p+1)*q+p), ((p+1)*q+1):((p+1)*q+p)] <- -H33
    Fn[lower.tri(Fn)] <- t(Fn)[lower.tri(Fn)]
  }

  idx.list <- lapply(seq_len(p+1), function(j) (j-1)*q + seq_len(q))
  list(b2.list = b2.list, Fn = Fn,
       Fnj = lapply(idx.list, function(idx) Fn[idx, idx]))
}


Cq <- function(theta, logarithm = TRUE) {
  q <- length(theta); NORM <- norm(theta, "2")
  val <- (q/2-1)*log(NORM) -
    (q/2*log(2*pi) + log(besselI(NORM, q/2-1, expon.scaled=TRUE)) + NORM)
  if (logarithm) val else exp(val)
}

Bq <- function(theta) {
  q <- length(theta); NORM <- norm(theta, "2")
  I0 <- besselI(NORM, nu=q/2-1, expon.scaled=TRUE)
  I1 <- besselI(NORM, nu=q/2,   expon.scaled=TRUE)
  if (I1 == 0) NORM / ((q/2-1) + sqrt(NORM^2 + (q/2-1)^2)) else I1/I0
}

Hq <- function(theta) {
  q <- length(theta); NORM <- norm(theta, "2")
  Bq.val <- Bq(theta)
  1 - Bq.val^2 - (q-1)/NORM * Bq.val
}

subgrad <- function(theta) {
  NORM <- norm(theta, "2")
  if (NORM == 0) { v <- runif(length(theta),-1,1); v/(norm(v,"2")*1.1) }
  else theta/NORM
}

#' @export b1.vMF
b1.vMF <- function(theta) Bq(theta) * subgrad(theta)

#' @export b2.vMF
b2.vMF <- function(theta) {
  q <- length(theta); NORM <- norm(theta, "2")
  if (NORM < 1e-14) return((1/q)*diag(1,q))
  s <- subgrad(theta)
  tcrossprod(s)*Hq(theta) + Bq(theta)/NORM*(diag(1,q)-tcrossprod(s))
}


diag.matrix <- function(X, lambda) {
  n <- dim(X)[1]; p <- dim(X)[2]
  if (n != p) stop("Non-square matrix.")
  diag(lambda, n, p)
}




# summary.vMFglm  (main function)

# #' @method summary vMFglm
# #' @export
# summary.vMFglm <- function(fit,
#                               type            = "wald",
#                               conf.level      = 0.95,
#                               mu_true         = NULL,
#                               beta_true       = NULL,
#                               reduce_jacobian = TRUE) {
#   
#   # reduce_jacobian = TRUE  (기본값, 권장):
#   #   축소 야코비안 사용 → 유한 표본에서도 type I error 정확 제어
#   #
#   # reduce_jacobian = FALSE (비교용):
#   #   원래 야코비안 사용 → ||β||/||μ|| 클 때 유한 표본 type I error 팽창
# 
#   if(!is.null(beta_true)){
#     beta_true <- as.matrix(beta_true)
#   }
#   
# 
#   if (inherits(fit, "try-error")) {
#     na_result <- list(beta=NA, norm=NA, Fisher=NA, orthogonal=NA,
#                       wald=NA, df=NA,
#                       statistic=list(wald=NA, LRT=NA, score=NA),
#                       pvalue   =list(wald=NA, LRT=NA, score=NA))
#     return(na_result)
#   }
# 
#   if (length(type) > 1) type <- "wald"
#   type.wald  <- type %in% c("all", "wald")
#   type.LRT   <- type %in% c("all", "LRT")
#   type.score <- type %in% c("all", "score")
# 
# 
#   # 0. 기본 추출 ----
#   X        <- fit$X;  Y <- fit$Y;  Offset <- fit$offset
#   p        <- ncol(X); q <- ncol(Y); n <- nrow(Y)
#   gamma    <- fit$gamma
#   beta_new <- fit$beta
#   mu       <- fit$mu
# 
#   orthogonal <- fit$params$orthogonal
#   df         <- ifelse(orthogonal, q - 1, q)
#   if (is.null(beta_true)) beta_true <- matrix(0, p, q)
# 
# 
#   # 1. 귀무 적합 (LRT / Score 전용) ----
#   # betahat0 <- if (type.LRT || type.score) {
#   #   lapply(seq_len(p), function(j) {
#   #     if (orthogonal && p > 1) {
#   #       fit0  <- vMFglm(X=X[,-j,drop=FALSE], Y=Y,
#   #                          orthogonal=orthogonal, gamma=gamma)
#   #       beta0 <- matrix(0, p+1, q)
#   #       beta0[1,]          <- fit0$beta[1,]
#   #       beta0[-c(1,j+1),] <- fit0$beta[-1,]
#   #       beta0
#   #     } else {
#   #       pf <- rep(0, p); pf[j] <- 1
#   #       vMFglm(X=X, Y=Y, orthogonal=orthogonal,
#   #                 gamma=gamma, penalty.factor=pf)$beta
#   #     }
#   #   })
#   # } else NULL
# 
# 
#   # 2. Fisher 정보 + 공분산 행렬 ----
#   fit.Fisher <- FisherMatrix(X=X, beta_new=beta_new,
#                              orthogonal=orthogonal, gamma=gamma)
#   Fn        <- fit.Fisher$Fn
#   F_nj.list <- fit.Fisher$Fnj
# 
#   dim_gamma <- (p + 1) * q
#   cov.mat   <- MASS::ginv(Fn)[seq_len(dim_gamma), seq_len(dim_gamma)]
#   idx.list  <- lapply(seq_len(p+1), function(j) (j-1)*q + seq_len(q))
# 
# 
#   # 3. 조건부 공분산 (Wald용) ----
#   idx_mu <- idx.list[[1]]
#   S11 <- cov.mat[idx_mu,  idx_mu]
#   S12 <- cov.mat[idx_mu,  -idx_mu]
#   S21 <- cov.mat[-idx_mu, idx_mu]
#   S22 <- cov.mat[-idx_mu, -idx_mu]
# 
# 
#   # 4. 사영 행렬 ----
#   mu_unit      <- mu / sqrt(sum(mu^2))
#   projmat.mu   <- tcrossprod(mu_unit)
#   projmat.perp <- diag(q) - projmat.mu
# 
# 
#   # 5. Jacobian 기반 공분산 ----
#   proj_jac_list     <- vector("list", p)
#   Omega_loc_n_list  <- vector("list", p)
#   Omega_loc_full_n_list  <- vector("list", p)
#   Omega_conc_n_list <- vector("list", p)
#   Omega_conc_full_n_list  <- vector("list", p)
#   U_perp_list <- vector("list", p)
#   
#   for (j in seq_len(p)) {
#     beta_j    <- beta_new[j + 1, ]
#     idx_j     <- idx.list[[j + 1]]
#     cov_joint <- cov.mat[c(idx_mu, idx_j), c(idx_mu, idx_j)]
# 
#     jac <- compute_proj_jacobians(mu, beta_j, cov_joint,
#                                   reduce = reduce_jacobian)
# 
#     proj_jac_list[[j]]     <- jac
#     Omega_loc_n_list[[j]]  <- jac$Omega_loc_n
#     Omega_loc_full_n_list[[j]]  <- jac$Omega_loc_full_n
#     Omega_conc_n_list[[j]] <- jac$Omega_conc_n
#     Omega_conc_full_n_list[[j]] <- jac$Omega_conc_full_n
#   }
# 
# 
#   # 6. Wald 검정 ----
#   rnames <- c("global", paste0("ind", seq_len(p)))
# 
#   make_test_mat <- function(stat_global, stat_ind, pv_global, pv_ind) {
#     mat <- cbind(stat   = c(stat_global, stat_ind),
#                  pvalue = c(pv_global,   pv_ind))
#     rownames(mat) <- rnames
#     mat
#   }
# 
#   TestResult <- NULL
# 
#   if (type.wald) {
#     Wald.global <- {
#       bhat     <- as.vector(t(beta_new[-1, ]))
#       idx_beta <- seq(q+1, (p+1)*q)
#       S_global <- cov.mat[idx_beta, idx_beta]
#       as.numeric(t(bhat) %*% solve(S_global) %*% bhat)
#     }
#     pv.wald.global <- 1 - pchisq(Wald.global, df = p * df)
# 
#     Wald.j.list <- sapply(seq_len(p), function(j) {
#       idx_j <- idx.list[[j+1]]
#       S_j   <- cov.mat[idx_j, idx_j]
#       as.numeric(t(beta_new[j+1,]) %*% solve(S_j) %*% beta_new[j+1,])
#     })
#     pv.wald.ind <- sapply(Wald.j.list, function(x) 1 - pchisq(x, df))
# 
#     TestResult$wald <- make_test_mat(Wald.global, Wald.j.list,
#                                      pv.wald.global, pv.wald.ind)
#   }
# 
# 
#   # 7. LRT ----
#   # if (type.LRT) {
#   #   beta0      <- rbind(fit$mu, matrix(0, p, q))
#   #   LRT.global <- -2*(loglik.vmf(X,Y,beta0,gamma=gamma) -
#   #                       loglik.vmf(X,Y,beta_new,gamma=gamma))
#   #   LRT.ind <- if (p==1) LRT.global else
#   #     sapply(betahat0, function(b0)
#   #       -2*(loglik.vmf(X,Y,b0,gamma=gamma) -
#   #             loglik.vmf(X,Y,beta_new,gamma=gamma)))
#   # 
#   #   TestResult$LRT <- make_test_mat(
#   #     LRT.global, LRT.ind,
#   #     1 - pchisq(LRT.global, p*df),
#   #     sapply(LRT.ind, function(x) 1 - pchisq(x, df)))
#   # }
# 
# 
#   # 8. Score 검정 ----
#   # if (type.score) {
#   #   beta0.score <- rbind(fit$mu, matrix(1e-12, p, q))
#   #   fsh0  <- FisherMatrix(X=X, beta_new=beta0.score,
#   #                         orthogonal=orthogonal, gamma=gamma)
#   #   sv0   <- score(X, Y, beta0.score, gamma=gamma, orthogonal=orthogonal)
#   #   COV0  <- suppressWarnings(
#   #     solve(fsh0$Fn + diag(1e-12, nrow(fsh0$Fn), nrow(fsh0$Fn))))
#   # 
#   #   Score.global <- as.numeric(
#   #     t(sv0[-(1:q)]) %*% COV0[(q+1):(p*q+q),(q+1):(p*q+q)] %*% sv0[-(1:q)])
#   # 
#   #   Score.ind <- if (p==1) Score.global else
#   #     sapply(seq_len(p), function(j) {
#   #       b0j   <- betahat0[[j]]
#   #       fshj  <- FisherMatrix(X=X, beta_new=b0j,
#   #                             orthogonal=orthogonal, gamma=gamma)
#   #       svj   <- score(X,Y,b0j, gamma=gamma, orthogonal=orthogonal)
#   #       idx   <- seq_len(q) + q*j
#   #       cov0j <- MASS::ginv(fshj$Fn)[idx, idx]
#   #       as.numeric(t(svj[idx]) %*% cov0j %*% svj[idx])
#   #     })
#   # 
#   #   TestResult$score <- make_test_mat(
#   #     Score.global, Score.ind,
#   #     1 - pchisq(Score.global, p*df),
#   #     sapply(Score.ind, function(x) 1 - pchisq(x, df)))
#   # }
# 
# 
#   # 9. Projection 검정 ----
#   #
#   # ── reduce_jacobian = TRUE  (축소) ────────────────────────
#   #  W_loc  = proj_red' solve(Omega_loc_n) proj_red  ~  chi²(q-1)
#   #    proj_red = U_perp' β̂_j      [(q-1) 벡터]
#   #    Omega_loc_n : (q-1)×(q-1),  solve() 사용
#   #    → ||β||/||μ|| 에 무관하게 유한 표본 type I error 제어
#   #      이유: U_perp' μ_unit = 0 (항등식) → μ 방향 오염 대수적 소거
#   #
#   #  W_conc = s_hat² / Omega_conc_scalar  ~  chi²(1)
#   #
#   # ── reduce_jacobian = FALSE (비축소) ──────────────────────
#   #  W_loc  = proj' ginv(Omega_loc_n) proj  ~  chi²(q-1) 점근적
#   #    proj    = P_perp β̂_j         [q 벡터]
#   #    Omega_loc_n : q×q rank-(q-1), ginv() 사용
#   #    → μ 방향 오염 크기 ~ ||β_j||²/||μ||²
#   #      s_loc=3, ||μ||=10: 오염 ~ 9/100 → type I error ≈ 0.10 (기대: 0.05)
#   #
#   #  W_conc = proj_mu' ginv(Omega_conc_n) proj_mu  ~  chi²(1) 점근적
#   {
#     if (reduce_jacobian) {
# 
#       # ── W_loc (축소) ───────────────────────────────────────
#       Wloc.ind <- sapply(seq_len(p), function(j) {
#         beta_j   <- beta_new[j+1, ]
#         U_perp   <- proj_jac_list[[j]]$U_perp       # q × (q-1)
#         proj_red <- drop(t(U_perp) %*% beta_j)      # (q-1) 벡터
#         Omega_loc_n  <- Omega_loc_n_list[[j]]            # (q-1)×(q-1)
#         as.numeric(t(proj_red) %*% solve(Omega_loc_n) %*% proj_red)
#       })
# 
#       # ── W_conc (스칼라) ────────────────────────────────────
#       Wconc.ind <- sapply(seq_len(p), function(j) {
#         beta_j   <- beta_new[j+1, ]
#         s_hat    <- sum(mu_unit * beta_j)
#         Omega_sc <- proj_jac_list[[j]]$Omega_conc_n
#         if (Omega_sc < 1e-14) return(0)
#         s_hat^2 / Omega_sc
#       })
# 
#       # ── 집중도 신뢰 구간 ────────────────────────────────────
#       CI.Wconc <- unlist(lapply(seq_len(p), function(j) {
#         beta_j   <- beta_new[j+1, ]
#         s_hat    <- sum(mu_unit * beta_j)
#         Omega_sc <- proj_jac_list[[j]]$Omega_conc_n
#         se_s     <- sqrt(max(Omega_sc, 0))
#         z        <- qnorm(1 - (1 - conf.level) / 2)
#         ci       <- c(s_hat - z*se_s, s_hat + z*se_s)
#         names(ci) <- c(paste0("lower.",j), paste0("upper.",j))
#         ci
#       }))
# 
#     } else {
# 
#       # ── W_loc (비축소) ─────────────────────────────────────
#       Wloc.ind <- sapply(seq_len(p), function(j) {
#         beta_j  <- beta_new[j+1, ]
#         proj    <- drop(projmat.perp %*% beta_j)    # q 벡터
#         Omega_loc_full_n <- Omega_loc_full_n_list[[j]]             # q×q, rank q
#         as.numeric(t(proj) %*% solve(Omega_loc_full_n) %*% proj)
#       })
# 
#       # ── W_conc (비축소) ────────────────────────────────────
#       Wconc.ind <- sapply(seq_len(p), function(j) {
#         beta_j  <- beta_new[j+1, ]
#         proj_mu <- drop(projmat.mu %*% beta_j)      # q 벡터
#         Omega_conc_n <- Omega_conc_full_n_list[[j]]            # q×q, rank 1
#         as.numeric(t(proj_mu) %*% solve(Omega_conc_n) %*% proj_mu)
#       })
# 
#       # ── 집중도 신뢰 구간 ────────────────────────────────────
#       CI.Wconc <- unlist(lapply(seq_len(p), function(j) {
#         beta_j  <- beta_new[j+1, ]
#         s_hat   <- sum(mu_unit * beta_j)
#         Omega_conc_full_n <- Omega_conc_full_n_list[[j]]
#         var_s   <- as.numeric(t(mu_unit) %*% Omega_conc_full_n %*% mu_unit)
#         se_s    <- sqrt(max(var_s, 0))
#         z       <- qnorm(1 - (1 - conf.level) / 2)
#         ci      <- c(s_hat - z*se_s, s_hat + z*se_s)
#         names(ci) <- c(paste0("lower.",j), paste0("upper.",j))
#         ci
#       }))
#     }
# 
#     Wloc.global    <- sum(Wloc.ind)
#     pv.Wloc.global <- 1 - pchisq(Wloc.global,  p*(q-1))
#     pv.Wloc.ind    <- sapply(Wloc.ind,  function(x) 1 - pchisq(x, q-1))
# 
#     Wconc.global    <- sum(Wconc.ind)
#     pv.Wconc.global <- 1 - pchisq(Wconc.global, p*1)
#     pv.Wconc.ind    <- sapply(Wconc.ind, function(x) 1 - pchisq(x, 1))
# 
#     TestResult$Wloc <- TestResult$Wloc <- make_test_mat(
#       Wloc.global,  Wloc.ind,  pv.Wloc.global,  pv.Wloc.ind)
#     TestResult$Wconc <- make_test_mat(Wconc.global, Wconc.ind, pv.Wconc.global, pv.Wconc.ind)
#     TestResult$CI.Wconc <- CI.Wconc
#     TestResult$s0 <- sapply(seq_len(p), function(j) proj_jac_list[[j]]$s0)
#   }
# 
# 
#   # 10. 신뢰 영역 포함 확률 판정 ----
#   #
#   # reduce_jacobian = TRUE:
#   #   T_loc = (U_hat'(β̂_j - P_perp0 β_true))' solve(Ω_loc) (U_hat'(β̂_j - P_perp0 β_true))
#   #   U_hat : Omega_loc_n 과 동일한 mu_hat 좌표계 → chi²(q-1) 보장
#   #   P_perp0 : mu_true 기준으로 β_true 의 ⊥ 성분 추출 (μ 방향 제거)
#   #
#   # reduce_jacobian = FALSE:
#   #   T_loc = (P_perp0(β̂_j - β_true))' ginv(P_perp0 Ω_loc P_perp0) (P_perp0(β̂_j - β_true))
#   conf.inside <- list(wald=NULL, Wconc=NULL, Wloc=NULL)
#   
#   z_a2 <- qnorm(1 - (1 - conf.level) / 2)
# 
#   for (j in seq_len(p)) {
#     beta_j <- beta_new[j+1, ]
#     idx_j  <- idx.list[[j+1]]
# 
#     if (!is.null(mu_true)) {
#       mu_true    <- as.vector(mu_true)
#       mu_unit_0  <- mu_true / sqrt(sum(mu_true^2))
#       P_mu_0     <- tcrossprod(mu_unit_0)
#       P_perp_0   <- diag(q) - P_mu_0
#       s_true_scl <- sum(mu_unit_0 * beta_true[j,])
#     } else {
#       P_mu_0     <- projmat.mu
#       P_perp_0   <- projmat.perp
#       s_true_scl <- sum(mu_unit * beta_true[j,])
#     }
#     
#     # ── Wald CI (공통) ────────────────────────────────────────
#     S_j     <- cov.mat[idx_j, idx_j]
#     S_j_inv <- solve(S_j)
#     conf.inside$wald[j] <- ifelse(
#       as.numeric(t(beta_j - beta_true[j,]) %*% S_j_inv %*% (beta_j - beta_true[j,])) <=
#         qchisq(conf.level, q),
#       "inside", "outside")
# 
#     # ── W_loc CI ─────────────────────────────────────────────
#     cv_loc  <- qchisq(conf.level, q - 1)
#     Omega_loc_n <- Omega_loc_n_list[[j]]
#     Omega_loc_full_n <- Omega_loc_full_n_list[[j]]
# 
#     if (reduce_jacobian) {
#       # U_perp_hat 으로 통일 → Omega_n 과 같은 좌표계 → chi²(q-1) 보장
#       U_perp_hat <- proj_jac_list[[j]]$U_perp
#       perp_true  <- drop(P_perp_0 %*% beta_true[j,])     # mu_true 기준 ⊥ 성분
#       diff_red   <- drop(t(U_perp_hat) %*% (beta_j - perp_true))
#       T_loc      <- as.numeric(t(diff_red) %*% solve(Omega_loc_n) %*% diff_red)
#     } else {
#       
#       P_perp_hat <- projmat.perp
#       P_perp_true  <- P_perp_0    # mu_true 기준 ⊥ 성분
#       diff   <- drop( P_perp_hat %*% beta_j - P_perp_true %*% beta_true[j,] )
#       T_loc      <- as.numeric(t(diff) %*% solve(Omega_loc_full_n) %*% diff)
#       
#     }
# 
#     conf.inside$Wloc[j] <- ifelse(T_loc <= cv_loc, "inside", "outside")
# 
#     # ── W_conc CI (공통: 스칼라 공식) ────────────────────────
#     s_hat <- sum(mu_unit * beta_j)
# 
#     if (reduce_jacobian) {
#       Omega_sc <- proj_jac_list[[j]]$Omega_conc_n
#       se_s     <- sqrt(max(Omega_sc, 0))
#     } else {
#       Omega_cn <- Omega_conc_full_n_list[[j]]
#       var_s    <- as.numeric(t(mu_unit) %*% Omega_cn %*% mu_unit)
#       se_s     <- sqrt(max(var_s, 0))
#     }
# 
#     conf.inside$Wconc[j] <- ifelse(
#       abs(s_hat - s_true_scl) <= z_a2 * se_s,
#       "inside", "outside")
#   }
# 
# 
#   # 11. 반환 ----
#   list(
#     beta            = beta_new,
#     norm            = apply(beta_new, 1, norm, "2"),
#     Fisher          = list(Fn=Fn, Fnj=F_nj.list),
#     orthogonal      = fit$params$orthogonal,
#     wald            = if (type.wald) Wald.j.list else NULL,
#     df              = df,
#     idx.list        = idx.list,
#     cov.mat         = cov.mat,
#     cov.mu          = Omega_conc_n_list,
#     cov.mu.full     = Omega_conc_full_n_list,
#     cov.ortho       = Omega_loc_n_list,
#     cov.ortho.full  = Omega_loc_full_n_list,
#     projmat.mu      = projmat.mu,
#     projmat.ortho   = projmat.perp,
#     s0              = TestResult$s0,
#     test            = TestResult,
#     conf.inside     = as.data.frame(conf.inside),
#     mu_true_used    = !is.null(mu_true),
#     U_perp          = lapply(seq_len(p), function(j) proj_jac_list[[j]]$U_perp),
#     reduce_jacobian = reduce_jacobian
#   )
# }














#  summary_fast.R : summary.vMFglm 개선 버전
#
#  원본 대비 변경점
#  (1) [버그] 인자에서 사라진 type이 본문에서 참조되어 즉시 에러
#      -> type 인자 복원 (기본 "wald")
#  (2) [버그] TestResult$Wloc 이중 대입 오타 제거
#  (3) FisherMatrix -> FisherMatrix.fast (루프/kronecker 리스트 제거)
#      + fit$offset을 Fisher 계산에 반영 (원본은 무시)
#  (4) MASS::ginv 전체 역행렬 -> 정칙이면 촐레스키(chol2inv),
#      실패 시에만 고유값 기반 pseudo-inverse로 폴백
#  (5) U_perp: SVD -> QR 여공간 (compute_proj_jacobians.fast)
#  (6) 1 - pchisq(x, df) -> pchisq(x, df, lower.tail = FALSE)
#      (극단 p값에서 수치적으로 안정)
#  (7) Omega_loc의 촐레스키를 검정통계량과 신뢰영역 판정에서 재사용
#  (8) 반환 구조는 원본과 동일 (하위 호환)
#
#  사전 조건: source("vmf_fast.R")  (FisherMatrix.fast,
#             compute_proj_jacobians.fast, vmf_moments 등)

# ---- 안전한 대칭 역행렬: chol 우선, 실패 시 eigen pseudo-inverse ----
sym_inv <- function(A, tol = 1e-12) {
  A <- (A + t(A)) / 2
  out <- tryCatch(chol2inv(chol(A)), error = function(e) NULL)
  if (!is.null(out)) return(out)
  ei <- eigen(A, symmetric = TRUE)
  pos <- ei$values > tol * max(abs(ei$values))
  ei$vectors[, pos, drop = FALSE] %*%
    ((1 / ei$values[pos]) * t(ei$vectors[, pos, drop = FALSE]))
}

# x' A^{-1} x 를 chol 캐시로 계산
quad_inv <- function(cholA, x) {
  z <- forwardsolve(t(cholA), x)
  sum(z^2)
}



compute_proj_jacobians.fast <- function(mu, beta_j, cov_joint, reduce = TRUE) {
  q       <- length(mu)
  mu2     <- sum(mu^2)
  mu_norm <- sqrt(mu2)
  mu_unit <- mu / mu_norm
  s0      <- sum(mu * beta_j) / mu2
  P_mu    <- tcrossprod(mu_unit)
  P_perp  <- diag(q) - P_mu
  
  J_perp_full <- cbind(-tcrossprod(mu, beta_j)/mu2 + s0*P_mu - s0*P_perp, P_perp)
  J_par_full  <- cbind( tcrossprod(mu, beta_j)/mu2 - s0*P_mu + s0*P_perp, P_mu)
  
  U_perp     <- orth_complement(mu)                        # SVD 대체
  J_perp_red <- crossprod(U_perp, J_perp_full)             # (q-1) x 2q
  J_par_red  <- c(drop(P_perp %*% beta_j)/mu_norm, mu_unit)
  
  CJt_perp <- tcrossprod(cov_joint, J_perp_red)            # 중간곱 재사용
  CJt_full <- tcrossprod(cov_joint, J_perp_full)
  
  list(reduced = reduce,
       J_perp = J_perp_red, J_perp_full = J_perp_full,
       J_par_red = J_par_red, J_par_full = J_par_full,
       U_perp = U_perp, P_mu = P_mu, P_perp = P_perp, s0 = s0,
       Omega_loc_n       = J_perp_red  %*% CJt_perp,
       Omega_loc_full_n  = J_perp_full %*% CJt_full,
       Omega_conc_n      = drop(crossprod(J_par_red, cov_joint %*% J_par_red)),
       Omega_conc_full_n = J_par_full %*% tcrossprod(cov_joint, J_par_full))
}


# ---- Fisher 행렬 (루프/kronecker 리스트 없음) ----
FisherMatrix.fast <- function(X, beta_new, orthogonal = FALSE, gamma = 0,
                              Offset = NULL) {
  n <- nrow(X); p <- ncol(X); q <- ncol(beta_new)
  if (length(gamma) == 1) gamma <- rep(gamma, p)
  X1 <- cbind(1, X)
  Theta <- X1 %*% beta_new
  if (!is.null(Offset)) Theta <- Theta + Offset
  mo <- vmf_moments(Theta, q)
  
  # F0 = sum_i (x1_i x1_i') ⊗ W_i,  W_i = c1_i I + c2_i s_i s_i'
  G1 <- crossprod(X1, mo$c1 * X1)                    # (p+1) x (p+1)
  M  <- rowKron(X1, mo$S)                            # n x (p+1)q
  F0 <- kronecker(G1, diag(q)) + crossprod(M, mo$c2 * M)
  
  idx.list <- lapply(seq_len(p + 1), function(j) (j - 1) * q + seq_len(q))
  
  if (!orthogonal) {
    Fn <- F0
  } else {
    # 라그랑주 확장: (p+1)q + p 차원, 원 코드와 동일 블록
    d  <- (p + 1) * q
    Fn <- matrix(0, d + p, d + p)
    Fn[1:d, 1:d] <- F0
    # H12 페널티 블록 (mu행 x beta_j열)에 gamma_j I 추가
    for (j in seq_len(p)) {
      Fn[1:q, idx.list[[j + 1]]] <- Fn[1:q, idx.list[[j + 1]]] + gamma[j] * diag(q)
    }
    bB <- beta_new[-1, , drop = FALSE]
    Fn[1:q, d + seq_len(p)] <- t(bB)                 # -H13 = t(B)
    for (j in seq_len(p)) Fn[idx.list[[j + 1]], d + j] <- beta_new[1, ]  # -H23
    Fn[lower.tri(Fn)] <- t(Fn)[lower.tri(Fn)]
  }
  
  list(Fn = Fn,
       Fnj = lapply(idx.list, function(idx) Fn[idx, idx]),
       moments = mo)
}




#' @method summary vMFglm
#' @export
summary.vMFglm <- function(object,
                           type            = "wald",
                           conf.level      = 0.95,
                           mu_true         = NULL,
                           beta_true       = NULL,
                           reduce_jacobian = FALSE,
                           ...) {
  fit <- object
  
  if (inherits(fit, "try-error")) {
    return(list(beta = NA, norm = NA, Fisher = NA, orthogonal = NA,
                wald = NA, df = NA,
                statistic = list(wald = NA, LRT = NA, score = NA),
                pvalue    = list(wald = NA, LRT = NA, score = NA)))
  }
  
  if (length(type) > 1) type <- "wald"
  type.wald <- type %in% c("all", "wald")
  
  # 0. 기본 추출 ----
  X <- fit$X; Y <- fit$Y; Offset <- fit$offset
  p <- ncol(X); q <- ncol(Y); n <- nrow(Y)
  gamma      <- fit$gamma
  beta_new   <- fit$beta
  mu         <- fit$mu
  orthogonal <- isTRUE(fit$params$orthogonal)
  df         <- if (orthogonal) q - 1 else q
  
  if (!is.null(beta_true)) beta_true <- as.matrix(beta_true)
  if (is.null(beta_true))  beta_true <- matrix(0, p, q)
  
  # 1. Fisher 정보 + 공분산 ----
  fit.Fisher <- FisherMatrix.fast(X = X, beta_new = beta_new,
                                  orthogonal = orthogonal, gamma = gamma,
                                  Offset = Offset)
  Fn <- fit.Fisher$Fn
  
  dim_gamma <- (p + 1) * q
  cov.full  <- sym_inv(Fn)                       # 정칙이면 chol, 아니면 pinv
  cov.mat   <- cov.full[seq_len(dim_gamma), seq_len(dim_gamma)]
  idx.list  <- lapply(seq_len(p + 1), function(j) (j - 1) * q + seq_len(q))
  idx_mu    <- idx.list[[1]]
  
  # 2. 사영 행렬 ----
  mu_unit      <- mu / sqrt(sum(mu^2))
  projmat.mu   <- tcrossprod(mu_unit)
  projmat.perp <- diag(q) - projmat.mu
  
  # 3. Jacobian 기반 공분산 (관측별 캐시) ----
  proj_jac_list <- vector("list", p)
  chol_loc_list <- vector("list", p)   # Omega_loc(축소/비축소)의 chol 캐시
  for (j in seq_len(p)) {
    beta_j    <- beta_new[j + 1, ]
    idx_j     <- idx.list[[j + 1]]
    cov_joint <- cov.mat[c(idx_mu, idx_j), c(idx_mu, idx_j)]
    jac <- compute_proj_jacobians.fast(mu, beta_j, cov_joint,
                                       reduce = reduce_jacobian)
    proj_jac_list[[j]] <- jac
    Om <- if (reduce_jacobian) jac$Omega_loc_n else jac$Omega_loc_full_n
    chol_loc_list[[j]] <- tryCatch(chol((Om + t(Om)) / 2),
                                   error = function(e) NULL)
  }
  Omega_loc_n_list       <- lapply(proj_jac_list, `[[`, "Omega_loc_n")
  Omega_loc_full_n_list  <- lapply(proj_jac_list, `[[`, "Omega_loc_full_n")
  Omega_conc_n_list      <- lapply(proj_jac_list, `[[`, "Omega_conc_n")
  Omega_conc_full_n_list <- lapply(proj_jac_list, `[[`, "Omega_conc_full_n")
  
  # W_loc용 이차형식 (chol 캐시 사용, 실패 시 pinv)
  qform_loc <- function(j, v) {
    if (!is.null(chol_loc_list[[j]])) quad_inv(chol_loc_list[[j]], v)
    else {
      Om <- if (reduce_jacobian) Omega_loc_n_list[[j]] else Omega_loc_full_n_list[[j]]
      drop(crossprod(v, sym_inv(Om) %*% v))
    }
  }
  
  rnames <- c("global", paste0("ind", seq_len(p)))
  make_test_mat <- function(sg, si, pg, pi) {
    mat <- cbind(stat = c(sg, si), pvalue = c(pg, pi))
    rownames(mat) <- rnames
    mat
  }
  TestResult <- NULL
  
  # 4. Wald 검정 ----
  Wald.j.list <- NULL
  if (type.wald) {
    idx_beta  <- seq(q + 1, (p + 1) * q)
    S_global  <- cov.mat[idx_beta, idx_beta]
    bhat      <- as.vector(t(beta_new[-1, , drop = FALSE]))
    Wald.global <- drop(crossprod(bhat, sym_inv(S_global) %*% bhat))
    
    Wald.j.list <- sapply(seq_len(p), function(j) {
      idx_j <- idx.list[[j + 1]]
      S_j   <- cov.mat[idx_j, idx_j]
      bj    <- beta_new[j + 1, ]
      drop(crossprod(bj, sym_inv(S_j) %*% bj))
    })
    
    TestResult$wald <- make_test_mat(
      Wald.global, Wald.j.list,
      pchisq(Wald.global, p * df, lower.tail = FALSE),
      pchisq(Wald.j.list, df,     lower.tail = FALSE))
  }
  
  # 5. Projection 검정 (W_loc, W_conc) ----
  z_a2 <- qnorm(1 - (1 - conf.level) / 2)
  
  Wloc.ind <- sapply(seq_len(p), function(j) {
    beta_j <- beta_new[j + 1, ]
    v <- if (reduce_jacobian)
      drop(crossprod(proj_jac_list[[j]]$U_perp, beta_j))
    else
      drop(projmat.perp %*% beta_j)
    qform_loc(j, v)
  })
  
  s_hat_vec <- drop(beta_new[-1, , drop = FALSE] %*% mu_unit)  # 길이 p
  
  Wconc.ind <- sapply(seq_len(p), function(j) {
    if (reduce_jacobian) {
      Om <- Omega_conc_n_list[[j]]
      if (Om < 1e-14) 0 else s_hat_vec[j]^2 / Om
    } else {
      proj_mu <- drop(projmat.mu %*% beta_new[j + 1, ])
      Om <- Omega_conc_full_n_list[[j]]
      drop(crossprod(proj_mu, sym_inv(Om) %*% proj_mu))
    }
  })
  
  se_conc <- sapply(seq_len(p), function(j) {
    if (reduce_jacobian) sqrt(max(Omega_conc_n_list[[j]], 0))
    else sqrt(max(drop(crossprod(mu_unit,
                                 Omega_conc_full_n_list[[j]] %*% mu_unit)), 0))
  })
  CI.Wconc <- unlist(lapply(seq_len(p), function(j) {
    ci <- s_hat_vec[j] + c(-1, 1) * z_a2 * se_conc[j]
    names(ci) <- paste0(c("lower.", "upper."), j)
    ci
  }))
  
  Wloc.global  <- sum(Wloc.ind)
  Wconc.global <- sum(Wconc.ind)
  
  TestResult$Wloc <- make_test_mat(
    Wloc.global, Wloc.ind,
    pchisq(Wloc.global, p * (q - 1), lower.tail = FALSE),
    pchisq(Wloc.ind,    q - 1,       lower.tail = FALSE))
  TestResult$Wconc <- make_test_mat(
    Wconc.global, Wconc.ind,
    pchisq(Wconc.global, p, lower.tail = FALSE),
    pchisq(Wconc.ind,    1, lower.tail = FALSE))
  TestResult$CI.Wconc <- CI.Wconc
  TestResult$s0 <- sapply(proj_jac_list, `[[`, "s0")
  
  # 6. 신뢰 영역 포함 여부 판정 ----
  conf.inside <- list(wald = character(p), Wconc = character(p),
                      Wloc = character(p))
  
  if (!is.null(mu_true)) {
    mu_true   <- as.vector(mu_true)
    mu_unit_0 <- mu_true / sqrt(sum(mu_true^2))
    P_mu_0    <- tcrossprod(mu_unit_0)
    P_perp_0  <- diag(q) - P_mu_0
  } else {
    mu_unit_0 <- mu_unit
    P_perp_0  <- projmat.perp
  }
  cv_wald <- qchisq(conf.level, q)
  cv_loc  <- qchisq(conf.level, q - 1)
  
  for (j in seq_len(p)) {
    beta_j <- beta_new[j + 1, ]
    idx_j  <- idx.list[[j + 1]]
    diff_w <- beta_j - beta_true[j, ]
    S_j    <- cov.mat[idx_j, idx_j]
    conf.inside$wald[j] <- ifelse(
      drop(crossprod(diff_w, sym_inv(S_j) %*% diff_w)) <= cv_wald,
      "inside", "outside")
    
    if (reduce_jacobian) {
      U_hat     <- proj_jac_list[[j]]$U_perp
      perp_true <- drop(P_perp_0 %*% beta_true[j, ])
      v <- drop(crossprod(U_hat, beta_j - perp_true))
    } else {
      v <- drop(projmat.perp %*% beta_j - P_perp_0 %*% beta_true[j, ])
    }
    conf.inside$Wloc[j] <- ifelse(qform_loc(j, v) <= cv_loc,
                                  "inside", "outside")
    
    s_true <- sum(mu_unit_0 * beta_true[j, ])
    conf.inside$Wconc[j] <- ifelse(
      abs(s_hat_vec[j] - s_true) <= z_a2 * se_conc[j],
      "inside", "outside")
  }
  
  # 7. 반환 (원본과 동일 구조) ----
  list(
    beta            = beta_new,
    norm            = apply(beta_new, 1, function(z) sqrt(sum(z^2))),
    Fisher          = list(Fn = Fn, Fnj = fit.Fisher$Fnj),
    orthogonal      = orthogonal,
    wald            = Wald.j.list,
    df              = df,
    idx.list        = idx.list,
    cov.mat         = cov.mat,
    cov.mu          = Omega_conc_n_list,
    cov.mu.full     = Omega_conc_full_n_list,
    cov.ortho       = Omega_loc_n_list,
    cov.ortho.full  = Omega_loc_full_n_list,
    projmat.mu      = projmat.mu,
    projmat.ortho   = projmat.perp,
    s0              = TestResult$s0,
    test            = TestResult,
    conf.inside     = as.data.frame(conf.inside),
    mu_true_used    = !is.null(mu_true),
    U_perp          = lapply(proj_jac_list, `[[`, "U_perp"),
    reduce_jacobian = reduce_jacobian
  )
}





vmf_moments <- function(Theta, q = ncol(Theta), tol = 1e-12) {
  r <- sqrt(rowSums(Theta^2))
  safe_r <- pmax(r, tol)
  S <- Theta / safe_r                       # n x q, 단위 방향
  A <- Aq_vec(r, q)
  a <- ifelse(r < tol, 1/q, A / safe_r)     # A/r  ->  1/q (r -> 0)
  h <- 1 - A^2 - (q - 1) * a                # H_q(r), r -> 0에서 1/q
  list(r = r, S = S, A = A, c1 = a, c2 = h - a)
}
