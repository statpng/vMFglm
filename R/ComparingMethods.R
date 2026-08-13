
#' @export .get_Y_euc 
.get_Y_euc <- function(Y, X1) {
  lv1 <- names(table(X1))[1]
  lv2 <- names(table(X1))[2]
  euc   <- euclideanization(t(Y[X1 == lv1, ]), t(Y[X1 == lv2, ]))
  Y_euc <- matrix(NA, nrow(Y), ncol(Y) - 1)
  Y_euc[X1 == lv1, ] <- euc$euclideanG1
  Y_euc[X1 == lv2, ] <- euc$euclideanG2
  Y_euc
}

.vmf_loglik <- function(Y, mu, kappa) {
  p <- ncol(Y); nu <- p / 2 - 1
  logI  <- log(besselI(kappa, nu = nu, expon.scaled = TRUE)) + kappa
  log_c <- nu * log(kappa) - (p / 2) * log(2 * pi) - logI
  nrow(Y) * log_c + kappa * sum(Y %*% mu)
}

.kappa_mle <- function(Rbar, p, tol = 1e-8, maxit = 100) {
  if (Rbar <= 1e-8) return(1e-6)
  if (Rbar >= 1 - 1e-8) return(1e6)
  kappa <- Rbar * (p - Rbar^2) / (1 - Rbar^2); nu <- p / 2 - 1
  for (i in seq_len(maxit)) {
    Ak  <- besselI(kappa, nu + 1, expon.scaled = TRUE) /
      besselI(kappa, nu,     expon.scaled = TRUE)
    dAk <- 1 - Ak^2 - (p - 1) / kappa * Ak
    if (abs(dAk) < 1e-14) break
    delta <- (Ak - Rbar) / dAk
    kappa <- max(1e-6, kappa - delta)
    if (abs(delta) < tol) break
  }
  kappa
}


comp_VMF_LRT_inner <- function(Y, X1) {
  
  lv <- names(table(X1))
  
  Y0 <- Y[X1 == lv[1], ]; Y1 <- Y[X1 == lv[2], ]
  p  <- ncol(Y); n0 <- nrow(Y0); n1 <- nrow(Y1)
  mle0 <- Directional::vmf.mle(Y0); mle1 <- Directional::vmf.mle(Y1)
  ll_full <- .vmf_loglik(Y0, mle0$mu, mle0$kappa) +
    .vmf_loglik(Y1, mle1$mu, mle1$kappa)
  
  mu_pool  <- colSums(Y0) + colSums(Y1)
  mu_pool  <- mu_pool / sqrt(sum(mu_pool^2))
  kappa0_h <- .kappa_mle(mean(Y0 %*% mu_pool), p)
  kappa1_h <- .kappa_mle(mean(Y1 %*% mu_pool), p)
  ll_h0_loc <- .vmf_loglik(Y0, mu_pool, kappa0_h) +
    .vmf_loglik(Y1, mu_pool, kappa1_h)
  p_loc <- pchisq(max(0, 2*(ll_full - ll_h0_loc)), df = p-1, lower.tail = FALSE)
  
  Ybar0 <- colMeans(Y0); mu0_f <- Ybar0 / sqrt(sum(Ybar0^2))
  Ybar1 <- colMeans(Y1); mu1_f <- Ybar1 / sqrt(sum(Ybar1^2))
  Rbar_pool  <- (n0*sqrt(sum(Ybar0^2)) + n1*sqrt(sum(Ybar1^2))) / (n0+n1)
  kappa_pool <- .kappa_mle(Rbar_pool, p)
  ll_h0_conc <- .vmf_loglik(Y0, mu0_f, kappa_pool) +
    .vmf_loglik(Y1, mu1_f, kappa_pool)
  p_conc <- pchisq(max(0, 2*(ll_full - ll_h0_conc)), df = 1, lower.tail = FALSE)
  
  c(p_loc_VMF = p_loc, p_conc_VMF = p_conc)
}

#' @export comp_VMF_LRT
comp_VMF_LRT <- function(Y, X1, X2=NULL, degree = 1, ortho=FALSE) {
  
  if(!is.null(X2)){
    if(ortho){
      X1_num <- as.numeric(X1 == names(table(X1))[2])
      X2_orth <- residuals(lm(X2 ~ X1_num))   # X1과 직교하는 X2 성분만
    } else {
      X2_orth <- X2
    }
    
    pp     <- .geod_preprocess(Y, X2_orth, degree, tol=1e-6, maxit=500)  # ← X2_orth 사용
    Y_loc  <- pp$Y_adj_loc
    Y_adj <- pp$Y_adj
    
  } else {
    
    Y_adj <- Y
    
  }
  
  comp_VMF_LRT_inner(Y_adj, X1)
}



log_map <- function(Y, m) {
  ip <- pmin(pmax(as.numeric(Y %*% m), -1), 1)
  th <- acos(ip)
  V  <- Y - outer(ip, m)
  nv <- sqrt(rowSums(V^2)); nv[nv < 1e-12] <- 1
  V * (th / nv)
}


#' @export comp_TanSpace
comp_TanSpace <- function(Y, X1, X2=NULL, type="binary", degree = 1, ortho=FALSE) {
  
  m     <- colMeans(Y); m <- m / sqrt(sum(m^2))
  Bmat  <- svd(diag(length(m)) - outer(m, m))$u[, 1:(length(m)-1), drop = FALSE]
  Y_euc <- log_map(Y, m) %*% Bmat
  
  # Y_euc <- Directional::euclid.inv(Y)
  
  
  
  if(!is.null(X2)){
    
    if(ortho){
      X1_num <- as.numeric(X1 == names(table(X1))[2])
      X2_orth <- residuals(lm(X2 ~ X1_num))   # X1과 직교하는 X2 성분만
    } else {
      X2_orth <- X2
    }
    
    
    # Step 1: location 제거 (X2_orth 사용)
    Y_adj <- as.matrix( residuals(lm(Y_euc ~ poly(X2_orth, degree))) )
    
    # Step 2: dispersion 제거 (X2_orth 사용)
    d_i       <- sqrt(rowSums(Y_adj^2))
    d_pred    <- exp(fitted(lm(log(d_i) ~ poly(X2_orth, degree))))
    Y_student <- Y_adj / d_pred
    
  } else {
    
    Y_student <- Y_euc
    
  }
  
  if(type == "binary"){
    
    p_loc  <- HotellingStat(Y=Y_student, group=X1, euclidean=FALSE)$pvalue
    bd     <- vegan::betadisper(dist(Y_student), group=as.factor(X1))
    p_conc <- anova(bd)$`Pr(>F)`[1]
    
  } else {
    
    man    <- summary(manova(Y_student ~ X1), test = "Pillai")
    p_loc  <- man$stats["X1", "Pr(>F)"]
    
    mu_hat <- fitted(lm(Y_student ~ X1))
    res    <- Y_student - mu_hat
    d_test <- sqrt(rowSums(res^2))
    p_conc <- anova(lm(log(d_test) ~ X1))$`Pr(>F)`[1]
    
  }
  
  c(p_loc_TS = p_loc, p_conc_TS = p_conc)
}



#' @import dglm
#' @export comp_dglm
comp_dglm <- function(Y, X1, X2=NULL) {
  library(dglm)
  
  m     <- colMeans(Y); m <- m / sqrt(sum(m^2))
  Bmat  <- svd(diag(length(m)) - outer(m, m))$u[, 1:(length(m)-1), drop = FALSE]
  Y_euc <- log_map(Y, m) %*% Bmat
  
  # Y_euc <- Directional::euclid.inv(Y)
  
  
  if(!is.null(X2)){
    df <- data.frame(X1=X1, X2=X2)
    
    dglm.mat <- apply(Y_euc, 2, function(y){
      fit=dglm(y ~ X1+X2, ~X1+X2, family=gaussian, data=df)
      c(loc=summary(fit)$coefficients[-1,4],
        disp=summary(fit)$dispersion.summary$coefficients[-1,4])
    })
    
  } else {
    df <- data.frame(X1=X1)
    
    dglm.mat <- apply(Y_euc, 2, function(y){
      fit=dglm(y ~ X1, ~X1, family=gaussian, data=df)
      c(loc=summary(fit)$coefficients[-1,4],
        disp=summary(fit)$dispersion.summary$coefficients[-1,4])
    })
    
  }
  
  dglm.comb <- apply(dglm.mat, 1, function(pv) 1-pchisq( sum(-2*log(na.omit(pv))), 2*length(na.omit(pv)) ) )
  
  if(is.null(X2)){
    names(dglm.comb) <- c("loc.X1", "disp.X1")
  }
  # dglm.comb <- dglm.comb[c("loc.X1", "disp.X1")]
  # names(dglm.comb) <- c("p_loc_dglm", "p_conc_dglm")
  
  dglm.comb
}


#' @export comp_dglm2
comp_dglm2 <- function(Y, X1, X2=NULL) {
  library(dglm)
  
  m     <- colMeans(Y); m <- m / sqrt(sum(m^2))
  Bmat  <- svd(diag(length(m)) - outer(m, m))$u[, 1:(length(m)-1), drop = FALSE]
  Y_euc <- log_map(Y, m) %*% Bmat
  
  # Y_euc <- Directional::euclid.inv(Y)
  
  
  if(!is.null(X2)){
    df <- data.frame(X1=X1, X2=X2)
    
  } else {
    df <- data.frame(X1=X1)
  }
  
  dglm.mat <- apply(Y_euc, 2, function(y){
    fit=dglm(y ~ ., ~ ., family=gaussian, data=df)
    out <- c(loc=summary(fit)$coefficients[-1,4],
             disp=summary(fit)$dispersion.summary$coefficients[-1,4])
    names(out) <- c(paste0("loc.X",1:ncol(df)), paste0("disp.X",1:ncol(df)))
    out
  })
  
  
  dglm.comb <- apply(dglm.mat, 1, function(pv) 1-pchisq( sum(-2*log(na.omit(pv))), 2*length(na.omit(pv)) ) )
  dglm.comb <- dglm.comb[c("loc.X1", "disp.X1")]
  names(dglm.comb) <- c("p_loc_dglm", "p_conc_dglm")
  dglm.comb
}




# ── 구면 기하 헬퍼 ─────────────────────────────────────────────
.log_base <- function(base, Y) {
  ct <- pmin(pmax(drop(Y %*% base), -1 + 1e-8), 1 - 1e-8)
  th <- acos(ct); V <- Y - outer(ct, base)
  V * (th / sin(th))
}
.exp_base <- function(base, V) {
  nrm <- sqrt(rowSums(V^2)); s <- pmax(nrm, 1e-10)
  Y <- cos(nrm) * matrix(base, nrow(V), length(base), byrow = TRUE) +
    (sin(s)/s) * V
  Y / sqrt(rowSums(Y^2))
}
.log_pw <- function(Y_hat, Y) {
  ct <- pmin(pmax(rowSums(Y * Y_hat), -1 + 1e-8), 1 - 1e-8)
  th <- acos(ct); V <- Y - ct * Y_hat
  V * (th / sin(th))
}
.par_trans <- function(p_mat, q_base, V) {
  n <- nrow(p_mat); pd <- ncol(p_mat)
  ct <- pmin(pmax(drop(p_mat %*% q_base), -1 + 1e-8), 1 - 1e-8)
  d <- acos(ct); s <- pmax(sin(d), 1e-10)
  Q <- matrix(q_base, n, pd, byrow = TRUE)
  u_p <- (Q - ct * p_mat) / s
  vdu <- rowSums(V * u_p)
  V - sin(d) * vdu * p_mat - (1 - ct) * vdu * u_p
}


# ── 측지 회귀 / 전처리 (원래 코드 그대로) ──────────────────────
#' @export .geod_reg
.geod_reg <- function(Y_sph, X2, degree = 1, tol = 1e-6, maxit = 300) {
  n <- nrow(Y_sph); p <- ncol(Y_sph)
  X2p <- poly(X2, degree)
  K   <- ncol(X2p)                    # ← (1) 실제 basis 열 수
  mu0 <- colMeans(Y_sph); mu0 <- mu0 / sqrt(sum(mu0^2))
  beta <- matrix(0, p, K)             # ← (2) degree → K
  lr <- 0.1; prev_loss <- Inf
  
  for (iter in seq_len(maxit)) {
    eta   <- X2p %*% t(beta)
    Y_hat <- .exp_base(mu0, eta)
    V_res <- .log_pw(Y_hat, Y_sph)
    loss  <- sum(rowSums(V_res^2))
    if (abs(prev_loss - loss) < tol * (1 + abs(prev_loss))) break
    prev_loss <- loss
    
    V_mu0 <- .par_trans(Y_hat, mu0, V_res)
    grad_loss <- -2 * t(crossprod(X2p, V_mu0))
    
    beta_new <- beta - lr * grad_loss
    for (k in seq_len(K)){             # ← (3) degree → K
      beta_new[, k] <- beta_new[, k] - sum(beta_new[, k] * mu0) * mu0
    }
    
    eta_new   <- X2p %*% t(beta_new)
    Y_hat_new <- .exp_base(mu0, eta_new)
    loss_new  <- sum(rowSums(.log_pw(Y_hat_new, Y_sph)^2))
    
    if (loss_new < loss) { 
      beta <- beta_new; lr <- lr * 1.1 
    } else { 
      lr <- lr * 0.5; if (lr < 1e-14) break 
    }
  }
  eta   <- X2p %*% t(beta)
  Y_hat <- .exp_base(mu0, eta)
  V_res <- .log_pw(Y_hat, Y_sph)
  V_mu0 <- .par_trans(Y_hat, mu0, V_res)
  list(mu0 = mu0, beta = beta, Y_hat = Y_hat, V_mu0 = V_mu0,
       iter = iter, loss = loss)
}

.geod_preprocess <- function(Y_sph, X2, degree = 1, tol = 1e-6, maxit = 300) {
  gr <- .geod_reg(Y_sph, X2, degree, tol, maxit)
  V_mu0 <- gr$V_mu0
  Y_adj_loc <- .exp_base(gr$mu0, V_mu0)
  d_i    <- pmax(sqrt(rowSums(V_mu0^2)), 1e-10)
  d_pred <- exp(fitted(lm(log(d_i) ~ poly(X2, degree))))
  V_std  <- V_mu0 / d_pred
  Y_adj_conc  <- .exp_base(gr$mu0, V_std)
  list(Y_adj_loc = Y_adj_loc, 
       Y_adj = Y_adj_conc, 
       mu0 = gr$mu0, gr = gr)
}







.safe <- function(expr) {
  out <- tryCatch(expr,
                  error   = function(e) NA_real_,
                  warning = function(w) suppressWarnings(expr))
  if (length(out) == 0 || is.null(out)) NA_real_ else out
}

## 리스트/객체에서 loglik 을 robust 하게 추출
## (Directional 회귀들은 버전에 따라 loglik / logLik / lik / logl  로 반환)
.get_loglik <- function(fit) {
  if (is.null(fit)) return(NA_real_)
  cand <- c("loglik", "logLik", "lik", "loglik1", "logl", "ll", "value")
  for (nm in cand) {
    if (!is.null(fit[[nm]])) {
      v <- as.numeric(fit[[nm]])
      return(v[length(v)])   ## 벡터면 마지막(수렴) 값
    }
  }
  NA_real_
}

## nested-model LRT:  H0 reduced  vs  H1 full
.nested_lrt <- function(ll_full, ll_red, df) {
  if (is.na(ll_full) || is.na(ll_red) || df <= 0) return(NA_real_)
  stat <- 2 * (ll_full - ll_red)
  if (!is.finite(stat)) return(NA_real_)
  if (stat < 0) stat <- 0
  pchisq(stat, df = df, lower.tail = FALSE)
}


## (1) Regression 비교 : vMF / IAG / ESAG  (GLOBAL effect, LRT)
##
##   Y      : n x q 단위벡터 행렬
##   X1, X2 : 길이 n covariate 벡터  (X1 = group indicator, X2 = nuisance)
##   q = ncol(Y)
##
##   주의: 이 함수들은 mu 벡터 전체를 covariate 에 연결하므로
##         location + concentration 을 합친 GLOBAL 효과만 검정한다.
##         반환 이름 p_glob_* 는 vMF-GLM 의 wald 행과 정렬된다.

comp_Directional_reg <- function(Y, X1, X2) {
  
  q <- ncol(Y)
  
  Xfull <- cbind(X1, X2)     ## H1 : X1 + X2
  Xred  <- cbind(X2)         ## H0 : X2 only  (X1 효과 제거)
  
  ## ---- vMF regression (vmfreg) : H0 자유도 = q -------------------
  start <- proc.time()
  p_glob_vMFreg <- .safe({
    f_full <- Directional::vmfreg(Y, Xfull)
    f_red  <- Directional::vmfreg(Y, Xred)
    .nested_lrt(.get_loglik(f_full), .get_loglik(f_red), df = q)
  })
  end <- proc.time()
  time.vMFreg <- end-start
  
  ## ---- IAG (projected normal) regression : 자유도 = q -----------
  start <- proc.time()
  p_glob_IAG <- .safe({
    f_full <- Directional::iag.reg(Y, Xfull)
    f_red  <- Directional::iag.reg(Y, Xred)
    .nested_lrt(.get_loglik(f_full), .get_loglik(f_red), df = q)
  })
  end <- proc.time()
  time.IAG <- end-start
  
  start <- proc.time()
  p_glob_SPCauchy <- .safe({
    f_full <- Directional::spcauchy.reg (Y, Xfull)
    f_red  <- Directional::spcauchy.reg (Y, Xred)
    .nested_lrt(.get_loglik(f_full), .get_loglik(f_red), df = q)
  })
  end <- proc.time()
  time.SPCauchy <- end-start
  
  ## ---- ESAG (Structure 2, anisotropic) : 자유도 = q + 2 ---------
  ## rotation(Q) 그리드 탐색을 생략한다:
  ##   - 시뮬레이션 데이터는 좌표계가 고정되어 있고,
  ##   - full/reduced 가 동일 좌표계를 공유해야 nested LRT 가 정직하다.
  ##     (esag.reg 를 그대로 부르면 full/reduced 가 서로 다른 Q 를
  ##      독립 최적화하여 l_red > l_full 같은 비정합이 생길 수 있음)
  ## Directional 내부의 .esag.reg2 를 Q=I 로 직접 호출한다.
  ## 버전에 따라 .esag.reg2 가 없으면 esag.reg(lati=1,longi=1) 로 대체.
  start <- proc.time()
  p_glob_ESAG <- .safe({
    Xf <- model.matrix(~., data.frame(Xfull))   ## intercept 포함
    Xr <- model.matrix(~., data.frame(Xred))
    
    f_full <- Directional:::.esag.reg2(Y, Xf, con = TRUE, xnew = NULL, tol = 1e-6)
    f_red  <- Directional:::.esag.reg2(Y, Xr, con = TRUE, xnew = NULL, tol = 1e-6)
    .nested_lrt(.get_loglik(f_full), .get_loglik(f_red), df = q + 2)
  })
  end <- proc.time()
  time.ESAG <- end-start
  
  ## ---- spherical Cauchy (spcauchy.reg) : 자유도 = q -------------
  ## heavy-tail, isotropic(rotationally symmetric) → mu(x)=Bx 만 회귀.
  ## be 행렬이 q열 한 블록뿐(anisotropy 모수 없음)이므로 X1 제거 시
  ## 사라지는 자유 모수는 q개 → df = q. rotation 탐색 불필요.
  
  start <- proc.time()
  p_glob_PKBD <- .safe({
    f_full <- Directional::pkbd.reg(Y, Xfull)
    f_red  <- Directional::pkbd.reg(Y, Xred)
    
    .nested_lrt(.get_loglik(f_full), .get_loglik(f_red), df = q)
  })
  end <- proc.time()
  time.PKBD <- end-start
  
  
  list(
    result = c(p_glob_VMFreg = p_glob_vMFreg,
               p_glob_IAG    = p_glob_IAG,
               p_glob_ESAG   = p_glob_ESAG,
               p_glob_SPCauchy   = p_glob_SPCauchy,
               p_glob_PKBD   = p_glob_PKBD),
    timing = c(vMFreg=time.vMFreg[3], 
               ESAG=time.ESAG[3],
               SPCauchy=time.SPCauchy[3], 
               IAG=time.IAG[3], 
               PKBD=time.PKBD[3])
  )
  
  
}


## (2) Location & concentration two-sample 검정
##
##   X1 부호로 두 그룹 분리, X2 효과는 studentize_fun 으로 사전 제거.
##   studentize_fun : (Y, X1, X2) -> Y*(n x q sphere). NULL 이면 raw Y.

studentize_fun <- function(Y, X1, X2=NULL, degree = 1, ortho=FALSE){
  
  ## X2 가 없으면 보정 없이 원본 반환 (comp_VMF_LRT 의 else 분기와 동일)
  if (missing(X2) || is.null(X2)) return(Y)
  
  if (ortho) {
    ## X1 을 0/1 수치형으로: 두 번째 레벨을 1 로
    lev    <- names(table(X1))
    X1_num <- as.numeric(X1 == lev[2])
    X2_orth <- residuals(lm(X2 ~ X1_num))   ## X1 과 직교하는 X2 성분
  } else {
    X2_orth <- X2
  }
  
  pp <- .geod_preprocess(Y, X2_orth, degree, tol = 1e-6, maxit = 500)
  pp$Y_adj   ## comp_VMF_LRT_inner 에 들어가는 스튜던트화 구면 잔차
}




comp_Directional_test <- function(Y, X1, X2=NULL, B = 999) {
  
  if( !is.null(X2) ){
    Ystar <- studentize_fun(Y, X1, X2)
  } else {
    Ystar <- Y
  }
  
  grp  <- ifelse(X1 >= 0, 1L, 2L)
  X1.levels <- names(table(X1))
  x1 <- Ystar[X1==X1.levels[1], , drop = FALSE]
  x2 <- Ystar[X1==X1.levels[2], , drop = FALSE]
  
  ## location : vMF LRT (lr.perm)
  start <- proc.time()
  p_loc_LRperm <- .safe({
    as.numeric(Directional::lr.perm(x1, x2, B = B)$p.value)
  })
  end <- proc.time()
  time.LRperm <- end-start
  
  ## location : 비등농도 LRT (het.perm)
  start <- proc.time()
  p_loc_HETperm <- .safe({
    as.numeric(Directional::het.perm(x1, x2, B = B)$p.value)
  })
  end <- proc.time()
  time.HETperm <- end-start
  
  ## concentration : 그룹 간 농도 동일성 (spherconc.test)
  start <- proc.time()
  p_conc_SPHER <- .safe({
    as.numeric(Directional::spherconc.test(Ystar, ina = grp)$res["p-value"])
  })
  end <- proc.time()
  time.SPHER <- end-start
  
  list(result = c(p_loc_LRperm  = p_loc_LRperm,
                  p_loc_HETperm = p_loc_HETperm,
                  p_conc_SPHER  = p_conc_SPHER),
       timing = c(LRperm=time.LRperm[3],
                  HETperm=time.HETperm[3],
                  SPHER=time.SPHER[3]))
}



