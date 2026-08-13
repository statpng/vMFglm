#' @importFrom manifold frechetMean
#' @importFrom manifold createM
#' @importFrom magic adiag
#' @importFrom dplyr `%>%`
#' @useDynLib vMFglm
#' 
#' @export vMFglm
# vMFglm <- function(X, Y, MU=NULL, orthogonal=FALSE, penalty.factor=rep(0,ncol(X)), standardize=TRUE, Offset=NULL, maxit=100, eps=1e-6, lambda=1e-3, gamma=NULL){
#   # 250206
#   
#   if(FALSE){
#     MU=NULL; orthogonal=FALSE; standardize=TRUE; Offset=NULL; maxit=100; eps=1e-4; lambda=1e-4; gamma=0; penalty.factor=rep(0,ncol(X))
#   }
#   
#   if(FALSE){
#     
#     devtools::load_all()
#     
#     set.seed(1)
#     simdata <- sim.vMFglm(n=100, p=2, q=3, mu=c(0,0,20), s=10, s0=0, type="vMF", seed.UDV=1, qr.U=FALSE, qr.V=TRUE)
#     
#     X=simdata$X; Y=simdata$Y; orthogonal=TRUE; standardize=TRUE; Offset=NULL; maxit=100; eps=1e-4; lambda=1e-4; penalty.factor=rep(0,ncol(X))
#   }
#   
# 
#   
#   if(FALSE){
#     devtools::load_all()
#     
#     
#     library(dplyr)
#     set.seed(1)
#     simdata <- vMFglm::sim.vMFglm(n=100, p=1, q=3, mu=c(0,10,0), snr=5, s=5, s0=0.0, type=c("vMF", "Proj", "ExpMap")[1], dir = "ortho")
#     
#     
#     system.time(fit0 <- vMFglm(X=simdata$X, Y=simdata$Y, eps=1e-6, maxit=100, lambda=1e-10, orthogonal=F))
#     system.time(fit1 <- vMFglm.R(X=simdata$X, Y=simdata$Y, eps=1e-6, maxit=100, lambda=1e-6, orthogonal=F))
#     
#     fit0
#     
#     fit0 %>% class
#     fit0 %>% summary
#     fit0 %>% plot
#     
#     fit1 %>% class
#     fit1 %>% summary
#     fit1 %>% plot
#     
#     
#     
#     #
#     
#     
#     fit0 %>% summary()
#     fit1 %>% summary()
#     
#     fit0$beta
#     fit1$beta
#     
#     plot(fit0, barblen = 0.001, lwd = 0.1, cex.tangent = 0.5)
#     plot(fit1)
#     
#     
#     
#     res.benchmark <- microbenchmark::microbenchmark(
#       fit0=vMFglm(X=simdata$X, Y=simdata$Y, eps=1e-6, maxit=100, lambda=5e-4, orthogonal=T),
#       fit1=vMFglm.R(X=simdata$X, Y=simdata$Y, eps=1e-6, maxit=100, lambda=5e-4, orthogonal=T), times=10
#     )
#     res.benchmark
#     
#     
#     
#     
#     
#     set.seed(1)
#     simdata <- vMFglm::sim.vMFglm(n=1000, p=1, q=3, mu=c(0,100,0), snr=50, s=10, s0=0.0, type=c("vMF", "Proj", "ExpMap")[1])
#     
#     res.benchmark.1000 <- microbenchmark::microbenchmark(
#       fit0=vMFglm_cpp(X=simdata$X, Y=simdata$Y, eps=1e-6, maxit=100, lambda=1e-6, orthogonal=F), times=10
#     )
#     
#     res.benchmark.1000
#     #
#     
#     
#     #
#     
#   }
#   
#   
#   
#   
#   
#   X <- as.matrix(X)
#   Y <- as.matrix(Y)
#   
#   
#   
#   X00 <- X
#   X0 <- scale(X, center=TRUE, scale=FALSE)
#   
#   if(standardize){
#     sdx.inv <- apply(X0,2,sd) %>% {diag(1/., length(.), length(.))}
#     X <- X0 %*% sdx.inv
#   } else {
#     X <- X0
#   }
#   
#   
#   n <- nrow(X)
#   p <- ncol(X)
#   q <- ncol(Y)
#   
# 
#   if(is.null(Offset)){
#     Offset <- matrix(0,n,q)
#   }
#   
#   
#   
#   {
#     
#     
#     MU <- FrechetMean(Y)
#     
#     muhat <- colSums(Y) %>% {./norm(.,"2")} # mle of mean direction in vMF
#     rbar <- colSums(Y) %>% {norm(.,"2")/nrow(Y)}
#     
#     kappa.type <- 4
#     
#     if(is.null(kappa.type)){
#       kappa0 <- 1
#     } else {
#       kappa0 <- list(
#         (q-1) / 2*(1-rbar),
#         q*rbar*( 1+q/(q+2)*rbar^2 + q^2*(q+8)/((q+2)^2*(q+4))*rbar^4 ),
#         (rbar*q) / (1-rbar^2),
#         (rbar*q - rbar^3) / (1-rbar^2)
#       )[[kappa.type]]
#     }
#     
#     
#     # kappa0 <- (q-1) / 2*(1-rbar) # approximation
#     # kappa0 <- q*rbar*( 1+q/(q+2)*rbar^2 + q^2*(q+8)/((q+2)^2*(q+4))*rbar^4 ) # approximation
#     # kappa0 <- (rbar*q) / (1-rbar^2) # approximation
#     # kappa0 <- (rbar*q - rbar^3) / (1-rbar^2) # approximation
#     
#     mu <- muhat * kappa0
#     
#     
#     
#     
#     
#     
#     {
#       B0 <- lm(Y ~ -1 + X, offset=tcrossprod(rep(1,n),MU)+Offset) %>% coef
#       if(orthogonal){
#         B0 <- project2orthogonal(B0, mu)
#       }
#       
#       beta <- as.vector(t(B0))
#       # beta <- as.vector(t(rbind(mu, B0)))
#     }
#     
#     
#     
#     
#     
#     l1 <- 1
#     beta_old <- beta
#     beta_new <- beta + 1
#     beta.list <- Fn.list <- loglik.list <- crit.list <- NULL
#     crit <- 1
#     
#     
#     if(orthogonal){
#       
#       if( !is.null(gamma) ){
#         
#         if(length(gamma)==1){
#           gamma <- rep(gamma, p)
#         }
#         
#         X1 <- cbind(1, X)
#         Xt_list <- apply(X1, 1, function(Xi) kronecker(Xi, diag(1,q)), simplify=FALSE) # [pq x q]
#         Xt <- do.call("cbind", Xt.list) # pq x qn
#         yt <- as.vector(t(Y))
#         beta <- as.vector(t(rbind(mu, B0)))
#         
#         fit0 <- suppressWarnings({
#           vMFglm_iteration_FixedGamma(X=X, Y=Y, Offset=Offset, beta=beta, Xt=Xt, Xt_list=Xt_list, eps=eps, maxit=maxit, lambda=lambda, orthogonal=orthogonal, gamma=gamma, zero_beta = which(penalty.factor==1))
#         })
#         
#       } else {
#         
#         X1 <- X
#         Xt.list <- apply(X1, 1, function(Xi) kronecker(Xi, diag(1,q)), simplify=FALSE) # [pq x q]
#         Xt <- do.call("cbind", Xt.list) # pq x qn
#         yt <- as.vector(t(Y))
#         
#         fit0 <- suppressWarnings({
#           vMFglm_iteration_KKT(X=X, Y=Y, Offset=Offset, mu=mu, beta=beta, Xt=Xt, Xt_list=Xt.list, eps=eps, maxit=maxit, lambda=lambda)
#         })
#         
#       }
#       
#     } else {
#       
#       gamma <- 0
#       if(length(gamma)==1){
#         gamma <- rep(gamma, p)
#       }
#       
#       X1 <- cbind(1, X)
#       Xt.list <- apply(X1, 1, function(Xi) kronecker(Xi, diag(1,q)), simplify=FALSE) # [pq x q]
#       Xt <- do.call("cbind", Xt.list) # pq x qn
#       yt <- as.vector(t(Y))
#       beta <- as.vector(t(rbind(mu, B0)))
#       
#       if( !is.null(penalty.factor) ){
#         fit0 <- vMFglm_iteration_FixedGamma(X=X, Y=Y, Offset=Offset, beta=beta, Xt=Xt, Xt_list=Xt.list, eps=eps, maxit=maxit, lambda=lambda, orthogonal=orthogonal, gamma=gamma, zero_beta = which(penalty.factor==1))
#       } else {
#         fit0 <- vMFglm_iteration_FixedGamma(X=X, Y=Y, Offset=Offset, beta=beta, Xt=Xt, Xt_list=Xt.list, eps=eps, maxit=maxit, lambda=lambda, orthogonal=orthogonal, gamma=gamma, zero_beta = NULL)
#       }
#       
#       
#     }
#     
#     #
#     
#     
#     beta.list <- fit0$beta.list
#     beta_new <- fit0$beta
#     gamma <- fit0$gamma
#     
#     beta.list <- lapply( beta.list, function(b){
#       matrix(ifelse(abs(b) < 1e-10, 0, b), p+1, q, byrow=TRUE)
#     })
#     beta_new <- ifelse(abs(beta_new) < 1e-10, 0, beta_new)
#     beta_new <- matrix(beta_new, p+1, q, byrow=TRUE)
#     
#     if(standardize){
#       beta.list <- lapply(beta.list, function(b) as.matrix( Matrix::bdiag(1,sdx.inv) %*% b) )
#       beta_new <- as.matrix( Matrix::bdiag(1,sdx.inv) %*% beta_new )
#     }
#     
#     
#     loglik.list <- fit0$loglik.list
#     crit.list <- fit0$crit.list
#     Fn.list <- fit0$Fn.list
#     
#     iterations <- fit0$iterations
#     
#     
# 
#     
#     
#     result <- list(X = X0,
#                    Y = Y,
#                    gamma = gamma,
#                    mu = beta_new[1,],
#                    beta = beta_new,
#                    beta.list = beta.list,
#                    offset = Offset,
#                    loglik.list = fit0$loglik.list[-1], 
#                    crit.list = fit0$crit.list[-1],
#                    n.iter = fit0$iterations-1)
#   }
#   
#   
#   
#   
#   params <- list( MU=muhat * kappa0, Offset=Offset, maxit=maxit, eps=eps, standardize=standardize, orthogonal=orthogonal, lambda=lambda )
#   
#   result$params <- params
#   
#   class(result) <- "vMFglm"
#   
#   return(result)
#   
# }






# ============================================================
#  vMFglm_main.R : 메인 함수 개선판 (drop-in 대체)
#
#  새 인자
#   method = c("newton", "irls")
#     "newton" (기본): 벌점 Fisher scoring + step-halving.
#        고집중(kappa 큰) 데이터에서도 집중도 효과(s_conc)를 올바르게 추정.
#     "irls": 기존 알고리즘 재현 (기존 결과와의 비교/재현용).
#   engine = c("auto", "cpp", "R")
#     "cpp": vMFglm_iteration_FixedGamma_fast (vMFglm_fast.cpp) 사용.
#     "R"  : 순수 R 경로 (vmf_fast.R) 사용. 컴파일 불가 환경 대비.
#     "auto"(기본): cpp 함수가 로드되어 있으면 cpp, 아니면 R.
#   kappa_max: newton에서 ||beta|| 상한. 초과 시 경고 후 중단 (MLE 발산 신호).
#
#  원본 대비 수정
#   - orthogonal+gamma 분기의 Xt_list/Xt.list 변수명 불일치 버그 제거
#     (fast 반복은 Xt/Xt_list 자체가 불필요)
#   - Xt (pq x nq) 및 Xt.list를 아예 만들지 않음 (메모리 절약)
#   - 사용자가 MU를 지정하면 무시하지 않고 초기값에 사용
#   - 반환값에 converged / status / method 추가
#
#  의존: vMFglm_fast.cpp (engine="cpp") 또는 vmf_fast.R (engine="R")
#  orthogonal & is.null(gamma) 의 KKT 경로는 기존
#  vMFglm_iteration_KKT 가 로드되어 있을 때만 그대로 위임.
# ============================================================

#' @export vMFglm
vMFglm <- function(X, Y, MU = NULL, orthogonal = FALSE,
                   penalty.factor = rep(0, ncol(X)),
                   standardize = TRUE, Offset = NULL,
                   maxit = 500, eps = 1e-6, lambda = 1e-5, gamma = NULL,
                   method = c("newton", "irls"),
                   engine = c("auto", "cpp", "R"),
                   kappa_max = 1e6) {
  
  method <- match.arg(method)
  engine <- match.arg(engine)
  if (engine == "auto")
    engine <- if (exists("vMFglm_iteration_FixedGamma_fast")) "cpp" else "R"
  
  X <- as.matrix(X); Y <- as.matrix(Y)
  n <- nrow(X); p <- ncol(X); q <- ncol(Y)
  
  # ── 표준화 (원본과 동일) ─────────────────────────────────
  X0 <- scale(X, center = TRUE, scale = FALSE)
  if (standardize) {
    sdx     <- apply(X0, 2, sd)
    sdx.inv <- diag(1 / sdx, p, p)
    Xs      <- X0 %*% sdx.inv
  } else {
    sdx.inv <- diag(1, p, p)
    Xs      <- X0
  }
  if (is.null(Offset)) Offset <- matrix(0, n, q)
  
  # ── 초기값 (원본과 동일한 kappa 근사 + LS) ────────────────
  csY   <- colSums(Y)
  muhat <- csY / sqrt(sum(csY^2))
  rbar  <- sqrt(sum(csY^2)) / n
  kappa0 <- min(9999, (rbar * q - rbar^3) / (1 - rbar^2))
  mu0 <- if (!is.null(MU)) as.vector(MU) else muhat * kappa0
  
  MU_dir <- if (exists("FrechetMean")) tryCatch(FrechetMean(Y), error = function(e) muhat)
  else muhat
  B0 <- qr.coef(qr(Xs), Y - tcrossprod(rep(1, n), MU_dir) - Offset)
  B0[!is.finite(B0)] <- 0
  if (orthogonal && exists("project2orthogonal"))
    B0 <- project2orthogonal(B0, mu0)
  beta_init <- as.vector(t(rbind(mu0, B0)))
  
  zero_beta <- {
    zi <- which(penalty.factor == 1)
    if (length(zi)) as.integer(zi) else NULL
  }
  
  # ── KKT 경로 (orthogonal & gamma 미지정): 기존 함수에 위임 ──
  if (orthogonal && is.null(gamma)) {
    if (!exists("vMFglm_iteration_KKT"))
      stop("orthogonal=TRUE, gamma=NULL 은 vMFglm_iteration_KKT 가 필요합니다. ",
           "gamma 값을 지정하거나 KKT 반복 함수를 로드하세요.")
    X1      <- cbind(1, Xs)  # KKT 원 코드 관례에 맞춰 구성
    Xt.list <- apply(Xs, 1, function(Xi) kronecker(Xi, diag(1, q)), simplify = FALSE)
    Xt      <- do.call(cbind, Xt.list)
    fit0 <- suppressWarnings(
      vMFglm_iteration_KKT(X = Xs, Y = Y, Offset = Offset, mu = mu0,
                           beta = as.vector(t(B0)), Xt = Xt, Xt_list = Xt.list,
                           eps = eps, maxit = maxit, lambda = lambda))
    gamma_out <- fit0$gamma
  } else {
    if (is.null(gamma)) gamma <- 0
    if (length(gamma) == 1) gamma <- rep(gamma, p)
    
    if (engine == "cpp") {
      fit0 <- vMFglm_iteration_FixedGamma_fast(
        X = Xs, Y = Y, Offset = Offset, beta = beta_init,
        Xt = matrix(0, 0, 0), Xt_list = list(),        # 미사용 (호환 인자)
        eps = eps, maxit = maxit, lambda = lambda,
        orthogonal = orthogonal, gamma = gamma,
        zero_beta = zero_beta, verbose = FALSE,
        method = method, kappa_max = kappa_max)
    } else {
      fit0 <- vMFglm_iteration_R(
        Xs = Xs, Y = Y, Offset = Offset, beta_init = beta_init,
        eps = eps, maxit = maxit, lambda = lambda,
        orthogonal = orthogonal, gamma = gamma,
        zero_beta = zero_beta, method = method, kappa_max = kappa_max)
    }
    gamma_out <- gamma
  }
  
  # ── 후처리 (원본과 동일: threshold + 역표준화) ────────────
  reshape_b <- function(b) {
    B <- matrix(ifelse(abs(b) < 1e-10, 0, b), p + 1, q, byrow = TRUE)
    rbind(B[1, ], sdx.inv %*% B[-1, , drop = FALSE])
  }
  beta_new  <- reshape_b(fit0$beta)
  beta.list <- lapply(fit0$beta.list, reshape_b)
  
  n.iter    <- fit0$iterations
  crit.last <- if (length(fit0$crit.list)) tail(fit0$crit.list, 1) else NA
  converged <- is.finite(crit.last) && crit.last < eps
  
  result <- list(
    X = X0, Y = Y,
    gamma = gamma_out,
    mu    = beta_new[1, ],
    beta  = beta_new,
    beta.list   = beta.list,
    offset      = Offset,
    loglik.list = fit0$loglik.list,
    crit.list   = fit0$crit.list,
    n.iter      = n.iter,
    converged   = converged,
    status      = if (!is.null(fit0$status)) fit0$status else "ok",
    params = list(MU = mu0, Offset = Offset, maxit = maxit, eps = eps,
                  standardize = standardize, orthogonal = orthogonal,
                  lambda = lambda, method = method, engine = engine)
  )
  if (!converged)
    warning(sprintf("수렴하지 않았습니다 (iters=%d, crit=%.2e, status=%s). ",
                    n.iter, crit.last, result$status),
            "maxit을 늘리거나 method를 바꿔 보세요.")
  class(result) <- "vMFglm"
  result
}


# ── 순수 R 반복 (engine="R" 폴백): cpp와 같은 알고리즘 ──────
# vmf_fast.R 의 vmf_moments / rowKron / Cq_vec 필요
vMFglm_iteration_R <- function(Xs, Y, Offset, beta_init, eps, maxit, lambda,
                               orthogonal, gamma, zero_beta, method, kappa_max) {
  n <- nrow(Xs); p <- ncol(Xs); q <- ncol(Y); d <- (p + 1) * q
  X1 <- cbind(1, Xs)
  ridge <- c(rep(0, q), rep(lambda, p * q))
  Bmat <- function(b) matrix(b, p + 1, q, byrow = TRUE)
  
  build_pen <- function(b) {
    pen <- matrix(0, d, d)
    if (orthogonal) {
      mu <- b[1:q]; nrm <- sqrt(sum(mu^2))
      md <- if (nrm > 1e-12) mu / nrm else mu
      MMt <- tcrossprod(md)
      for (j in seq_len(p))
        pen[j*q + 1:q, j*q + 1:q] <- gamma[j] * MMt
    }
    if (!is.null(zero_beta))
      for (idx in zero_beta)
        pen[idx*q + 1:q, idx*q + 1:q] <- pen[idx*q + 1:q, idx*q + 1:q] + diag(1e12, q)
    pen
  }
  llfun <- function(b) {
    Th <- X1 %*% Bmat(b) + Offset
    r <- sqrt(rowSums(Th^2))
    sum(Cq_vec(pmax(r, 1e-300), q)) + sum(Th * Y)
  }
  
  beta <- beta_init
  pen <- build_pen(beta)
  penll <- function(b) llfun(b) - 0.5*drop(crossprod(b, pen %*% b)) -
    0.5*lambda*sum(b[-(1:q)]^2)
  ll_curr <- llfun(beta); pl_curr <- penll(beta)
  beta.list <- list(); loglik.list <- crit.list <- numeric(0)
  status <- "ok"; crit <- Inf; l1 <- 1
  
  while (crit > eps && l1 <= maxit) {
    beta.list[[l1]] <- beta
    beta_old <- beta
    Th <- X1 %*% Bmat(beta) + Offset
    mo <- vmf_moments_cpp(Th, q)
    G1 <- crossprod(X1, mo$c1 * X1)
    M  <- rowKron(X1, mo$S)
    XWX <- kronecker(G1, diag(q)) + crossprod(M, mo$c2 * M)
    R <- Y - mo$A * mo$S
    pen <- build_pen(beta)
    
    if (method == "irls") {
      d1   <- mo$c1 / (mo$c1 + lambda)
      dpar <- (mo$c1 + mo$c2) / (mo$c1 + mo$c2 + lambda)
      V <- d1 * R + ((dpar - d1) * rowSums(mo$S * R)) * mo$S
      rhs <- XWX %*% beta + as.vector(t(crossprod(X1, V)))
      cand <- drop(solve(XWX + pen + diag(1e-12, d), rhs))
      dir <- cand - beta_old; step <- 1; bt <- cand; lt <- llfun(bt); h <- 0
      while ((!is.finite(lt) || lt < ll_curr - 1e-10) && h < 20) {
        step <- step/2; bt <- beta_old + step*dir; lt <- llfun(bt); h <- h+1 }
      if (!is.finite(lt)) { status <- "step-halving failed"; break }
      beta <- bt; ll_curr <- lt
      crit <- sqrt(sum((beta_old - beta)^2))
    } else {
      sc <- as.vector(t(crossprod(X1, R))) - drop(pen %*% beta) - ridge*beta
      H <- XWX + pen; diag(H) <- diag(H) + ridge
      step_v <- drop(solve(H, sc))
      fac <- 1; bt <- beta; pl <- pl_curr
      for (h in 1:40) {
        bt <- beta_old + fac*step_v; pl <- penll(bt)
        if (is.finite(pl) && pl >= pl_curr - 1e-10) break
        fac <- fac/2 }
      if (!is.finite(pl)) { status <- "step-halving failed"; break }
      beta <- bt; pl_curr <- pl; ll_curr <- llfun(beta)
      crit <- sqrt(sum((beta_old - beta)^2)) / (1 + sqrt(sum(beta_old^2)))
      if (sqrt(sum(beta^2)) > kappa_max) {
        status <- "kappa_max reached (MLE may diverge)"
        warning("||beta|| > kappa_max: MLE가 발산할 수 있습니다.")
        loglik.list <- c(loglik.list, ll_curr); crit.list <- c(crit.list, crit)
        l1 <- l1 + 1; break
      }
    }
    loglik.list <- c(loglik.list, ll_curr)
    crit.list <- c(crit.list, crit)
    l1 <- l1 + 1
  }
  list(beta = beta, beta.list = beta.list, loglik.list = loglik.list,
       crit.list = crit.list, iterations = l1 - 1, status = status)
}









diag.matrix <- function(X, lambda){
  n=dim(X)[1]; p=dim(X)[2]
  if(n != p) stop("The matrix is not square; it has unequal dimensions.")
  
  diag(lambda, n, p)
}





Cq <- function(theta, logarithm=TRUE){
  
  q <- length(theta)
  NORM <- norm(theta, "2")
  
  
  # (NORM)^(q/2-1) / ( (2*pi)^(q/2) * besselI(NORM, q/2-1) )
  
  if(logarithm){
    
    (q/2-1) * log(NORM) - ( (q/2) * log(2*pi) + log( besselI(NORM, q/2-1, expon.scaled=TRUE) ) + NORM )
    
  } else {
    
    exp(  (q/2-1) * log(NORM) - ( (q/2) * log(2*pi) + log( besselI(NORM, q/2-1, expon.scaled=TRUE) ) + NORM )  )
    
  }
  
}




Bq <- function(theta){
  q <- length(theta)
  NORM <- norm(theta, "2")
  
  # besselI(NORM, nu=q/2, expon.scaled=FALSE) == besselI(NORM, nu=q/2, expon.scaled=TRUE)/exp(-NORM)
  # besselI(NORM, nu=q/2-1, expon.scaled=FALSE) == besselI(NORM, nu=q/2-1, expon.scaled=TRUE)/exp(-NORM)
  
  I0 <- besselI(NORM, nu=q/2-1, expon.scaled=TRUE)
  I1 <- besselI(NORM, nu=q/2, expon.scaled=TRUE)
  
  if(I1==0){
    return( NORM / ( (q/2-1) + sqrt( NORM^2 + (q/2-1)^2 ) ) )
  } else {
    return( I1 / I0 )
  }
    
  
}



Hq <- function(theta){
  
  q <- length(theta)
  NORM <- norm(theta, "2")
  
  I0 <- besselI(NORM, nu=q/2-1, expon.scaled=TRUE)
  I1 <- besselI(NORM, nu=q/2, expon.scaled=TRUE)
  
  I1_divided_by_I0 <- NORM / ( (q/2-1) + sqrt( NORM^2 + (q/2-1)^2 ) )
  
  if(I1==0){
    return( 1 - (I1_divided_by_I0)^2 - (q-1)/NORM * (I1_divided_by_I0) )
  } else {
    return( 1 - (I1/I0)^2 - (q-1)/NORM * (I1/I0) )
  }
  
  
  # numerator <- I0^2 - I1^2 - (q-1)/(NORM) * I0 * I1
  # denumerator <- I0^2
  # 
  # numerator / denumerator
}




subgrad <- function(theta){
  q <- length(theta)
  NORM <- norm(theta, "2")
  
  if(NORM == 0){
    v <- runif(q,-1,1)
    v / (norm(v, "2")*1.1)
  } else {
    theta / NORM
  }
}



#' @export b1.vMF
b1.vMF <- function(theta){
  Bq(theta) * subgrad(theta)
}



#' @export b2.vMF
b2.vMF <- function(theta){
  q <- length(theta)
  NORM <- norm(theta, "2")
  s <- subgrad(theta)
  
  # if( NORM < 1e-10 ){
  #   tcrossprod(s) * Hq(theta) + Bq(theta)
  # } else {
  tcrossprod( theta/NORM ) * Hq(theta) + Bq(theta) / NORM * (diag(1,q) - tcrossprod(theta/NORM) )
  # }
  
}






Aq_vec <- function(r, q) {
  I0 <- besselI(r, q/2 - 1, expon.scaled = TRUE)
  I1 <- besselI(r, q/2,     expon.scaled = TRUE)
  out <- ifelse(I0 > 0 & is.finite(I0) & is.finite(I1),
                I1 / I0,
                r / ((q/2 - 1) + sqrt(r^2 + (q/2 - 1)^2)))  # 근사(대규모 r 안전)
  out[r == 0] <- 0
  out
}



Cq_vec <- function(r, q) {
  (q/2 - 1) * log(r) -
    (q/2 * log(2*pi) + log(besselI(r, q/2 - 1, expon.scaled = TRUE)) + r)
}


rowKron <- function(X1, S) {
  p1 <- ncol(X1); q <- ncol(S)
  X1[, rep(seq_len(p1), each = q), drop = FALSE] *
    S[,  rep(seq_len(q), times = p1), drop = FALSE]
}


orth_complement <- function(mu) {
  mu_unit <- mu / sqrt(sum(mu^2))
  Qfull <- qr.Q(qr(matrix(mu_unit, ncol = 1)), complete = TRUE)
  Qfull[, -1, drop = FALSE]                  # q x (q-1),  t(U) %*% mu = 0
}
