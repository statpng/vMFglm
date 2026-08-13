# ============================================================
#  plot_vMFglm.R
#
#  변경점 (vs. 이전 버전):
#  - fit.summary$cov.mu[[j]]   : 이제 Omega_conc_n (rank 1)  사용
#  - fit.summary$cov.ortho[[j]]: 이제 Omega_loc_n  (rank q-1) 사용
#  - plot.vMFglm.2d(): 새 summary 출력 구조(test$Wloc, test$Wconc)에 맞춰 수정
#  - 나머지 rgl 코드 및 2D 투영 코드는 변경 없음
# ============================================================

#' @import rgl
#' @method plot vMFglm
#' @export
plot.vMFglm <- function(fit, col.value = NULL, plot.mu = TRUE,
                           plot.conf = TRUE, plot.beta = TRUE,
                           plot.restricted_conf = TRUE,
                           plot.tangent = TRUE, use.PNS = FALSE,
                           conf.level = 0.95, scale = TRUE,
                           cex = 0.015, opacity = FALSE, add = FALSE,
                           main = NULL, arrow.type = "rotation",
                           alpha = 1.0, scale.factor = 0.8,
                           grid.lwd=0.4,
                           barblen = 0.015, lwd = 0.5, cex.tangent = 1.5,
                           lit = TRUE, n.mesh = 100,
                           plane.type = "plane",
                           plane.size = NULL,    # NULL → beta norm 기준 자동 설정
                           plane.alpha = 0.15,   # 투명도 (기본값 낮춤)
                           useNULL=FALSE,
                           ...) {
  
  # ---- extract basics ----
  X          <- fit$X;  Y <- fit$Y
  mu         <- fit$beta[1, , drop = TRUE]
  beta       <- fit$beta[-1, , drop = FALSE]
  orthogonal <- fit$params$orthogonal
  n  <- nrow(Y); p <- ncol(X); q <- ncol(Y)
  
  norm.mu <- ifelse(abs(norm(mu, "2") - 1) < 0.1,
                    norm(mu, "2"),
                    norm(mu, "2") * scale.factor)
  
  if (scale) { mu   <- mu / norm.mu }
  if (scale) { beta <- apply(beta, 1, function(x) x / norm.mu) %>% t() }
  
  # plane.size: beta arrow 최대 norm 기준으로 자동 설정
  beta_max_norm <- max(apply(beta, 1, norm, "2"))
  if (is.null(plane.size)) {
    plane.size <- beta_max_norm * 1.2   # beta 크기의 1.2배
  }
  
  MU <- tcrossprod(rep(1, nrow(beta)), mu)
  
  # ---- colour ----
  if (is.null(col.value)) {
    Xtmp <- as.matrix(X)[, ncol(X)]
    value_to_color <- function(v, lo, hi) {
      sc <- (v - lo) / (hi - lo); rgb(1 - sc, 0, sc)
    }
    col1 <- value_to_color(Xtmp, min(Xtmp), max(Xtmp))
  } else {
    value_to_color <- function(v, lo, hi) {
      sc <- (v - lo) / (hi - lo); rgb(1 - sc, 0, sc)
    }
    col1 <- value_to_color(col.value, min(col.value), max(col.value))
  }
  col2 <- "blue"; col3 <- "purple"; col4 <- "green4"
  arrow.num <- ifelse(arrow.type == "lines", 100, 3)
  
  # ---- base sphere ----
  plot.sphere(Y, col = col1, lwd=grid.lwd, opacity = opacity, add = add,
              cex = cex, main = main, alpha = alpha, lit = lit, useNULL=useNULL, ...)
  
  # ---- tangent plane ----
  if (plot.tangent) {
    if (plane.type == "plane") {
      nv <- mu; d <- -sum(nv * mu)
      a <- nv[1]; b <- nv[2]; c_coef <- nv[3]
      
      # mu 방향 기저 두 벡터 (평면 내 축)를 구해 plane.size 범위로 격자 생성
      # → mu ± plane.size 범위의 직사각형이 아니라
      #   mu에 직교하는 두 벡터로 span한 사각 패치를 그림
      e1 <- if (abs(mu[1]) < 0.9) c(1, 0, 0) else c(0, 1, 0)
      e1 <- e1 - sum(e1 * mu) * mu / sum(mu^2)
      e1 <- e1 / sqrt(sum(e1^2))
      e2 <- c(mu[2] * e1[3] - mu[3] * e1[2],
              mu[3] * e1[1] - mu[1] * e1[3],
              mu[1] * e1[2] - mu[2] * e1[1])
      e2 <- e2 / sqrt(sum(e2^2))
      
      s_seq <- seq(-plane.size, plane.size, length.out = 30)
      grid  <- expand.grid(s1 = s_seq, s2 = s_seq)
      pts   <- outer(grid$s1, e1, `*`) + outer(grid$s2, e2, `*`) +
        matrix(mu, nrow = nrow(grid), ncol = 3, byrow = TRUE)
      
      surface3d(
        x     = matrix(pts[, 1], 30, 30),
        y     = matrix(pts[, 2], 30, 30),
        z     = matrix(pts[, 3], 30, 30),
        color = "gray", alpha = plane.alpha, lit = FALSE
      )
      
    } else if (plane.type == "mesh") {
      half <- plane.size
      create_plane_mesh(a = mu, d = -sum(mu * mu),
                        alpha = plane.alpha * 0.7,
                        depth_mask = TRUE,
                        front = "lines", back = "lines",
                        n = n.mesh, col = "blue")
    }
  }
  
  # ---- beta arrows ----
  if (plot.beta) {
    for (k in seq_len(nrow(MU))) {
      rgl::arrow3d(MU[k, ], MU[k, ] + beta[k, ],
                   type = arrow.type, col = col2, s = 0.05,
                   barblen = barblen, width = lwd, add = TRUE, n = arrow.num)
      # dashed3d(rep(0, q), mu - beta[k, ], n = 100, col = col2, lwd = 2)
      # dashed3d(rep(0, q), mu + beta[k, ], n = 100, col = col2, lwd = 2)
    }
  }
  
  # ---- mu arrow ----
  if (plot.mu) {
    rgl::arrow3d(c(0, 0, 0), mu, type = arrow.type, col = col3,
                 s = 0.02, barblen = barblen, width = lwd,
                 add = TRUE, n = arrow.num)
  }
  
  # ---- constrained model: skip confidence regions ----
  if (orthogonal) return(invisible(NULL))
  
  # ---- confidence regions ----
  if (plot.conf || plot.restricted_conf) {
    # fit.summary <- summary(fit)
    fit.summary <- summary(fit, conf.level = conf.level, reduce_jacobian = TRUE)
    df_val      <- fit.summary$df
  }
  
  conf.inside <- list(wald = NULL, ProjMu = NULL, ProjOrtho = NULL)
  
  idx.list <- lapply(seq_len(p+1), function(j) (j-1)*q + seq_len(q))
  cov.list <- lapply(idx.list[-1], function(idx) fit.summary$cov.mat[idx, idx])
  
  # -- Wald CI (unchanged: uses CondSigma_beta.list) --
  if (plot.conf) {
    for (j in seq_len(p)) {
      S <- if (scale){
        cov.list[[j]] / norm.mu^2
      } else{
        cov.list[[j]]
      }
      
      for (sign in c(-1, 1)[2]) {
        png.ellipse3d(S, centre = mu + sign * beta[j, ],
                      level = conf.level,
                      t = sqrt(qchisq(conf.level, q))) %>%
          wire3d(col = col2, lit = FALSE, alpha = 0.4, lwd = grid.lwd*0.2)
      }
      
      conf.inside$wald[j] <- ifelse(
        as.numeric(t(beta[j, ]) %*% MASS::ginv(S) %*% beta[j, ]) <=
          qchisq(conf.level, q),
        "inside", "outside")
    }
  }
  
  
  
  
  if (plot.restricted_conf) {
    for (j in seq_len(p)) {
      
      ## ---- Location  (rank q-1, chi^2_{q-1}) ----
      U_perp_j <- fit.summary$U_perp[[j]]        # q x (q-1)
      Omega_dd <- fit.summary$cov.ortho[[j]]     # (q-1)x(q-1), 이미 n 으로 나뉨
      if (scale) Omega_dd <- Omega_dd / norm.mu^2
      
      S.ortho <- U_perp_j %*% Omega_dd %*% t(U_perp_j)   # q x q, rank q-1
      centre.loc <- mu + beta[j, ]
      
      png.ellipse3d(S.ortho, centre = centre.loc,
                    t = sqrt(qchisq(conf.level, q - 1))) %>%
        wire3d(col = col4, lit = FALSE, alpha = 0.4, lwd = grid.lwd*0.5)
      
      proj_red <- drop(t(U_perp_j) %*% beta[j, ])
      T_loc    <- drop(t(proj_red) %*% solve(Omega_dd) %*% proj_red)
      conf.inside$ProjOrtho[j] <- ifelse(T_loc <= qchisq(conf.level, q - 1),
                                         "inside", "outside")
      
      ## ---- Concentration  (rank 1, chi^2_1) ----
      S.mu <- fit.summary$cov.mu.full[[j]]
      if (scale) S.mu <- S.mu / norm.mu^2
      centre.conc <- mu + beta[j, ]
      
      png.ellipse3d(S.mu, centre = centre.conc,
                    t = sqrt(qchisq(conf.level, 1))) %>%
        wire3d(col = col4, lit = FALSE, alpha = 0.4, lwd = grid.lwd*0.8)
      
      s_hat <- sum((mu / sqrt(sum(mu^2))) * beta[j, ])
      var_s <- drop(t(mu / sqrt(sum(mu^2))) %*% S.mu %*% (mu / sqrt(sum(mu^2))))
      conf.inside$ProjOrtho  # (위와 동일 패턴)
      conf.inside$ProjMu[j] <- ifelse(s_hat^2 / var_s <= qchisq(conf.level, 1),
                                      "inside", "outside")
      
      
    }
  }
  
  
  
  
  # 
  # # -- Projection-based CI  [UPDATED: uses new Omega_conc_n / Omega_loc_n] --
  # if (plot.restricted_conf) {
  # 
  #   for (j in seq_len(p)) {
  # 
  #     # ---- Concentration CI  (rank 1, chi^2_1) ----
  #     S.mu <- if (scale){
  #       fit.summary$cov.mu[[j]] / norm.mu^2   # Omega_conc_n / norm.mu^2
  #     } else{
  #       fit.summary$cov.mu[[j]]
  #     }
  # 
  #     for (sign in c(-1, 1)[2]) {
  #       png.ellipse3d(S.mu, centre = mu + sign * beta[j, ],
  #                     level = conf.level,
  #                     t = sqrt(qchisq(conf.level, 1))) %>%
  #         wire3d(col = col4, lit = FALSE, alpha = 0.4, lwd = 3.5)
  #     }
  # 
  #     conf.inside$ProjMu[j] <- ifelse(
  #       as.numeric(
  #         t(fit.summary$projmat.mu %*% beta[j, ]) %*%
  #           MASS::ginv(S.mu) %*%
  #           (fit.summary$projmat.mu %*% beta[j, ])) <=
  #         qchisq(conf.level, 1),
  #       "inside", "outside")
  # 
  # 
  #     # ---- Location CI  (rank q-1, chi^2_{q-1}) ----
  #     S.ortho <- if (scale){
  #       fit.summary$cov.ortho.full[[j]] / norm.mu^2  # Omega_loc_n / norm.mu^2
  #     } else{
  #       fit.summary$cov.ortho.full[[j]]
  #     }
  # 
  #     for (sign in c(-1, 1)[2]) {
  #       png.ellipse3d(S.ortho, centre = mu + sign * beta[j, ],
  #                     level = conf.level,
  #                     t = sqrt(qchisq(conf.level, q - 1))) %>%
  #         wire3d(col = col4, lit = FALSE, alpha = 0.4, lwd = 2.5)
  #     }
  # 
  #     conf.inside$ProjOrtho[j] <- ifelse(
  #       as.numeric(
  #         t(fit.summary$projmat.ortho %*% beta[j, ]) %*%
  #           MASS::ginv(S.ortho) %*%
  #           (fit.summary$projmat.ortho %*% beta[j, ])) <=
  #         qchisq(conf.level, q - 1),
  #       "inside", "outside")
  #   }
  # }
  
  rglwidget()
  
  as.data.frame(conf.inside)
}




# ------------------------------------------------------------------
#  plot.vMFglm.2d  (updated to use Wloc / Wconc naming)
# ------------------------------------------------------------------
#' @export plot.vMFglm.2d
plot.vMFglm.2d <- function(fit, conf.level = 0.95,
                              scale.factor = 1,
                              xylim.scale  = 1.5,
                              proj.type    = c("geodesic", "euclidean"),
                              ...) {
  
  proj.type  <- match.arg(proj.type)
  X  <- fit$X;  Y <- fit$Y
  mu <- fit$mu
  B  <- fit$beta[-1, , drop = FALSE]
  p  <- ncol(X); q <- ncol(Y)
  
  sm     <- summary(fit, conf.level = conf.level)
  mu_hat <- mu / sqrt(sum(mu^2))
  
  # ---- tangent-space projection ----
  if (proj.type == "geodesic") {
    # log map at mu_hat
    project_pt <- function(y) {
      dotp <- sum(y * mu_hat)
      dotp <- pmax(pmin(dotp, 1), -1)
      v    <- y - dotp * mu_hat
      nv   <- sqrt(sum(v^2))
      if (nv < 1e-14) return(c(0, 0))
      theta <- acos(dotp)
      head(v / nv * theta, 2)
    }
  } else {
    # orthogonal projection (Euclidean)
    e1 <- sm$projmat.ortho %*% c(1, rep(0, q - 1))
    e1 <- e1 / sqrt(sum(e1^2))
    e2 <- sm$projmat.ortho %*% c(0, 1, rep(0, q - 2))
    e2 <- e2 - sum(e2 * e1) * e1
    if (sqrt(sum(e2^2)) < 1e-14) {
      e2 <- sm$projmat.ortho %*% c(0, 0, 1, rep(0, q - 3))
      e2 <- e2 - sum(e2 * e1) * e1
    }
    e2 <- e2 / sqrt(sum(e2^2))
    project_pt <- function(y) c(sum(y * e1), sum(y * e2))
  }
  
  Y_2d <- t(apply(Y, 1, project_pt))
  xlim <- range(Y_2d[, 1]) * xylim.scale
  ylim <- range(Y_2d[, 2]) * xylim.scale
  
  plot(Y_2d, pch = 20, cex = 0.6, col = "grey60",
       xlim = xlim, ylim = ylim, asp = 1,
       xlab = "Tangent axis 1", ylab = "Tangent axis 2",
       main = "Projection-based Confidence Regions", ...)
  
  cols <- c("blue", "darkorange", "green4", "red3")
  
  for (j in seq_len(p)) {
    beta_j    <- B[j, ] * scale.factor
    col_j     <- cols[(j - 1) %% length(cols) + 1]
    
    # coefficient arrow
    b2d <- project_pt(beta_j)
    arrows(0, 0, b2d[1], b2d[2], length = 0.12, lwd = 2, col = col_j)
    
    # p-values from updated test
    wl_pv  <- sm$test$Wloc["ind.1",  "pvalue"]
    wc_pv  <- sm$test$Wconc["ind.1", "pvalue"]
    leg_txt <- sprintf("beta_%d  Wloc p=%.3f  Wconc p=%.3f",
                       j, wl_pv, wc_pv)
    legend("topright", legend = leg_txt, col = col_j, lwd = 2, bty = "n",
           cex = 0.85)
  }
  abline(h = 0, v = 0, lty = 2, col = "grey70")
  invisible(NULL)
}




# ------------------------------------------------------------------
#  Helpers (unchanged from original)
# ------------------------------------------------------------------

compute_confidence_ellipse <- function(S, df, MU, beta, conf.level = 0.95) {
  fit.eig   <- eigen(S)
  chisq.val <- qchisq(conf.level, df = df)
  theta     <- seq(0, 2 * pi, length.out = 100)
  
  rotated_data2 <- list()
  for (idx in list(c(1, 2), c(1, 3), c(2, 3))) {
    eig1 <- max(fit.eig$values[idx[1]], 0)
    eig2 <- max(fit.eig$values[idx[2]], 0)
    a    <- sqrt(chisq.val * eig1)
    b    <- sqrt(chisq.val * eig2)
    ell  <- cbind(a * cos(theta), b * sin(theta))
    rot  <- t(fit.eig$vectors[, idx] %*% t(ell))
    rotated_data2 <- append(rotated_data2,
                            list(tcrossprod(rep(1, 100), MU[1, ] + beta[1, ]) + rot))
  }
  rotated_data2
}


make_pd_robust <- function(mat) {
  eig    <- eigen(mat)
  vals   <- pmax(eig$values, 1e-16)
  new_m  <- eig$vectors %*% diag(vals) %*% t(eig$vectors)
  (new_m + t(new_m)) / 2
}


png.ellipse3d <- function(x, scale = c(1, 1, 1), centre = c(0, 0, 0),
                          level = 0.95,
                          t = sqrt(qchisq(level, 3)),
                          which = 1:3, subdivide = 3,
                          smooth = TRUE, ...) {
  stopifnot(is.matrix(x))
  cov  <- x[which, which]
  cov  <- make_pd_robust(cov)
  chol_cov <- chol(cov)
  
  sphere <- subdivision3d(cube3d(...), subdivide)
  nrm    <- sqrt(sphere$vb[1, ]^2 + sphere$vb[2, ]^2 + sphere$vb[3, ]^2)
  for (i in 1:3) sphere$vb[i, ] <- sphere$vb[i, ] / nrm
  sphere$vb[4, ] <- 1
  if (smooth) sphere$normals <- sphere$vb
  
  result <- scale3d.mesh3d(transform3d(sphere, chol_cov), t, t, t)
  if (!missing(scale))  result <- scale3d.mesh3d(result, scale[1], scale[2], scale[3])
  if (!missing(centre)) result <- translate3d.mesh3d(result, centre[1], centre[2], centre[3])
  result
}
