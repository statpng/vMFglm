function(){
  
  ## =============================================================
  ## Toy visualisation: θ(x) = μ + x β₁ + x² β₂  on S²
  ## =============================================================
  library(rgl)
  library(movMF)
  
  # ------------------------------------------------------------------
  # 1. Set true parameters
  # ------------------------------------------------------------------
  mu    <- c(0, 0, 50)
  beta1 <- c(10, 0, 0)
  beta2 <- c(0, 10, 5)
  
  # ------------------------------------------------------------------
  # 2. Compute θ(x) for a grid of x values
  # ------------------------------------------------------------------
  x_grid  <- seq(-2.1, 2.1, length.out = 501)
  theta_mat <- matrix(0, nrow = length(x_grid), ncol = 3)
  for (i in seq_along(x_grid)) {
    x <- x_grid[i]
    theta_mat[i, ] <- mu + x * beta1 + x^2 * beta2
  }
  
  kappa_vec <- sqrt(rowSums(theta_mat^2))
  zeta_mat  <- theta_mat / kappa_vec
  
  # ------------------------------------------------------------------
  # 3. Colour ramp: blue → purple → red
  # ------------------------------------------------------------------
  col_ramp <- colorRampPalette(c("blue", "purple", "red"))
  n_cols   <- length(x_grid)
  cols     <- col_ramp(n_cols)
  
  # ------------------------------------------------------------------
  # 4. Draw wire-frame sphere
  # ------------------------------------------------------------------
  open3d(useNULL = TRUE, windowRect = c(50, 50, 850, 850))
  bg3d("white")
  
  n_lat <- 18; n_lon <- 36
  lat <- seq(-pi/2, pi/2, length.out = n_lat + 1)
  lon <- seq(0, 2 * pi,   length.out = n_lon + 1)
  
  for (la in lat) {
    pts <- cbind(cos(la) * cos(lon), cos(la) * sin(lon), sin(la))
    lines3d(pts, col = "grey80", lwd = 0.5)
  }
  for (lo in lon) {
    pts <- cbind(cos(lat) * cos(lo), cos(lat) * sin(lo), sin(lat))
    lines3d(pts, col = "grey80", lwd = 0.5)
  }
  
  # ------------------------------------------------------------------
  # 5. Plot the mean-direction curve ζ(x) on the sphere
  # ------------------------------------------------------------------
  for (i in 1:(n_cols - 1)) {
    segments3d(zeta_mat[c(i, i + 1), ], col = cols[i], lwd = 3)
  }
  
  # Mark three key points
  idx_neg  <- 1
  idx_zero <- which.min(abs(x_grid))
  idx_pos  <- n_cols
  
  pts_key <- rbind(zeta_mat[idx_neg, ],
                   zeta_mat[idx_zero, ],
                   zeta_mat[idx_pos, ])
  spheres3d(pts_key, radius = 0.04,
            col = c("blue", "purple", "red"))
  text3d(pts_key + 0.07,
         texts = c(sprintf("x = %.1f", x_grid[idx_neg]),
                   sprintf("x = %.1f", x_grid[idx_zero]),
                   sprintf("x = %.1f", x_grid[idx_pos])),
         col = c("blue", "purple", "red"), cex = 1.2)
  
  # ------------------------------------------------------------------
  # 6. Continuous vMF samples along x (colour gradient)
  # ------------------------------------------------------------------
  set.seed(2025)
  x_sample <- seq(-2.1, 2.1, length.out = 15)  # 15 x values
  n_per    <- 30                                 # samples per x value
  
  for (k in seq_along(x_sample)) {
    x  <- x_sample[k]
    th <- mu + x * beta1 + x^2 * beta2
    ka <- sqrt(sum(th^2))
    ze <- th / ka
    
    samp <- rmovMF(n_per, theta = ka * ze)
    
    # map x to colour index
    col_idx <- max(1, min(n_cols,
                          round((x - min(x_grid)) / diff(range(x_grid)) * (n_cols - 1)) + 1))
    points3d(samp, col = adjustcolor(cols[col_idx], alpha.f = 0.4), size = 5)
  }
  
  # ------------------------------------------------------------------
  # 7. Camera angle & save
  # ------------------------------------------------------------------
  view3d(theta = 25, phi = 25, zoom = 0.75)
  rgl.snapshot("quadratic_vMFglm.png", fmt = "png")
  
  cat("\n=== Done! ===\n")
  
  
}