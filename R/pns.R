#' @import shapes
#' @export pns
pns <- function (x, sphere.type = "seq.test", alpha = 0.1, R = 100, 
nlast.small.sphere = 1, output = TRUE){
  
  if(FALSE){
    x=t(sphere4d); output=FALSE; alpha=0.1; R=100; nlast.small.sphere=1
  }
  
  n = ncol(x)
  k = nrow(x)
  PNS = list()
  if (abs(sum(apply(x^2, 2, sum)) - n) > 1e-08) {
    stop("Error: Each column of x should be a unit vector, ||x[ , i]|| = 1.")
  }
  svd.x = svd(x, nu = nrow(x))
  uu = svd.x$u
  maxd = which(svd.x$d < 1e-15)[1]
  if (is.na(maxd) | k > n) {
    maxd = min(k, n) + 1
  }
  nullspdim = k - maxd + 1
  d = k - 1
  if (output) {
    cat("Message from pns() : dataset is on ", d, "-sphere. \n", 
        sep = "")
  }
  if (nullspdim > 0) {
    if (output) {
      cat(" .. found null space of dimension ", nullspdim, 
          ", to be trivially reduced. \n", sep = "")
    }
  }
  if (d == 2) {
    PNS$spherePNS <- t(x)
  }
  resmat = matrix(NA, d, n)
  orthaxis = list()
  orthaxis[[d - 1]] = NA
  dist = rep(NA, d - 1)
  pvalues = matrix(NA, d - 1, 2)
  ratio = rep(NA, d - 1)
  currentSphere = x
  if (nullspdim > 0) {
    for (i in 1:nullspdim) {
      oaxis = uu[, ncol(uu) - i + 1]
      r = pi/2
      pvalues[i, ] = c(NaN, NaN)
      res = acos(t(oaxis) %*% currentSphere) - r
      orthaxis[[i]] = oaxis
      dist[i] = r
      resmat[i, ] = res
      NestedSphere = rotMat(oaxis) %*% currentSphere
      currentSphere = NestedSphere[1:(k - i), ]/repmat(matrix(sqrt(1 - 
                                                                     NestedSphere[nrow(NestedSphere), ]^2), nrow = 1), 
                                                       k - i, 1)
      uu = rotMat(oaxis) %*% uu
      uu = uu[1:(k - i), ]/repmat(matrix(sqrt(1 - uu[nrow(uu), 
      ]^2), nrow = 1), k - i, 1)
      if (output) {
        cat(d - i + 1, "-sphere to ", d - i, "-sphere, by ", 
            "NULL space \n", sep = "")
      }
    }
  }
  if (sphere.type == "seq.test") {
    if (output) {
      cat(" .. sequential tests with significance level ", 
          alpha, "\n", sep = "")
    }
    isIsotropic = FALSE
    for (i in (nullspdim + 1):(d - 1)) {
      if (!isIsotropic) {
        sp = getSubSphere(x = currentSphere, geodesic = "small")
        center.s = sp$center
        r.s = sp$r
        resSMALL = acos(t(center.s) %*% currentSphere) - 
          r.s
        sp = getSubSphere(x = currentSphere, geodesic = "great")
        center.g = sp$center
        r.g = sp$r
        resGREAT = acos(t(center.g) %*% currentSphere) - 
          r.g
        pval1 = LRTpval(resGREAT, resSMALL, n)
        pvalues[i, 1] = pval1
        if (pval1 > alpha) {
          center = center.g
          r = r.g
          pvalues[i, 2] = NA
          if (output) {
            cat(d - i + 1, "-sphere to ", d - i, "-sphere, by GREAT sphere, p(LRT) = ", 
                pval1, "\n", sep = "")
          }
        }
        else {
          pval2 = vMFtest(currentSphere, R)
          pvalues[i, 2] = pval2
          if (pval2 > alpha) {
            center = center.g
            r = r.g
            if (output) {
              cat(d - i + 1, "-sphere to ", d - i, "-sphere, by GREAT sphere, p(LRT) = ", 
                  pval1, ", p(vMF) = ", pval2, "\n", sep = "")
            }
            isIsotropic = TRUE
          }
          else {
            center = center.s
            r = r.s
            if (output) {
              cat(d - i + 1, "-sphere to ", d - i, "-sphere, by SMALL sphere, p(LRT) = ", 
                  pval1, ", p(vMF) = ", pval2, "\n", sep = "")
            }
          }
        }
      }
      else if (isIsotropic) {
        sp = getSubSphere(x = currentSphere, geodesic = "great")
        center = sp$center
        r = sp$r
        if (output) {
          cat(d - i + 1, "-sphere to ", d - i, "-sphere, by GREAT sphere, restricted by testing vMF distn", 
              "\n", sep = "")
        }
        pvalues[i, 1] = NA
        pvalues[i, 2] = NA
      }
      res = acos(t(center) %*% currentSphere) - r
      orthaxis[[i]] = center
      dist[i] = r
      resmat[i, ] = res
      cur.proj = project.subsphere(x = currentSphere, 
                                   center = center, r = r)
      NestedSphere = rotMat(center) %*% currentSphere
      currentSphere = NestedSphere[1:(k - i), ]/repmat(matrix(sqrt(1 - 
                                                                     NestedSphere[nrow(NestedSphere), ]^2), nrow = 1), 
                                                       k - i, 1)
      if (nrow(currentSphere) == 3) {
        PNS$spherePNS = t(currentSphere)
      }
      if (nrow(currentSphere) == 2) {
        PNS$circlePNS = t(cur.proj)
      }
    }
  }
  else if (sphere.type == "BIC") {
    if (output) {
      cat(" .. with BIC \n")
    }
    for (i in (nullspdim + 1):(d - 1)) {
      sp = getSubSphere(x = currentSphere, geodesic = "small")
      center.s = sp$center
      r.s = sp$r
      resSMALL = acos(t(center.s) %*% currentSphere) - 
        r.s
      sp = getSubSphere(x = currentSphere, geodesic = "great")
      center.g = sp$center
      r.g = sp$r
      resGREAT = acos(t(center.g) %*% currentSphere) - 
        r.g
      BICsmall = n * log(mean(resSMALL^2)) + (d - i + 
                                                1 + 1) * log(n)
      BICgreat = n * log(mean(resGREAT^2)) + (d - i + 
                                                1) * log(n)
      if (output) {
        cat("BICsm: ", BICsmall, ", BICgr: ", BICgreat, 
            "\n", sep = "")
      }
      if (BICsmall > BICgreat) {
        center = center.g
        r = r.g
        if (output) {
          cat(d - i + 1, "-sphere to ", d - i, "-sphere, by ", 
              "GREAT sphere, BIC \n", sep = "")
        }
      }
      else {
        center = center.s
        r = r.s
        if (output) {
          cat(d - i + 1, "-sphere to ", d - i, "-sphere, by ", 
              "SMALL sphere, BIC \n", sep = "")
        }
      }
      res = acos(t(center) %*% currentSphere) - r
      orthaxis[[i]] = center
      dist[i] = r
      resmat[i, ] = res
      cur.proj = project.subsphere(x = currentSphere, 
                                   center = center, r = r)
      NestedSphere = rotMat(center) %*% currentSphere
      currentSphere = NestedSphere[1:(k - i), ]/repmat(matrix(sqrt(1 - 
                                                                     NestedSphere[nrow(NestedSphere), ]^2), nrow = 1), 
                                                       k - i, 1)
      if (nrow(currentSphere) == 3) {
        PNS$spherePNS = t(currentSphere)
      }
      if (nrow(currentSphere) == 2) {
        PNS$circlePNS = t(cur.proj)
      }
    }
  }
  else if (sphere.type == "small" | sphere.type == "great") {
    pvalues = NaN
    for (i in (nullspdim + 1):(d - 1)) {
      sp = getSubSphere(x = currentSphere, geodesic = sphere.type)
      center = sp$center
      r = sp$r
      res = acos(t(center) %*% currentSphere) - r
      orthaxis[[i]] = center
      dist[i] = r
      resmat[i, ] = res
      cur.proj = project.subsphere(x = currentSphere, 
                                   center = center, r = r)
      NestedSphere = rotMat(center) %*% currentSphere
      currentSphere = NestedSphere[1:(k - i), ]/repmat(matrix(sqrt(1 - 
                                                                     NestedSphere[nrow(NestedSphere), ]^2), nrow = 1), 
                                                       k - i, 1)
      if (nrow(currentSphere) == 3) {
        PNS$spherePNS = t(currentSphere)
      }
      if (nrow(currentSphere) == 2) {
        PNS$circlePNS = t(cur.proj)
      }
    }
  }
  else if (sphere.type == "bi.sphere") {
    if (nlast.small.sphere < 0) {
      cat("!!! Error from pns(): \n")
      cat("!!! nlast.small.sphere should be >= 0. \n")
      return(NULL)
    }
    mx = (d - 1) - nullspdim
    if (nlast.small.sphere > mx) {
      cat("!!! Error from pns(): \n")
      cat("!!! nlast.small.sphere should be <= ", mx, 
          " for this data. \n", sep = "")
      return(NULL)
    }
    pvalues = NaN
    if (nlast.small.sphere != mx) {
      for (i in (nullspdim + 1):(d - 1 - nlast.small.sphere)) {
        sp = getSubSphere(x = currentSphere, geodesic = "great")
        center = sp$center
        r = sp$r
        res = acos(t(center) %*% currentSphere) - r
        orthaxis[[i]] = center
        dist[i] = r
        resmat[i, ] = res
        cur.proj = project.subsphere(x = currentSphere, 
                                     center = center, r = r)
        NestedSphere = rotMat(center) %*% currentSphere
        currentSphere = NestedSphere[1:(k - i), ]/repmat(matrix(sqrt(1 - 
                                                                       NestedSphere[nrow(NestedSphere), ]^2), nrow = 1), 
                                                         k - i, 1)
        if (nrow(currentSphere) == 3) {
          PNS$spherePNS = t(currentSphere)
        }
        if (nrow(currentSphere) == 2) {
          PNS$circlePNS = t(cur.proj)
        }
      }
    }
    if (nlast.small.sphere != 0) {
      for (i in (d - nlast.small.sphere):(d - 1)) {
        sp = getSubSphere(x = currentSphere, geodesic = "small")
        center = sp$center
        r = sp$r
        res = acos(t(center) %*% currentSphere) - r
        orthaxis[[i]] = center
        dist[i] = r
        resmat[i, ] = res
        cur.proj = project.subsphere(x = currentSphere, 
                                     center = center, r = r)
        NestedSphere = rotMat(center) %*% currentSphere
        currentSphere = NestedSphere[1:(k - i), ]/repmat(matrix(sqrt(1 - 
                                                                       NestedSphere[nrow(NestedSphere), ]^2), nrow = 1), 
                                                         k - i, 1)
        if (nrow(currentSphere) == 3) {
          PNS$spherePNS = t(currentSphere)
        }
        if (nrow(currentSphere) == 2) {
          PNS$circlePNS = t(cur.proj)
        }
      }
    }
  }
  else {
    print("!!! Error from pns():")
    print("!!! sphere.type must be 'seq.test', 'small', 'great', 'BIC', or 'bi.sphere'")
    print("!!!   Terminating execution     ")
    return(NULL)
  }
  S1toRadian = atan2(currentSphere[2, ], currentSphere[1, 
  ])
  meantheta = geodmeanS1(S1toRadian)$geodmean
  orthaxis[[d]] = meantheta
  resmat[d, ] = mod(S1toRadian - meantheta + pi, 2 * pi) - 
    pi
  if (output) {
    par(mfrow = c(1, 1), mar = c(4, 4, 1, 1), mgp = c(2.5, 
                                                      1, 0), cex = 0.8)
    plot(currentSphere[1, ], currentSphere[2, ], xlab = "", 
         ylab = "", xlim = c(-1, 1), ylim = c(-1, 1), asp = 1)
    abline(h = 0, v = 0)
    points(cos(meantheta), sin(meantheta), pch = 1, cex = 3, 
           col = "black", lwd = 5)
    abline(a = 0, b = sin(meantheta)/cos(meantheta), lty = 3)
    l = mod(S1toRadian - meantheta + pi, 2 * pi) - pi
    points(cos(S1toRadian[which.max(l)]), sin(S1toRadian[which.max(l)]), 
           pch = 4, cex = 3, col = "blue")
    points(cos(S1toRadian[which.min(l)]), sin(S1toRadian[which.min(l)]), 
           pch = 4, cex = 3, col = "red")
    legend("topright", legend = c("Geodesic mean", "Max (+)ve from mean", 
                                  "Min (-)ve from mean"), col = c("black", "blue", 
                                                                  "red"), pch = c(1, 4, 4))
    {
      cat("\n")
      cat("length of BLUE from geodesic mean : ", max(l), 
          " (", round(max(l) * 180/pi), " degree)", "\n", 
          sep = "")
      cat("length of RED from geodesic mean : ", min(l), 
          " (", round(min(l) * 180/pi), " degree)", "\n", 
          sep = "")
      cat("\n")
    }
  }
  radii = 1
  for (i in 1:(d - 1)) {
    radii = c(radii, prod(sin(dist[1:i])))
  }
  resmat = flipud0(repmat(matrix(radii, ncol = 1), 1, n) * 
                     resmat)
  if (d > 1) {
    if (output) {
      rgl.sphgrid1()
      sphere1.f(col = "white", alpha = 0.6)
      sphrad <- 0.015
      spheres3d(-PNS$circlePNS[, 2], PNS$circlePNS[, 1], 
                PNS$circlePNS[, 3], radius = sphrad, col = 4)
      spheres3d(-PNS$spherePNS[, 2], PNS$spherePNS[, 1], 
                PNS$spherePNS[, 3], radius = sphrad, col = 2)
    }
    yy <- orthaxis[[d - 1]]
    xx <- c(-yy[2], yy[1], yy[3])
    c1 <- Enorm(c(xx[1], xx[2], xx[3]) - c(-PNS$circlePNS[1, 
                                                          2], PNS$circlePNS[1, 1], PNS$circlePNS[1, 3]))
    costheta <- 1 - c1^2/2
    angle <- (1:201)/(200) * 2 * pi
    centre <- xx * costheta
    A <- xx - centre
    B <- diag(3) - A %*% t(A)/Enorm(A)^2
    bv <- eigen(B)$vectors
    b1 <- bv[, 1]
    b2 <- bv[, 2]
    cc <- sin(acos(costheta)) * (cos(angle) %*% t(b1) + 
                                   sin(angle) %*% t(b2)) + rep(1, times = 201) %*% 
      t(centre)
    if (output) {
      lines3d(cc, col = 3, lwd = 2)
    }
    if (output) {
      lines3d(cc, col = 3, lwd = 2)
      sum <- 0
      for (i in 1:n) {
        sum = sum + (acos(cc %*% c(-PNS$circlePNS[i, 
                                                  2], PNS$circlePNS[i, 1], PNS$circlePNS[i, 
                                                                                         3])))^2
      }
      mean0angle <- which.min(sum[1:200])/200 * 2 * pi
      meanpt <- sin(acos(costheta)) * (cos(mean0angle) %*% 
                                         t(b1) + sin(mean0angle) %*% t(b2)) + t(centre)
      spheres3d(meanpt, radius = sphrad * 1.5, col = 7, 
                alpha = 0.8)
    }
  }
  PNS$scores = t(resmat)
  PNS$radii = radii
  PNS$pnscircle <- cbind(cbind(cc[, 2], -cc[, 1]), cc[, 3])
  PNS$orthaxis = orthaxis
  PNS$dist = dist
  PNS$pvalues = pvalues
  PNS$ratio = ratio
  PNS$basisu = NULL
  PNS$mean = c(PNSe2s(matrix(0, d, 1), PNS))
  meanplot <- c(-PNS$mean[2], PNS$mean[1], PNS$mean[3])
  if (output) {
    spheres3d(meanplot, radius = sphrad * 1.5, col = 7, 
              alpha = 0.8)
  }
  if (sphere.type == "seq.test") {
    PNS$sphere.type = "seq.test"
  }
  else if (sphere.type == "small") {
    PNS$sphere.type = "small"
  }
  else if (sphere.type == "great") {
    PNS$sphere.type = "great"
  }
  else if (sphere.type == "BIC") {
    PNS$sphere.type = "BIC"
  }
  else if (sphere.type == "bi.sphere") {
    PNS$sphere.type = "bi.sphere"
  }
  varPNS = apply(abs(resmat)^2, 1, sum)/n
  total = sum(varPNS)
  propPNS = varPNS/total * 100
  return(list(resmat = resmat, PNS = PNS, percent = propPNS))
}



# png.pns <- function (x, alpha = 0.1, R = 100, nlast.small.sphere = 1, output = TRUE){
#   
#   if(FALSE){
#     x=t(sphere4d); output=FALSE; alpha=0.1; R=100; nlast.small.sphere=1
#   }
#   
#   n = ncol(x)
#   k = nrow(x)
#   PNS = list()
#   
#   svd.x = svd(scale(x, center=T, scale=T), nu = nrow(x))
#   uu = svd.x$u
#   maxd = which(svd.x$d < 0.1)[1]
#   if( is.na(maxd) | k > n ) {
#     maxd = min(k, n) + 1
#   }
#   nullspdim = k - maxd + 1
#   d = k - 1
#   
#   if (d == 2) {
#     PNS$spherePNS <- t(x)
#   }
#   
#   resmat = matrix(NA, d, n)
#   orthaxis = list()
#   orthaxis[[d - 1]] = NA
#   dist = rep(NA, d - 1)
#   pvalues = matrix(NA, d - 1, 2)
#   ratio = rep(NA, d - 1)
#   currentSphere = x
#   
#   if (nullspdim > 0) {
#     for (i in 1:nullspdim) {
#       oaxis = uu[, ncol(uu) - i + 1]
#       r = pi/2
#       pvalues[i, ] = c(NaN, NaN)
#       res = acos(t(oaxis) %*% currentSphere) - r
#       orthaxis[[i]] = oaxis
#       dist[i] = r
#       resmat[i, ] = res
#       NestedSphere = rotMat(oaxis) %*% currentSphere
#       currentSphere = NestedSphere[1:(k - i), ]/repmat(matrix(sqrt(1 - NestedSphere[nrow(NestedSphere), ]^2), nrow = 1), k - i, 1)
#       uu = rotMat(oaxis) %*% uu
#       uu = uu[1:(k - i), ]/repmat(matrix(sqrt(1 - uu[nrow(uu), 
#       ]^2), nrow = 1), k - i, 1)
#       if (output) {
#         cat(d - i + 1, "-sphere to ", d - i, "-sphere, by ", 
#             "NULL space \n", sep = "")
#       }
#     }
#   }
#   
#   
#   
#   
#   {
#     
#     for (i in (nullspdim + 1):(d - 1)) {
#       
#       sp = getSubSphere(x = currentSphere, geodesic = "small")
#       center.s = sp$center
#       r.s = sp$r
#       resSMALL = acos(t(center.s) %*% currentSphere) - r.s
#       sp = getSubSphere(x = currentSphere, geodesic = "great")
#       center.g = sp$center
#       r.g = sp$r
#       resGREAT = acos(t(center.g) %*% currentSphere) - r.g
#       
#       pval1 = LRTpval(resGREAT, resSMALL, n)
#       pvalues[i, 1] = pval1
#       if (pval1 > alpha) {
#         center = center.g
#         r = r.g
#         pvalues[i, 2] = NA
#         if (output) {
#           cat(d - i + 1, "-sphere to ", d - i, "-sphere, by GREAT sphere, p(LRT) = ", 
#               pval1, "\n", sep = "")
#         }
#       } else {
#         pval2 = vMFtest(currentSphere, R)
#         pvalues[i, 2] = pval2
#         if (pval2 > alpha) {
#           center = center.g
#           r = r.g
#           if (output) {
#             cat(d - i + 1, "-sphere to ", d - i, "-sphere, by GREAT sphere, p(LRT) = ", 
#                 pval1, ", p(vMF) = ", pval2, "\n", sep = "")
#           }
#           isIsotropic = TRUE
#         } else {
#           center = center.s
#           r = r.s
#           if (output) {
#             cat(d - i + 1, "-sphere to ", d - i, "-sphere, by SMALL sphere, p(LRT) = ", 
#                 pval1, ", p(vMF) = ", pval2, "\n", sep = "")
#           }
#         }
#       }
#       
#       
#       
#       res = acos(t(center) %*% currentSphere) - r
#       orthaxis[[i]] = center
#       dist[i] = r
#       resmat[i, ] = res
#       cur.proj = project.subsphere(x = currentSphere, 
#                                    center = center, r = r)
#       NestedSphere = rotMat(center) %*% currentSphere
#       currentSphere = NestedSphere[1:(k - i), ]/repmat(matrix(sqrt(1 - 
#                                                                      NestedSphere[nrow(NestedSphere), ]^2), nrow = 1), 
#                                                        k - i, 1)
#       if (nrow(currentSphere) == 3) {
#         PNS$spherePNS = t(currentSphere)
#       }
#       if (nrow(currentSphere) == 2) {
#         PNS$circlePNS = t(cur.proj)
#       }
#     }
#     
#     
#   }
#   
#   
#   S1toRadian = atan2(currentSphere[2, ], currentSphere[1, 
#   ])
#   meantheta = geodmeanS1(S1toRadian)$geodmean
#   orthaxis[[d]] = meantheta
#   resmat[d, ] = mod(S1toRadian - meantheta + pi, 2 * pi) - pi
#   
#   if( output ){
#     par(mfrow = c(1, 1), mar = c(4, 4, 1, 1), mgp = c(2.5, 
#                                                       1, 0), cex = 0.8)
#     plot(currentSphere[1, ], currentSphere[2, ], xlab = "", 
#          ylab = "", xlim = c(-1, 1), ylim = c(-1, 1), asp = 1)
#     abline(h = 0, v = 0)
#     points(cos(meantheta), sin(meantheta), pch = 1, cex = 3, 
#            col = "black", lwd = 5)
#     abline(a = 0, b = sin(meantheta)/cos(meantheta), lty = 3)
#     l = mod(S1toRadian - meantheta + pi, 2 * pi) - pi
#     points(cos(S1toRadian[which.max(l)]), sin(S1toRadian[which.max(l)]), 
#            pch = 4, cex = 3, col = "blue")
#     points(cos(S1toRadian[which.min(l)]), sin(S1toRadian[which.min(l)]), 
#            pch = 4, cex = 3, col = "red")
#     legend("topright", legend = c("Geodesic mean", "Max (+)ve from mean", 
#                                   "Min (-)ve from mean"), col = c("black", "blue", 
#                                                                   "red"), pch = c(1, 4, 4))
#     {
#       cat("\n")
#       cat("length of BLUE from geodesic mean : ", max(l), 
#           " (", round(max(l) * 180/pi), " degree)", "\n", 
#           sep = "")
#       cat("length of RED from geodesic mean : ", min(l), 
#           " (", round(min(l) * 180/pi), " degree)", "\n", 
#           sep = "")
#       cat("\n")
#     }
#   }
#   
#   
#   radii = 1
#   for (i in 1:(d - 1)) {
#     radii = c(radii, prod(sin(dist[1:i])))
#   }
#   resmat = flipud0(repmat(matrix(radii, ncol = 1), 1, n) * 
#                      resmat)
#   
#   if (d > 1) {
#     if (output) {
#       rgl.sphgrid1()
#       sphere1.f(col = "white", alpha = 0.6)
#       sphrad <- 0.015
#       # spheres3d(-PNS$circlePNS[, 2], PNS$circlePNS[, 1], 
#       # PNS$circlePNS[, 3], radius = sphrad, col = 4)
#       spheres3d(-PNS$spherePNS[, 2], PNS$spherePNS[, 1], 
#                 PNS$spherePNS[, 3], radius = sphrad, col = 2)
#     }
#     
#     yy <- orthaxis[[d - 1]]
#     xx <- c(-yy[2], yy[1], yy[3])
#     c1 <- Enorm(c(xx[1], xx[2], xx[3]) - 
#                   c(-PNS$circlePNS[1,2], PNS$circlePNS[1,1], PNS$circlePNS[1,3]))
#     costheta <- 1 - c1^2/2
#     angle <- (1:201)/(200) * 2 * pi
#     centre <- xx * costheta
#     A <- xx - centre
#     B <- diag(3) - A %*% t(A)/Enorm(A)^2
#     bv <- eigen(B)$vectors
#     b1 <- bv[, 1]
#     b2 <- bv[, 2]
#     cc <- sin(acos(costheta)) * 
#       (cos(angle) %*% t(b1) + sin(angle) %*% t(b2)) + 
#       rep(1, times = 201) %*% t(centre)
#     
#     if (output) {
#       lines3d(cc, col = 3, lwd = 2)
#     }
#     if (output) {
#       lines3d(cc, col = 3, lwd = 2)
#       sum <- 0
#       for (i in 1:n) {
#         sum = sum + ( acos(cc %*% c(-PNS$circlePNS[i,2], PNS$circlePNS[i,1], PNS$circlePNS[i,3]) ) )^2
#       }
#       mean0angle <- which.min(sum[1:200])/200 * 2 * pi
#       meanpt <- sin(acos(costheta)) * (cos(mean0angle) %*% t(b1) + sin(mean0angle) %*% t(b2)) + t(centre)
#       spheres3d(meanpt, radius = sphrad * 1.5, col = 7, alpha = 0.8)
#     }
#   }
#   PNS$scores = t(resmat)
#   PNS$radii = radii
#   PNS$pnscircle <- cbind(cbind(cc[, 2], -cc[, 1]), cc[, 3])
#   PNS$orthaxis = orthaxis
#   PNS$dist = dist
#   PNS$pvalues = pvalues
#   PNS$ratio = ratio
#   PNS$basisu = NULL
#   PNS$mean = c(PNSe2s(matrix(0, d, 1), PNS))
#   meanplot <- c(-PNS$mean[2], PNS$mean[1], PNS$mean[3])
#   if (output) {
#     spheres3d(meanplot, radius = sphrad * 1.5, col = 7, 
#               alpha = 0.8)
#   }
#   
#   
#   varPNS = apply(abs(resmat)^2, 1, sum)/n
#   total = sum(varPNS)
#   propPNS = varPNS/total * 100
#   
#   
#   result <- list(resmat = resmat, PNS = PNS, percent = propPNS, sphere=PNS$spherePNS, circle=PNS$circlePNS)
#   class(result) <- "png.pns"
#   
#   return(result)
# }







getSubSphere <- function (x, geodesic = "small"){
  svd.x = svd(x)
  initialCenter = svd.x$u[, ncol(svd.x$u)]
  c0 = initialCenter
  TOL = 1e-10
  cnt = 0
  err = 1
  n = ncol(x)
  d = nrow(x)
  Gnow = 1e+10
  while (err > TOL) {
    c0 = c0/norm(c0, type = "2")
    rot = rotMat(c0)
    TpX = LogNPd(rot %*% x)
    fit = sphereFit(x = TpX, initialCenter = rep(0, d - 
                                                   1), geodesic = geodesic)
    newCenterTp = fit$center
    r = fit$r
    if (r > pi) {
      r = pi/2
      svd.TpX = svd(TpX)
      newCenterTp = svd.TpX$u[, ncol(svd.TpX$u)] * pi/2
    }
    newCenter = ExpNPd(newCenterTp)
    center = solve(rot, newCenter)
    Gnext = objfn(center, r, x)
    err = abs(Gnow - Gnext)
    Gnow = Gnext
    c0 = center
    cnt = cnt + 1
    if (cnt > 30) {
      break
    }
  }
  i1save = list()
  i1save$Gnow = Gnow
  i1save$center = center
  i1save$r = r
  U = princomp(t(x))$loadings[, ]
  initialCenter = U[, ncol(U)]
  c0 = initialCenter
  TOL = 1e-10
  cnt = 0
  err = 1
  n = ncol(x)
  d = nrow(x)
  Gnow = 1e+10
  while (err > TOL) {
    c0 = c0/norm(c0, type = "2")
    rot = rotMat(c0)
    TpX = LogNPd(rot %*% x)
    fit = sphereFit(x = TpX, initialCenter = rep(0, d - 
                                                   1), geodesic = geodesic)
    newCenterTp = fit$center
    r = fit$r
    if (r > pi) {
      r = pi/2
      svd.TpX = svd(TpX)
      newCenterTp = svd.TpX$u[, ncol(svd.TpX$u)] * pi/2
    }
    newCenter = ExpNPd(newCenterTp)
    center = solve(rot, newCenter)
    Gnext = objfn(center, r, x)
    err = abs(Gnow - Gnext)
    Gnow = Gnext
    c0 = center
    cnt = cnt + 1
    if (cnt > 30) {
      break
    }
  }
  if (i1save$Gnow == min(Gnow, i1save$Gnow)) {
    center = i1save$center
    r = i1save$r
  }
  if (r > pi/2) {
    center = -center
    r = pi - r
  }
  return(list(center = c(center), r = r))
}






#' @method plot plot.pns
#' @export plot.pns
plot.pns <- function(fit.pns, proj.col="grey", ...){
  
  pnscircle <- fit.pns$PNS$pnscircle
  
  fit.pns$PNS$spherePNS %>% plot.sphere(...)
  fit.pns$PNS$circlePNS %>% plot.sphere(col=proj.col, cex=0.005, add=TRUE, ...)
  
  lines3d(pnscircle[,1], pnscircle[,2], pnscircle[,3], 
          lwd=1, barblen=0.03, width=0.5, add=TRUE, n=200, col="green")
  
  df <- fit.pns$PNS$spherePNS
  angles <- calculate_view_angles(df)
  view3d(theta = angles['theta']*0.9, phi = angles['phi']*0.9, zoom = 0.75)
  
}




#' @export pns.circle.3d_to_2d
pns.circle.3d_to_2d <- function(circlePNS){
  
  t( apply( prcomp(circlePNS, scale=FALSE, center=FALSE)$x[,1:2], 1, function(x) x/norm(x,"2") ) )
  
}


