#' @export project2vector
project2vector <- function(B, mu, out.projmat=FALSE){
  p <- length(mu)
  Proj <- tcrossprod( mu / sqrt(sum(mu^2)) )  # Projection matrix
  
  if(out.projmat){
    return( Proj )
  } else {
    return( B %*% Proj )
  }
  
}


#' @export project2orthogonal
project2orthogonal<- function(B, mu, out.projmat=FALSE) {
  p <- length(mu)
  Proj <- tcrossprod( mu / sqrt(sum(mu^2)) )  # Projection matrix
  
  if(out.projmat){
    return( diag(p) - Proj )
  } else {
    return( B %*% (diag(p) - Proj) )
  }
  
}



# Function to project a point onto a plane
#' @export project2plane
project2plane <- function(point, plane_normal, plane_point) {
  # Extract the components of the point
  x <- point[1]
  y <- point[2]
  z <- point[3]
  
  # Extract the components of the plane normal vector
  a <- plane_normal[1]
  b <- plane_normal[2]
  c <- plane_normal[3]
  
  # Extract the components of a point on the plane
  x0 <- plane_point[1]
  y0 <- plane_point[2]
  z0 <- plane_point[3]
  
  # Calculate the projection
  t <- (a * (x0 - x) + b * (y0 - y) + c * (z0 - z)) / (a^2 + b^2 + c^2)
  
  # Calculate the projected point
  x_proj <- x + t * a
  y_proj <- y + t * b
  z_proj <- z + t * c
  
  return(c(x_proj, y_proj, z_proj))
}



#' @export cross_product
cross_product <- function(a, b) {
  c(
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )
}







vMF.MuKappa <- function(Y, kappa.type=4){
  n <- nrow(Y)
  q <- ncol(Y)
  
  muhat <- colSums(Y) %>% {./norm(.,"2")} # mle of mean direction in vMF
  rbar <- colSums(Y) %>% {norm(.,"2")/n}
  
  if(is.null(kappa.type)){
    kappa0 <- 1
  } else {
    kappa0 <- list(
      (q-1) / 2*(1-rbar),
      q*rbar*( 1+q/(q+2)*rbar^2 + q^2*(q+8)/((q+2)^2*(q+4))*rbar^4 ),
      (rbar*q) / (1-rbar^2),
      (rbar*q - rbar^3) / (1-rbar^2)
    )[[kappa.type]]
  }
  
  muhat * kappa0
  
}


#


rotation.x2y <- function(x, y) {
  if(FALSE){
    x <- rnorm(5) %>% {./norm(.,"2")}
    y <- rnorm(5) %>% {./norm(.,"2")}
    rotation.x2y(x, y) %*% x
    y
  }
  # 반사시키는 벡터 v 계산
  v <- x - y
  # 벡터 v를 정규화
  v <- v / sqrt(sum(v * v))
  # 단위 행렬 생성
  I <- diag(length(x))
  # 하우스홀더 변환을 사용하여 회전 행렬 Q 계산
  Q <- I - 2 * (v %*% t(v))
  return(Q)
}


some.useful.functions <- function(x,y){
  library(rotations)
  
  rotations::project.SO3()
  
  
  #
  rotations::rot.dist()
  #
  r <- rotations::rfisher(10000, kappa = 0.01)
  #
  r <- rotations::rvmises(10000, kappa = 0.01)
  range(r)
  
  #
  rotations::genR(10)
  #
  
  
  library(rotasym)
  
}





#' @import shapes
#' @import RiemBase
# for rotMat

#' @export euclideanization
euclideanization <- function(directions1, directions2, type="tangent space") {
  if(FALSE){
    library(dplyr)
    simdata <- sim.vMFglm.TwoSampleTest(n=50, p=1, q=3)
    simdata$Y %>% head
    Y1 <- simdata$Y[simdata$Z==-0.5,]
    Y2 <- simdata$Y[simdata$Z==0.5,]
    
    euc <- euclideanization(t(Y1), t(Y2), type="tangent space")
    rbind(euc$euclideanG1, euc$euclideanG2) %>% plot(pch=18, col=as.numeric(as.factor(simdata$Z)))
    
    rbind(Y1[,1:2], Y2[,1:2]) %>% plot(pch=18, col=as.numeric(as.factor(simdata$Z)))
  }
  
  
  
  nSamplesG1<-dim(directions1)[2]
  nSamplesG2<-dim(directions2)[2]
  d<-dim(directions1)[1]
  
  pooledDirection<-cbind(directions1,directions2)
  
  # For extremely concentrated data we use Mardia mean direction
  pcaData1<-prcomp(t(pooledDirection))
  if(pcaData1$sdev[1]<1e-02 | pcaData1$sdev[2]<1e-02){
    mu_g<-convertVec2unitVec(colMeans(t(pooledDirection)))
    
    R <- rotMat(mu_g,c(rep(0,d-1),1))
    shiftedG1<-R%*%directions1
    shiftedG2<-R%*%directions2
    
    #log transfer to the tangent space
    logG1<-t(LogNPd(shiftedG1))
    logG2<-t(LogNPd(shiftedG2))
    
    result<-list(euclideanG1=logG1, euclideanG2=logG2)
    
  }else if(type=="tangent space"){
    
    # mean by Fre'chet mean
    allDirTemp<-t(pooledDirection)
    data1 <- list()
    for (j in 1:dim(allDirTemp)[1]){
      data1[[j]] <-allDirTemp[j,]
    }
    data2 <- riemfactory(data1, name="sphere")
    # Fre'chet Mean
    out1<- rbase.mean(data2)
    mu_g<-as.vector(out1$x)
    
    #rotate data to the north pole 
    R <- rotMat(mu_g,c(rep(0,d-1),1))
    shiftedG1<-R%*%directions1
    shiftedG2<-R%*%directions2
    
    #log transfer to the tangent space
    logG1<-t(LogNPd(shiftedG1))
    logG2<-t(LogNPd(shiftedG2))
    
    result<-list(euclideanG1=logG1, euclideanG2=logG2)
    
  }else if(type=="PNS"){
    
    # use pns instead of Fre'chet mean
    typeOfSphere <- kurtosisTestFunction(pooledDirection)
    pnsDirection <- pns(pooledDirection,sphere.type = typeOfSphere) #pooled directions
    res_G1<-t(pnsDirection$resmat[,1:nSamplesG1])
    res_G2<-t(pnsDirection$resmat[,(nSamplesG1+1):(nSamplesG1+nSamplesG2)])
    
    result<-list(euclideanG1=res_G1, euclideanG2=res_G2)
    
  }else{
    stop("Please choose the type of analysis i.e., 'PNS' or 'tangent space'!")
  }
  
  return(result)
}




# convert vectors to unit vectors
convertVec2unitVec <- function(vec) {
  if(norm(vec,type = "2")==0){
    stop("vector is zero!")
  }
  return(vec/norm(vec,type = "2"))
}






# In PNS
kurtosisTestFunction <- function(sphericalData, alpha=0.1) {
  ndata<-dim(sphericalData)[2]
  
  subsphereSmall<-getSubSphere(sphericalData,geodesic = "small")
  subsphereGreat<-getSubSphere(sphericalData,geodesic = "great")
  
  currentSphere<-sphericalData
  
  rSmall<-subsphereSmall$r                     # rSmall is rs in matlab
  centerSmall<-subsphereSmall$center           # NB! center is the centerSmall is pnsSmall$PNS$orthaxis[[1]]
  # and centers in matlab
  resSmall <- acos(t(centerSmall)%*%currentSphere)-rSmall  # NB!!! resSmall==(pnsSmall$resmat)[2,] i.e., residuals are second coordinates of PNS
  
  rGreat<-subsphereGreat$r                     # rGreat is rg in matlab
  centerGreat<-subsphereGreat$center           # centerGreat is centers in matlab
  resGreat <- acos(t(centerGreat)%*%currentSphere)-rGreat  # NB!!! resGreat==(pnsGreat$resmat)[2,] i.e., residuals are second coordinates of PNS
  
  # LRTpval is the likelihood ratio test from 'shapes' package
  # Chi-squared statistic for a likelihood test
  pval1 <- LRTpval(resGreat,resSmall,n = ndata)
  pval1
  
  if(pval1>alpha){
    print('great by likelihood ratio test')
    return('great')
    break
  }
  
  # # equivalently we can find pval by pns function
  # pnsTest2<-pns(sphericalData)
  # pnsTest2$PNS$pvalues
  # sum(pnsTest2$resmat[2,]==resSmall)
  
  # kurtosis test routine
  X <- LogNPd(rotMat(centerSmall) %*% currentSphere)
  
  # plot3d(t(sphericalData),type="p",expand=10, add=TRUE)
  # plot3d(t(rbind(X,rep(1,dim(X)[2]))),type="p",col = "blue",expand=10, add=TRUE)
  
  # Note that the tangential point is the center of the small circle
  d<-dim(X)[1]
  n<-dim(X)[2]
  normX2 <- colSums(X^2)
  kurtosis <- sum( normX2^2 ) / n / ( sum( normX2 ) / (d * (n-1)) )^2
  M_kurt <- d * (d+2)^2 / (d+4)
  V_kurt <- (1/n) * (128*d*(d+2)^4) / ((d+4)^3*(d+6)*(d+8))
  pval2 <- pnorm((kurtosis - M_kurt) / sqrt(V_kurt))
  
  if(pval2>alpha){
    return('great')
  }else{
    # drawCircleS2(normalVec = centerSmall,radius = rSmall)
    return('small')
  }
}



















FrechetMean <- function(Y, max_iter = 100, tol = 1e-8) {
  # manifold::frechetMean(manifold::createM("Sphere"), t(Y), ...)
  
  mu <- colMeans(Y)
  mu <- mu / sqrt(sum(mu^2))
  for (iter in seq_len(max_iter)) {
    V      <- log_map(mu, Y)
    v_mean <- colMeans(V)
    if (sqrt(sum(v_mean^2)) < tol) break
    mu <- exp_map(mu, v_mean)
  }
  return(mu)
}

# log_map <- function(mu, Y) {
#   if (is.vector(Y)) Y <- matrix(Y, nrow = 1)
#   cos_dist  <- as.numeric(Y %*% mu)
#   cos_dist  <- pmin(pmax(cos_dist, -1 + 1e-10), 1 - 1e-10)
#   theta     <- acos(cos_dist)
#   resid     <- Y - outer(cos_dist, mu)
#   resid_norm <- sqrt(rowSums(resid^2))
#   scale     <- ifelse(resid_norm < 1e-10, 1, theta / resid_norm)
#   return(resid * scale)
# }
# Logmu <- function(mu, X){
#   # library(GeodRegr)
#   # t(apply(X, 1, function(x){
#   #   GeodRegr::log_map("sphere", mu, x)
#   # }))
#   
#   t(apply(X, 1, function(x){
#     LogmapSphere(mu, x)
#   }))
#   
# }
# LogmapSphere <- function(mu, v) {
#   # 두 점 사이의 각도 계산
#   theta <- acos(sum(mu * v))
#   
#   # 0으로 나누는 것을 방지
#   if (theta < .Machine$double.eps) {
#     return(rep(0, length(mu)))
#   }
#   
#   # 구면 로그 맵 계산
#   return(theta / sin(theta) * (v - cos(theta) * mu))
#   
#   # 
#   # log_map('sphere', c(0, 0, 1), c(0, 1/sqrt(2), 1/sqrt(2)))
#   # LogmapSphere(c(0, 0, 1), c(0, 1/sqrt(2), 1/sqrt(2)))
#   
# }



exp_map <- function(mu, v) {
  norm_v <- sqrt(sum(v^2))
  if (norm_v < 1e-10) return(mu)
  y <- cos(norm_v) * mu + sin(norm_v) * (v / norm_v)
  y / sqrt(sum(y^2))
}
# Expmu <- function(mu, V) {
#   n <- dim(V)[1]
#   # lv <- mapply(function(x1, x2, x3) norm(c(x1, x2, x3), type="2"), V[, 1], V[, 2], V[, 3])
#   lv <- as.vector( apply(V, 1, norm, "2") )
#   cos_lv <- diag(cos(lv), n, n)
#   sin_lv <- diag(sin(lv) / lv, n, n)
#   im <- cos_lv %*% kronecker(matrix(mu,nrow=1), matrix(1,n,1)) + sin_lv %*% V
#   return(im)
# }


# mu의 직교 보공간 정규직교기저 (q x (q-1))
orthogonal_basis <- function(mu) {
  q <- length(mu)
  Q <- qr.Q(qr(cbind(mu, diag(q))))
  Q[, 2:q, drop = FALSE]
}