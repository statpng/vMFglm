#' @importFrom Directional rvmf
#' @importFrom rvMF rvMF
#' @export sim.vMFglm
sim.vMFglm <- function(n       = 100,
                          q       = 3,
                          mu      = c(0, 0, 10),
                          rho     = 0,
                          s_loc1  = 0, s_conc1 = 0,
                          s_loc2  = 0, s_conc2 = 0,
                          quad_loc1 = 0, quad_conc1 = 0,
                          quad_loc2 = 0, quad_conc2 = 0,
                          x_type  = c("binary", "continuous"),
                          seed    = NULL) {
  
  s_loc = c(s_loc1, s_loc2);  s_conc = c(s_conc1, s_conc2)
  quad_loc = c(quad_loc1, quad_loc2);  quad_conc = c(quad_conc1, quad_conc2)
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  mu_norm <- sqrt(sum(mu^2));  mu_unit <- mu / mu_norm
  
  B <- NULL
  for( j in 1:2 ){
    v_raw   <- rnorm(q); v_raw <- v_raw - sum(v_raw * mu_unit) * mu_unit
    nv      <- sqrt(sum(v_raw^2))
    v_perp  <- if (nv < 1e-14) {
      e <- diag(q)[, which(abs(mu_unit) < 0.9)[1]]
      e <- e - sum(e * mu_unit) * mu_unit; e / sqrt(sum(e^2))
    } else v_raw / nv
    beta_true <- s_loc[j] * v_perp + s_conc[j] * mu_unit * sample(c(-1,1),1)
    
    B <- rbind(B, beta_true)
  }
  row.names(B) <- 1:nrow(B)
  
  beta1_true <- B[1,]
  beta2_true <- B[2,]
  
  
  
  B_quad <- NULL
  for( j in 1:2 ){
    v_raw   <- rnorm(q); v_raw <- v_raw - sum(v_raw * mu_unit) * mu_unit
    nv      <- sqrt(sum(v_raw^2))
    v_perp  <- if (nv < 1e-14) {
      e <- diag(q)[, which(abs(mu_unit) < 0.9)[1]]
      e <- e - sum(e * mu_unit) * mu_unit; e / sqrt(sum(e^2))
    } else v_raw / nv
    beta_quad_true <- quad_loc[j] * v_perp + quad_conc[j] * mu_unit * sample(c(-1,1),1)
    
    B_quad <- rbind(B_quad, beta_quad_true)
  }
  row.names(B_quad) <- 1:nrow(B_quad)
  
  beta1_quad_true <- B_quad[1,]
  beta2_quad_true <- B_quad[2,]
  
  
  
  L  <- t(chol(matrix(c(1, rho, rho, 1), 2, 2)))
  Z  <- matrix(rnorm(n * 2), n, 2) %*% t(L)
  mk <- function(z, tp) if (tp == "binary") as.numeric(z >= 0) - 0.5 else z
  X1 <- mk(Z[, 1], x_type[1]); X2 <- mk(Z[, 2], x_type[2])
  X  <- cbind(X1, X2)
  
  
  X_quad   <- apply(X, 2, function(x) x^2 - mean(x^2))   # centered quadratic term
  
  
  
  Theta <- matrix(mu, n, q, byrow = TRUE) +
    outer(X[, 1], beta1_true) + outer(X[, 2], beta2_true) +
    outer(X_quad[, 1], beta1_quad_true) + outer(X_quad[, 2], beta2_quad_true)
  
  
  Y <- t(apply(Theta, 1, function(th) {
    kap <- sqrt(sum(th^2))
    if (kap < 1e-14) {
      v <- c(1, rep(0, q - 1))
      return(v)
    }
    md <- th / kap + 1e-10
    md <- md / sqrt(sum(md^2))
    rvMF::rvMF(1, md, k = kap)
  }))
  
  list(Y = Y, X = X, B=B, B1 = beta1_true, B2 = beta2_true,
       mu = mu, mu_unit = mu_unit, v_perp = v_perp, rho = rho,
       s_loc1 = s_loc1, s_conc1 = s_conc1, 
       s_loc2 = s_loc2, s_conc2 = s_conc2, x_type=x_type)
}



















#' if(FALSE){
#'
#'   #' #' @importFrom Directional rvmf
#'   #' #' @export sim.vMFglm
#'   #' sim.vMFglm <- function(n=50, p=1, q=3, p2=1, mu=c(0,0,10), snr=NULL, s=1, s0=0, t=1, t0=0, type=c("vMF", "Proj", "ExpMap"), seed.UDV=1, seed.E=NULL, dir=c("random", "ortho", "mu")[1] ){
#'   #'
#'   #'   if(dir == "ortho"){
#'   #'     orthogonal=TRUE
#'   #'     parallel=FALSE
#'   #'   } else if(dir == "mu") {
#'   #'     orthogonal=FALSE
#'   #'     parallel=TRUE
#'   #'   } else {
#'   #'     orthogonal=FALSE
#'   #'     parallel=FALSE
#'   #'   }
#'   #'
#'   #'
#'   #'   # typeA = "sparse", typeB = "all", n = 50, p = 10, q = 10,
#'   #'   #  d = 3, rvec = NULL, nuA = 0.2, nuB = 0.5, d0 = 3, es = "1",
#'   #'   #  es.B = 1, snr = 1, simplify = TRUE, sigma = NULL, rho_X = 0.5,
#'   #'   #  rho_E = 0
#'   #'   if(FALSE){
#'   #'     n=50; p=1; q=3; r=1; mu=c(0,0,1); type="vMF"; orthogonal=TRUE; parallel=FALSE; seed.UDV=123; seed.E=123
#'   #'     t=1; t0=0
#'   #'   }
#'   #'
#'   #'   if(FALSE){
#'   #'     devtools::load_all()
#'   #'
#'   #'     simdata <- sim.vMFglm(p=1, p2=1, mu=c(0,0,50), s=10, t=50)
#'   #'     plot.sphere(simdata$Y, X=simdata$Z)
#'   #'
#'   #'     with( sim.vMFglm(p=1, p2=1, mu=c(0,0,50), s=10, t=50), plot.sphere(Y))
#'   #'     with( sim.vMFglm(p=1, mu=c(0,0,1), snr=200, s=5, type="Proj"), plot.sphere(Y))
#'   #'
#'   #'
#'   #'     simdata <- sim.vMFglm(p=1, p2=1, mu=c(0,0,10), s=10, t=50)
#'   #'     fit <- with(simdata, vMFglm(Z, Y))
#'   #'     plot(fit)
#'   #'
#'   #'
#'   #'
#'   #'   }
#'   #'
#'   #'   if(FALSE){
#'   #'     n=50; p=1; q=3; mu=c(0,0,1); snr=10; s=2; s0=0; use.ExpMap=FALSE
#'   #'   }
#'   #'
#'   #'
#'   #'
#'   #'   if(is.null(seed.E)){
#'   #'     seed.E <- sample( setdiff(1:1000, seed.UDV), 1 )
#'   #'   }
#'   #'
#'   #'   if(length(type)>1){
#'   #'     type <- "vMF"
#'   #'   }
#'   #'
#'   #'   if(length(s)==1){
#'   #'     s <- rep(s, p)
#'   #'   }
#'   #'
#'   #'
#'   #'   set.seed(seed.UDV)
#'   #'   D <- sapply(1:p, function(k) s[k] + s0)
#'   #'
#'   #'   {
#'   #'     V <- do.call("cbind", lapply(1:p, function(x){ x=rnorm(q); x/norm(x,"2") }))
#'   #'     if(orthogonal){
#'   #'
#'   #'       # original_lengths <- sqrt(colSums(V^2))
#'   #'       V <- qr.Q(qr( cbind(mu, V) ))[,-1] #%*% diag(original_lengths, ncol(V), ncol(V))
#'   #'
#'   #'     } else if(parallel){
#'   #'
#'   #'       V <- apply(V, 2, function(v) {
#'   #'         ProjMat <- (mu + runif(q,-norm(mu,"2"), norm(mu,"2"))*0.01) %>% {
#'   #'           tcrossprod(.) / c(crossprod(.))
#'   #'         }
#'   #'         ProjMat %*% v
#'   #'       })
#'   #'       V <- apply(V, 2, function(v) v/norm(v, "2"))
#'   #'
#'   #'     } else {
#'   #'       V <- apply(V, 2, function(v) v/norm(v, "2"))
#'   #'     }
#'   #'
#'   #'     V <- V %*% diag(D,p,p)
#'   #'   }
#'   #'
#'   #'   {
#'   #'     A <- do.call("cbind", lapply(1:p2, function(x) rnorm(q)))
#'   #'     if(orthogonal){
#'   #'       # original_lengths <- sqrt(colSums(V^2))
#'   #'       A <- qr.Q(qr( cbind(mu, A) ))[,-1] #%*% diag(original_lengths, ncol(V), ncol(V))
#'   #'     } else if(parallel){
#'   #'       A <- apply(A, 2, function(a) {
#'   #'         ProjMat <- (mu + runif(q,-norm(mu,"2"), norm(mu,"2"))*0.01) %>% {
#'   #'           tcrossprod(.) / c(crossprod(.))
#'   #'         }
#'   #'         ProjMat %*% a
#'   #'       })
#'   #'       A <- apply(A, 2, function(a) a/norm(a, "2"))
#'   #'     } else {
#'   #'       A <- apply(A, 2, function(a) a/norm(a, "2"))
#'   #'     }
#'   #'
#'   #'     A <- A %*% diag( sapply(1:p2, function(k) t[k] + t0), p2, p2 )
#'   #'   }
#'   #'
#'   #'
#'   #'   U <- do.call("cbind", lapply(seq_len(p), function(x) scale(rnorm(n),center=TRUE,scale=FALSE) ))
#'   #'   Z <- do.call("cbind", lapply(seq_len(p2), function(x) rbinom(n, 1, 0.5)-0.5 ))
#'   #'
#'   #'
#'   #'   if( p2 > 0 ){
#'   #'     Theta0 <- tcrossprod(U, V) + tcrossprod(Z, A)
#'   #'   } else {
#'   #'     Theta0 <- tcrossprod(U, V)
#'   #'   }
#'   #'
#'   #'   Theta <- tcrossprod(rep(1,n), mu) + Theta0
#'   #'
#'   #'
#'   #'
#'   #'
#'   #'   set.seed(seed.E)
#'   #'
#'   #'   E <- NULL
#'   #'   if( type %in% c("ExpMap", "Proj") ){
#'   #'     if(is.null(snr)){
#'   #'       snr <- 10
#'   #'     }
#'   #'
#'   #'     E0 <- matrix(rnorm(n*(q-1),0,1),n,q-1)
#'   #'     sigma <- sqrt(sum(as.numeric(Theta0)^2)/sum(as.numeric(E0)^2)/snr)
#'   #'     E0 <- E0 * sigma
#'   #'
#'   #'     E.basis <- qr.Q(qr( cbind(mu, do.call("cbind", lapply(1:(q-1), function(x) rnorm(q))) ) ))[,-1]
#'   #'     E <- E0 %*% t(E.basis)
#'   #'
#'   #'     if(FALSE){
#'   #'       rgl.sphgrid(radaxis=F, radlab=F)
#'   #'       rgl::points3d(E, col="red")
#'   #'       rgl::points3d(tcrossprod(rep(1,n),mu)+Theta0+E)
#'   #'     }
#'   #'
#'   #'   }
#'   #'
#'   #'
#'   #'
#'   #'
#'   #'   if( type == "ExpMap" ){
#'   #'
#'   #'     Y <- Expmu(mu, Theta0 + E )
#'   #'
#'   #'   } else if( type == "Proj" ){
#'   #'
#'   #'     Y <- (Theta + E) %>% apply(1, function(theta) theta/norm(theta,"2")) %>% t
#'   #'
#'   #'   } else if( type == "vMF" ){
#'   #'
#'   #'     Y <- Theta %>% apply(1, function(theta) rvmf(1, theta, k=norm(theta,"2")) ) %>% t
#'   #'
#'   #'     # Yproj <- (Theta + E) %>% apply(1, function(x) x/norm(x,"2")) %>% t
#'   #'     # Yvmf <- Theta %>% apply(1, function(x) rvmf(1, x, k=norm(x,"2")) ) %>% t
#'   #'     # {
#'   #'     #   plot.sphere(Yvmf, opacity = FALSE, add=FALSE, cex=1)
#'   #'     #   plot.sphere(Yproj, opacity = FALSE, add=TRUE, col="blue")
#'   #'     #   rgl::points3d(Theta+E, add=TRUE)
#'   #'     # }
#'   #'
#'   #'   }
#'   #'
#'   #'
#'   #'   params <- list(n=n, p=p, q=q, p2=p2, mu=mu, snr=snr, s=s, s0=s0, t=t, t0=t0, type=type, orthogonal=orthogonal, parallel=parallel)
#'   #'   result <- list(mu=mu, X=U, B=t(V), Z=Z, A=A, Y=Y, Theta=Theta, E=E, params=params)
#'   #'
#'   #'   result
#'   #' }
#'
#'
#'
#'
#'
#'   #' #' @export sim.vMFglm.TwoSampleTest
#'   #' sim.vMFglm.TwoSampleTest <- function(n=50, p=1, q=3, mu=c(0,0,10), s1=1, s2=1, dir=c("random", "mu", "ortho")){
#'   #'   library(Directional)
#'   #'
#'   #'   if(FALSE){
#'   #'     n=100; p=1; q=3; mu=c(0,0,100); snr=10; s1=100; s2=0
#'   #'     dir="ortho"
#'   #'
#'   #'     simdata <- sim.vMFglm.TwoSampleTest(n=50, p=1, q=3, mu=c(0,0,100), s1=10, s2=20, dir="ortho")
#'   #'
#'   #'     Y <- simdata$Y
#'   #'     X <- cbind(simdata$X, simdata$Z)
#'   #'     B <- rbind(simdata$mu, simdata$B, simdata$A)
#'   #'
#'   #'     plot.sphere(simdata$Y, simdata$Z)
#'   #'
#'   #'     fit <- vMFglm(simdata$X, simdata$Y)
#'   #'     B
#'   #'     fit
#'   #'
#'   #'     summary(fit)
#'   #'
#'   #'
#'   #'     plot(fit)
#'   #'     B
#'   #'
#'   #'
#'   #'     fit.perm <- png.hotelling.perm(simdata$Y, simdata$Z)
#'   #'
#'   #'
#'   #'
#'   #'   }
#'   #'
#'   #'
#'   #'
#'   #'
#'   #'   if(length(dir)>1){
#'   #'     dir = "random"
#'   #'   }
#'   #'
#'   #'   if(length(s1)==1) s1 <- rep(s1, p)
#'   #'
#'   #'
#'   #'   U <- do.call("cbind", lapply(1:p, function(j) rnorm(n)))
#'   #'   V <- do.call("cbind", lapply(1:p, function(j) rnorm(q, sd=s1[j])))
#'   #'   if(dir == "ortho"){ V <- apply(V, 2, function(v) project2orthogonal(v, mu)) }
#'   #'   if(dir == "mu"){ V <- apply(V, 2, function(v) project2vector(v, mu)) }
#'   #'
#'   #'   Z <- rbinom(n, 1, 0.5)-0.5
#'   #'   A <- rnorm(q, sd=s2)
#'   #'   if(dir == "ortho"){ A <- project2orthogonal(A, mu) %>% c() }
#'   #'   if(dir == "mu"){ A <- project2vector(A, mu) %>% c() }
#'   #'
#'   #'   Theta <- tcrossprod(rep(1,n), mu) + tcrossprod(U, V) + tcrossprod(Z, A)
#'   #'   Y <- Theta %>% apply(1, function(theta) Directional::rvmf(1, theta, k=norm(theta,"2")) ) %>% t
#'   #'
#'   #'
#'   #'   params <- list(n=n, q=q, p=p, mu=mu, s1=s1, s2=s2, dir=dir)
#'   #'   result <- list(mu=mu, X=U, B=t(V), Z=Z, A=t(A), Y=Y, Theta=Theta, params=params)
#'   #'
#'   #'   result
#'   #' }
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'   #' #' @export sim.sphere.runif
#'   #' sim.sphere.runif <- function(n){
#'   #'   aa <- seq(0, 360, length.out = n)
#'   #'   bb <- seq(0, 180, length.out = n)
#'   #'   grid <- expand.grid(theta = aa, phi = bb)
#'   #'
#'   #'   x <- sin(grid$theta * pi / 180) * cos(grid$phi * pi / 180)
#'   #'   y <- sin(grid$theta * pi / 180) * sin(grid$phi * pi / 180)
#'   #'   z <- cos(grid$theta * pi / 180)
#'   #'
#'   #'   list(x=x,y=y,z=z)
#'   #' }
#'
#'
#'
#'   #' @export sim.sphere.circle
#'   sim.sphere.circle <- function(n){
#'
#'     ndata = n
#'     theta = seq(from=-0.99,to=0.99,length.out=ndata)*pi
#'     tmpx  = cos(theta) + rnorm(ndata,sd=0.1)
#'     tmpy  = sin(theta) + rnorm(ndata,sd=0.1)
#'
#'     data <- cbind.data.frame(x=tmpx, y=tmpy, z=1)
#'     data <- t( apply(data, 1, function(x) x/norm(x,"2")) )
#'
#'     # data = riemfactory(data, name="sphere")
#'
#'     return(data)
#'   }
#'
#'
#'
#'
#'
#'
#'   #' @export sim.sphere.nonlinear
#'   sim.sphere.nonlinear <- function(n){
#'
#'     ndata = n
#'     theta = seq(from=-0.99,to=0.99,length.out=ndata)*pi
#'     tmpx  = cos(theta) + rnorm(ndata,sd=0.1)
#'     tmpy  = sin(theta) + rnorm(ndata,sd=0.1)
#'
#'     data <- cbind.data.frame(x=tmpx, y=tmpy, z=1)
#'     data <- t( apply(data, 1, function(x) x/norm(x,"2")) )
#'
#'     # data = riemfactory(data, name="sphere")
#'
#'     return(data)
#'   }
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'   #' @export orthogonalize_vectors
#'   orthogonalize_vectors <- function(mu, vectors) {
#'     # Ensure mu is a column vector
#'     mu <- as.matrix(mu)
#'     if (ncol(mu) != 1) mu <- t(mu)
#'
#'     # Combine mu and vectors into a matrix
#'     X <- cbind(mu, vectors)
#'
#'     # Perform QR decomposition
#'     qr_result <- qr(X)
#'     Q <- qr.Q(qr_result)
#'
#'     # Extract the orthogonalized vectors (excluding the first column which corresponds to mu)
#'     orthogonalized <- Q[, -1, drop = FALSE]
#'
#'     # Rescale the orthogonalized vectors to preserve their original lengths
#'     original_lengths <- sqrt(colSums(vectors^2))
#'     rescaled <- orthogonalized %*% diag(original_lengths)
#'
#'     return(rescaled)
#'   }
#'
#'
#'
#'
#'
#'
#'
#'
#'   #' @export sim.vMFglm.TwoSampleTest2
#'   sim.vMFglm.TwoSampleTest2 <- function(n=50, q=3, mu=c(0,0,10), s1=1, s2=1, dir1=c("random", "mu", "ortho"), dir2=c("random", "mu", "ortho"), rho=0.5, B1=NULL, B2=NULL){
#'     library(Directional)
#'
#'     if(FALSE){
#'
#'       devtools::load_all()
#'
#'       n=100; p=1; q=3; mu=c(0,0,100); snr=10; s1=100; s2=0;
#'       rho=0.5
#'       dir="ortho"
#'
#'       simdata <- sim.vMFglm.TwoSampleTest2(n=50, q=3, mu=c(0,0,100), s1=0, s2=100, dir="ortho", dir2="ortho", rho=0.8)
#'
#'       Y <- simdata$Y
#'       X <- cbind(simdata$X1, simdata$X2)
#'       B <- rbind(simdata$mu, simdata$B1, simdata$B2)
#'
#'       plot.sphere(simdata$Y, simdata$X2)
#'
#'       fit <- vMFglm(cbind(simdata$X1, simdata$X2), simdata$Y)
#'       plot(fit)
#'       summary(fit)
#'
#'       fit.perm <- HotellingPermTest(simdata$Y, simdata$X1)
#'       fit.perm
#'       #
#'
#'     }
#'
#'
#'
#'     if(length(dir)>1){
#'       dir = "random"
#'     }
#'
#'     Z <- mnormt::rmnorm(n, mean=rep(0,2), varcov=matrix(c(1,rho,rho,1),2,2))
#'
#'     X1 <- ifelse(Z[,1]>0,1,-1)
#'     X2 <- Z[,2]
#'
#'
#'
#'     if( is.null(B1) ){
#'       B1 <- rnorm(q)
#'       B1 <- B1/norm(B1,"2")
#'
#'       len1 <- runif(1, 0.5*s1, s1) * sample(c(-1,1), 1)
#'       if(dir1 == "ortho"){ B1 <- project2orthogonal(B1, mu) }
#'       if(dir1 == "mu"){ B1 <- project2vector(B1, mu) }
#'       B1 <- as.vector(B1/norm(B1,"2")) * len1
#'     }
#'
#'     if( is.null(B2) ){
#'       B2 <- rnorm(q)
#'       B2 <- B2/norm(B2,"2")
#'
#'       len2 <- runif(1, 0.5*s2, s2) * sample(c(-1,1), 1)
#'       if(dir2 == "ortho"){ B2 <- project2orthogonal(B2, mu) }
#'       if(dir2 == "mu"){ B2 <- project2vector(B2, mu) }
#'       B2 <- as.vector(B2/norm(B2,"2")) * len2
#'     }
#'
#'
#'
#'     Theta <- tcrossprod(rep(1,n), mu) + tcrossprod(X1, B1) + tcrossprod(X2, B2) + 1e-16
#'     Y <- Theta %>% apply(1, function(theta) Directional::rvmf(1, theta, k=norm(theta,"2")) ) %>% t
#'
#'
#'     params <- list(n=n, q=q, mu=mu, s1=s1, s2=s2, dir1=dir1, dir2=dir2)
#'     result <- list(mu=mu, X1=X1, X2=X2, B1=t(B1), B2=t(B2), Y=Y, Theta=Theta, params=params)
#'
#'     result
#'   }
#'
#'
#'
#' }
#'
