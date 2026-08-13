

if(FALSE){
  library(vMFglm)
  # detach("package:vMFglm", unload=TRUE)
  devtools::load_all()
  
  set.seed(1)
  simdata <- sim.vMFglm(n=200, mu=c(0,0,50), s_loc1=50, s_conc1=0, x_type = c("cont", "cont"))
  plot.sphere(simdata$Y, X=simdata$X[,1,drop=F], alpha=1)
  calculate_view_angles(simdata$Y) %>% {
    view3d(theta = .['theta']+2, phi = .['phi']-40, zoom = 0.75)
  }
  
  
  set.seed(2)
  simdata <- sim.vMFglm(n=200, mu=c(0,0,50), s_loc1=-50, s_conc1=0, x_type = c("binary", "cont"))
  plot.sphere(simdata$Y, X=simdata$X[,1,drop=F], alpha=1)
  view3d(theta = 0, phi = -35)
  rgl.postscript("fig01-1.pdf", fmt = "pdf")
  fit <- vMFglm(X=simdata$X[,1,drop=F], Y=simdata$Y)
  plot(fit, simdata$X[,1,drop=F], plot.conf=TRUE, plane.size=1)
  view3d(theta = 0, phi = -35)
  rgl.postscript("fig01-2.pdf", fmt = "pdf")
  
  set.seed(1)
  simdata <- sim.vMFglm(n=200, mu=c(0,0,50), s_loc1=0, s_conc1=-50, x_type = c("binary", "cont"))
  plot.sphere(simdata$Y, X=simdata$X[,1,drop=F], alpha=1)
  view3d(theta = 0, phi = -35)
  rgl.postscript("fig02-1.pdf", fmt = "pdf")
  fit <- vMFglm(X=simdata$X[,1,drop=F], Y=simdata$Y)
  plot(fit, simdata$X[,1,drop=F], plot.conf=TRUE, plane.size=1)
  view3d(theta = 0, phi = -35)
  rgl.postscript("fig02-2.pdf", fmt = "pdf")
  
  set.seed(1)
  simdata <- sim.vMFglm(n=200, mu=c(0,0,50), s_loc1=-40, s_conc1=-60, x_type = c("binary", "cont"))
  plot.sphere(simdata$Y, X=simdata$X[,1,drop=F], alpha=1)
  view3d(theta = 0, phi = -35)
  rgl.postscript("fig03-1.pdf", fmt = "pdf")
  fit <- vMFglm(X=simdata$X[,1,drop=F], Y=simdata$Y)
  plot(fit, simdata$X[,1,drop=F], plot.conf=TRUE, plane.size=1)
  view3d(theta = 0, phi = -35)
  rgl.postscript("fig03-2.pdf", fmt = "pdf")
  
  set.seed(3)
  simdata <- sim.vMFglm(n=200, mu=c(0,0,50), s_loc1=0, s_conc1=-125, x_type = c("binary", "cont"))
  plot.sphere(simdata$Y, X=simdata$X[,1,drop=F], alpha=1)
  view3d(theta = 0, phi = -35)
  rgl.postscript("fig04-1.pdf", fmt = "pdf")
  fit <- vMFglm(X=simdata$X[,1,drop=F], Y=simdata$Y)
  plot(fit, simdata$X[,1,drop=F], plot.conf=TRUE, plane.size=1)
  view3d(theta = 0, phi = -35)
  rgl.postscript("fig04-2.pdf", fmt = "pdf")
  
  set.seed(1)
  simdata <- sim.vMFglm(n=200, mu=c(0,0,50), s_loc1=-50, s_conc1=0, x_type = c("cont", "cont"))
  plot.sphere(simdata$Y, X=simdata$X[,1,drop=F], alpha=1)
  view3d(theta = 0, phi = -35)
  rgl.postscript("fig05-1.pdf", fmt = "pdf")
  fit <- vMFglm(X=simdata$X[,1,drop=F], Y=simdata$Y)
  plot(fit, simdata$X[,1,drop=F], plot.conf=TRUE, plane.size=1)
  view3d(theta = 0, phi = -35)
  rgl.postscript("fig05-2.pdf", fmt = "pdf")
  
  set.seed(2)
  simdata <- sim.vMFglm(n=200, mu=c(0,0,50), s_loc1=0, s_conc1=-20, x_type = c("cont", "cont"))
  plot.sphere(simdata$Y, X=png.utils::png.cut(simdata$X[,1,drop=F], n=5), alpha=1)
  view3d(theta = 0, phi = -35)
  rgl.postscript("fig06-1.pdf", fmt = "pdf")
  fit <- vMFglm(X=simdata$X[,1,drop=F], Y=simdata$Y)
  plot(fit, png.utils::png.cut(simdata$X[,1,drop=F], n=5), plot.conf=TRUE, plane.size=1)
  view3d(theta = 0, phi = -35)
  rgl.postscript("fig06-2.pdf", fmt = "pdf")
  
  set.seed(2)
  simdata <- sim.vMFglm(n=200, mu=c(0,0,50), s_loc1=-25, s_conc1=-25, x_type = c("cont", "cont"))
  plot.sphere(simdata$Y, X=simdata$X[,1,drop=F], alpha=1)
  view3d(theta = 0, phi = -25)
  rgl.postscript("fig07-1.pdf", fmt = "pdf")
  fit <- vMFglm(X=simdata$X[,1,drop=F], Y=simdata$Y)
  plot(fit, simdata$X[,1,drop=F], plot.conf=TRUE, plane.size=1)
  view3d(theta = 0, phi = -25)
  rgl.postscript("fig07-2.pdf", fmt = "pdf")
  
  set.seed(2)
  simdata <- sim.vMFglm(n=200, mu=c(0,0,50), s_loc1=-25, s_conc1=0, x_type = c("cont", "cont"), quad_loc1 = -25)
  plot.sphere(simdata$Y, X=simdata$X[,1,drop=F], alpha=1)
  view3d(theta = 0, phi = -25)
  rgl.postscript("fig08-1.pdf", fmt = "pdf")
  fit <- vMFglm(X=cbind.data.frame(simdata$X[,1], simdata$X[,1]^2 - mean(simdata$X[,1]^2)), Y=simdata$Y)
  plot(fit, simdata$X[,1,drop=F], plot.conf=TRUE, plane.size=1)
  view3d(theta = 0, phi = -25)
  rgl.postscript("fig08-2.pdf", fmt = "pdf")
  
  
  
}



# if(FALSE){
#   library(vMFglm)
#   
#   
#   # 1. Simulation data ----
#   set.seed(1)
#   simdata <- vMFglm::sim.vMFglm(n=50, p=1, q=3, mu=c(0,100,0), snr=50, s=100, s0=0.0, type=c("vMF", "Proj", "ExpMap")[1], dir="ortho")
#   
#   
#   # the estimated coefficient vectors (beta_j) are orthogonal (or non-orthogonal) to the mean vector (mu)
#   fit1 <- with(simdata, vMFglm(X=X, Y=Y, orthogonal=TRUE))
#   fit2 <- with(simdata, vMFglm(X=X, Y=Y, orthogonal=FALSE))
#   
#   n=nrow(simdata$Y);  q=ncol(simdata$Y)
#   
#   
#   {
#     library(dplyr)
#     X=cbind(1,simdata$X); Y=simdata$Y; MU <- simdata$mu; 
#     n=nrow(Y); q=ncol(Y); Offset <- matrix(0, n, q)
#     beta <- lm(Y ~ -1 + X, offset=tcrossprod(rep(1,n),MU)+Offset) %>% coef %>% 
#       { as.vector(t(rbind(MU, .))) }
#     
#     
#     
#   }
#   
#   
#   vMFglm_iteration(X=X, Y=Y, beta=beta, orthogonal=TRUE, Offset=matrix(0, n, q), maxit=100, eps=1e-6, lambda=5e-3)
#   
#   
#   plot.vMFglm(fit1)
#   plot.vMFglm(fit2)
#   
#   # Color mapping:
#   # - Larger values: Blue
#   # - Smaller values: Red
#   # - Intermediate values: Purple
#   
#   
#   
#   
#   # 2. Real data ----
#   
#   RData <- R.utils::loadToEnv("/Users/png/Library/CloudStorage/Dropbox/1. JSK/3. GLM using vMF/RealData/FinalData.RData")
#   
#   
#   print( RData$idx$Aligned )
#   
#   idx <- 1 # 1~5  based on length(RData$idx$Aligned)
#   # >> the 5th spoke seems to look better than the others
#   
#   
#   group <- RData$group
#   id <- RData$idx$Aligned[idx] # spoke
#   
#   X <- RData$connectionsLengths[group==1,id]
#   Y <- RData$framesBasedOnParentsUnitQuaternion[group==1,,id]
#   table( apply(Y,1,norm,"2") )
#   # >> All have the unit norm
#   
#   print( apply(Y, 2, sd) ) 
#   # >> The fourth dimension has the smallest variation >> visualize the data without it.
#   
#   Ynew <- t( apply(Y[,1:3], 1, function(x) x/norm(x, "2")) )
#   plot.sphere(Ynew)
#   
#   
#   fit1 <- vMFglm(X=X, Y=Ynew, orthogonal=TRUE) # this takes about 5.5 seconds
#   fit2 <- vMFglm(X=X, Y=Ynew, orthogonal=FALSE) # this takes about 2.6 seconds
#   
#   plot.vMFglm(fit1)
#   plot.vMFglm(fit2)
#   
#   
#   
#   
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# function(){
#   
#   library(vMFglm)
#   library(dplyr)
#   
#   
#   # set.seed(1)
#   # sim.vMFglm(n=50, p=1, q=3, p2=1, mu=c(0,0,10), snr=NULL, s=1, s0=0, t=1, t0=0, type=c("vMF", "Proj", "ExpMap"), seed.UDV=1, seed.E=NULL, orthogonal=FALSE, parallel=FALSE)
#     
#   set.seed(1)
#   simdata <- sim.vMFglm(n=500, p=1, q=3, mu=c(0,0,200), s=70, type="vMF", dir="ortho")
#   simdata$Y %>% plot.sphere()
#   
#   {
#     set.seed(1)
#     simdata <- sim.vMFglm(n=1000, p=1, q=3, mu=c(0,0,200), s=70, type="vMF", dir="ortho")
#     plot.sphere.MuBeta(df=simdata$Y, mu=simdata$mu, beta=simdata$B, X=simdata$X, scale=TRUE, scale.factor=0.5)
#     calculate_view_angles(simdata$Y) %>% {
#       view3d(theta = .['theta'], phi = .['phi']-40, zoom = 0.75)
#     }
#     
#     set.seed(1)
#     simdata <- sim.vMFglm(n=1000, p=1, q=3, mu=c(0,0,200), s=70, type="vMF", dir="random")
#     plot.sphere.MuBeta(df=simdata$Y, mu=simdata$mu, beta=simdata$B, X=simdata$X, scale=TRUE, scale.factor=0.5)
#     calculate_view_angles(simdata$Y) %>% {
#       view3d(theta = .['theta']+5, phi = .['phi']-40, zoom = 0.75)
#     }
#     
#     
#     set.seed(1)
#     simdata <- sim.vMFglm(n=1000, p=1, q=3, mu=c(0,0,200), s=70, type="vMF", dir="mu")
#     plot.sphere.MuBeta(df=simdata$Y, mu=simdata$mu, beta=simdata$B, X=png.utils::png.cut(simdata$X, n=3), scale=TRUE, scale.factor=0.5, plot.CondMean=FALSE)
#     calculate_view_angles(simdata$Y) %>% {
#       view3d(theta = .['theta']+2, phi = .['phi']-40, zoom = 0.75)
#     }
#     
#     
#   }
#   
#   
#   
#   
#   
#   
#   #
#   
#   simdata <- sim.vMFglm.TwoSampleTest2(n=50, q=3, mu=c(0,0,100), s1=100, s2=0, dir1=c("random", "mu", "ortho")[3], dir2=c("random", "mu", "ortho")[3], rho=0.5, B1=NULL, B2=NULL)
#   
#   simdata$Y %>% plot.sphere()
#   
#   
#   
#   
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Location and dispersion problems =====
# function(){
#   
#   # detach("package:vMFglm", unload=TRUE)
#   library(vMFglm)
#   
#   
#   {
#     set.seed(10)
#     simdata.loc <- sim.vMFglm.TwoSampleTest2(n=200, q=3, mu=c(0,0,100), s1=100, s2=0, dir="ortho", dir2="ortho", rho=0.0)
#     simdata.loc %>% { plot.sphere(.$Y, X=.$X1, alpha=1, cex=0.02) }
#     
#     rgl.postscript(paste0("./ToyExample.3_Loc_X1.pdf"), fmt = "pdf")
#     close3d()
#   }
#   
#   {
#     set.seed(10)
#     simdata.disp <- sim.vMFglm.TwoSampleTest2(n=200, q=3, mu=c(0,0,100), s1=200, s2=0, dir="mu", dir2="mu", rho=0.0)
#     plot.sphere(simdata.disp$Y, X=simdata.disp$X1, alpha=1, cex=0.02)
#     
#     rgl.postscript(paste0("./ToyExample.3_Disp_X1.pdf"), fmt = "pdf")
#     close3d()
#   }
#   
#   {
#     set.seed(13)
#     simdata.both <- sim.vMFglm.TwoSampleTest2(n=200, q=3, mu=c(0,0,100), s1=200, s2=0, dir="random", dir2="mu", rho=0.0)
#     plot.sphere(simdata.both$Y, X=simdata.both$X1, alpha=1, cex=0.02)
#     
#     rgl.postscript(paste0("./ToyExample.3_Both_X1.pdf"), fmt = "pdf")
#     close3d()
#   }
#   
#   
#   
#   {
#     set.seed(10)
#     simdata.loc.TwoVars<- sim.vMFglm.TwoSampleTest2(n=200, q=3, mu=c(0,0,100), s1=100, s2=100, dir="ortho", dir2="ortho", rho=0.9)
#     
#     fit1 <- simdata.loc.TwoVars %>% {
#       vMFglm.R(cbind(.$X1), .$Y, orthogonal = FALSE, gamma=0, lambda=0)
#     }
#     fit1 %>% 
#       plot(plot.conf=FALSE, plot.restricted_conf=FALSE, cex = 0.02, lwd=0.5)
#     
#     rgl.postscript(paste0("./Figures/RealData-ToyExample.3_Loc_TwoVars_X1.pdf"), fmt = "pdf")
#     close3d()
#     
#     
#     fit2 <- simdata.loc.TwoVars %>% {
#       vMFglm.R(cbind(.$X1, .$X2), .$Y, orthogonal = FALSE, gamma=0, lambda=0)
#     }
#     fit2 %>% 
#       plot(plot.conf=FALSE, plot.restricted_conf=FALSE, cex = 0.02, lwd=0.5)
#     
#     rgl.postscript(paste0("./Figures/RealData-ToyExample.3_Loc_TwoVars_X1X2.pdf"), fmt = "pdf")
#     close3d()
#     
#   }
#   
# }
# 
# 
# 
# 
# 
# 
# # An example of vMF dist =====
# function(){
#   
#   sim.vMFglm(n=500, p=1, q=3, mu=c(0,0,5^2)) %>% {
#     .$Y %>% plot.sphere()
#     angle = calculate_view_angles(.$Y)
#     {
#       view3d(theta = angle['theta'], phi = angle['phi']-30, zoom = 0.75)
#     }
#     rgl.postscript(paste0("./ToyExample-vMF_dist-kappa=5.pdf"), fmt = "pdf")
#     close3d()
#   }
#   sim.vMFglm(n=500, p=1, q=3, mu=c(0,0,10^2)) %>% {
#     .$Y %>% plot.sphere()
#     angle = calculate_view_angles(.$Y)
#     {
#       view3d(theta = angle['theta'], phi = angle['phi']-30, zoom = 0.75)
#     }
#     rgl.postscript(paste0("./ToyExample-vMF_dist-kappa=10.pdf"), fmt = "pdf")
#     close3d()
#   }
#   sim.vMFglm(n=500, p=1, q=3, mu=c(0,0,15^2)) %>% {
#     .$Y %>% plot.sphere()
#     angle = calculate_view_angles(.$Y)
#     {
#       view3d(theta = angle['theta'], phi = angle['phi']-30, zoom = 0.75)
#     }
#     rgl.postscript(paste0("./ToyExample-vMF_dist-kappa=15.pdf"), fmt = "pdf")
#     close3d()
#   }
#   
#   
#   
# }
# 







function(){
  
  set.seed(2)
  out.list <- NULL
  for(i in 1:100){
    simdata <- sim.vMFglm(n=10000, mu=c(0,0,50), s_loc1=1, s_conc1=3, x_type = c("cont", "cont"))
    fit <- vMFglm(X=simdata$X[,1:2,drop=F], Y=simdata$Y, maxit=1000, eps=1e-10)
    sm <- summary(fit)
    simdata$B
    fit$beta
    sm$test
    
    
    out.list[[i]] <- sm$conf.inside
  }
  mean(do.call("rbind", out.list)[,3] == "inside")
  
  
  
  
  
  devtools::load_all()
  
  
  
  simdata <- sim.vMFglm(n=1000, mu=c(0,0,100), s_loc1=20, s_conc1=20, x_type = c("cont", "cont"))
  microbenchmark::microbenchmark(
    # vMFglm(X=simdata$X[,1:2,drop=F], Y=simdata$Y, method="newton", maxit=1000) %>% {print(.); print(summaryfast.vMFglm(.,reduce_jacobian=FALSE)$test)},
    vMFglm(X=simdata$X[,1:2,drop=F], Y=simdata$Y, method="irls", maxit=500) %>% {summaryfast.vMFglm(.,reduce_jacobian=FALSE)$test},
    
    # vMFglm(X=simdata$X[,1:2,drop=F], Y=simdata$Y, method="newton", maxit=1000) %>% {print(.); print(summary.vMFglm(.,reduce_jacobian=FALSE)$test)},
    vMFglm(X=simdata$X[,1:2,drop=F], Y=simdata$Y, method="irls", maxit=500) %>% {summary.vMFglm(.,reduce_jacobian=FALSE)$test},
    times = 20
  )
  
  
  
}