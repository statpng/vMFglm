#' #' @importFrom scales alpha
#' #' @export plot.vMFglm.2d
#' plot.vMFglm.2d <- function(fit, proj.type="central", plot.conf=TRUE, scale.factor=1.0, xylim.scale=1.0, xylim.identical=FALSE, cex=1, col=NULL, add=FALSE, alpha=0.8, lwd=2, fit.summary=NULL, conf.level=0.95, conf.color="limegreen", conf.alpha=1.0, ...){
#'   
#'   # if( !arrow.type %in% c("lines", "rotation") ) stop( 'arrow.type %in% c("lines", "rotation")' )
#'   
#'   if(FALSE){
#'     proj.type="central"
#'     col=NULL
#'     plot.conf=TRUE
#'     add=FALSE
#'     alpha=0.8
#'     cex=1
#'     lwd=2
#'     conf.level=0.95
#'     scale.factor=0.5
#'     xylim.scale=1.5
#'     xylim.identical=TRUE
#'     conf.color="limegreen"; conf.alpha=1.0
#'   }
#'   
#'   
#'   # Toy Example Figure
#'   if(FALSE){
#'     devtools::load_all()
#'     
#'     set.seed(12)
#'     simdata <- sim.vMFglm(n=50, p=2, q=3, mu=c(0,0,20), s=c(2,1), seed.UDV=5, qr.V=TRUE, qr.U=FALSE)
#'     fit <- vMFglm(simdata$X, simdata$Y, orthogonal=TRUE)
#'     summary(fit)
#'     
#'     plot(fit, scale=FALSE, conf.level = 0.5, scale.factor=0.01)
#'     
#'     plot.vMFglm.2d(fit, proj.type="central", conf.level = 0.95, scale.factor=0.1)
#'     plot.vMFglm.2d(fit, proj.type="euclidean", conf.level = 0.5, scale.factor=0.4)
#'   }
#'   
#'   
#'   
#'   # start ----
#'   
#'   
#'   {
#'     result <- list()
#'     
#'     X <- fit$X
#'     Y <- fit$Y
#'     mu <- fit$beta[1,,drop=T] * scale.factor
#'     beta <- fit$beta[-1,,drop=FALSE] * scale.factor
#'     
#'     # if(scale.factor!=1){
#'     #   norm.mu <- norm(mu, "2") * scale.factor
#'     #   mu <- fit$beta[1,,drop=T] / norm.mu
#'     #   beta <- fit$beta[-1,,drop=FALSE] / norm.mu
#'     # }
#'     
#'     n <- nrow(Y)
#'   }
#'   
#'   
#'   
#'   
#'   
#'   
#'   
#'   if(is.null(col)){
#'     col1 <- "red"
#'     col2 <- "blue"
#'     col3 <- "purple"
#'     
#'     
#'     # blue: high; red: low
#'     
#'     col1 <- value_to_color(X, min(X), max(X))
#'     
#'   } else {
#'     
#'     
#'     {
#'       a1 <- col2rgb(col)
#'       a2 <- rgb2hsv(a1)
#'       pastel_col1 <- hsv(a2[1,], a2[2,]*0.9, a2[3,])
#'       pastel_col2 <- hsv(a2[1,], a2[2,]*0.7, a2[3,])
#'       pastel_col3 <- hsv(a2[1,], a2[2,]*0.5, a2[3,])
#'       hue <- a2["h",]*360
#'       n <- 0.4
#'       pastel_col4 <- hcl(hue, 35, 85)
#'     }
#'     
#'     col1 <- pastel_col1
#'     col2 <- pastel_col2
#'     col3 <- pastel_col3
#'   }
#'   
#'   
#'   
#'   
#'   
#'   
#'   
#'   normal_vector <- mu
#'   # if( nrow(beta) == 1 ){
#'   #   normal_vector <- mu
#'   # } else {
#'   #   normal_vector <- project2plane(c(0,0,0), cross_product(beta[1,], beta[2,]), mu)
#'   # }
#'   
#'   dd <- -sum(normal_vector * mu)
#'   
#'   
#'   
#'   
#'   proj.function <- switch(proj.type,
#'                           euclidean=euclidean_project_to_plane,
#'                           central=central_project_to_plane
#'   )
#'   
#'   Y.proj <- proj.function(Y, normal_vector, dd)
#'   result$Y.proj <- Y.proj
#'   fit.pca <- prcomp( Y.proj, center=TRUE, scale=FALSE )
#'   
#'   MU <- tcrossprod(rep(1,nrow(beta)), mu)
#'   
#'   
#'   
#'   proj2pc <- function(points, fit.pca){
#'     points.centered <- points - tcrossprod( rep(1,nrow(points)), fit.pca$center )
#'     (points.centered %*% fit.pca$rotation)
#'   }
#'   
#'   
#'   
#'   # draw_ellipse2d <- function(cov_2d, center = c(0, 0), n = 100) {
#'   #   eig <- eigen(cov_2d)
#'   #   a <- sqrt(eig$values[1])
#'   #   b <- sqrt(eig$values[2])
#'   #   theta <- seq(0, 2*pi, length.out = n)
#'   #   x <- a * cos(theta)
#'   #   y <- b * sin(theta)
#'   #   points <- cbind(x, y)
#'   #   points_rotated <- t(eig$vectors %*% t(points))
#'   #   points_final <- sweep(points_rotated, 2, fit.pca$center[1:2], "+")
#'   #   return(points_final)
#'   # }
#'   
#'   
#'   
#'   from.list <- to.list <- NULL
#'   for( k in 1:nrow(beta) ){
#'     
#'     from.list[[k]] <- proj.function( matrix(MU[k,],nrow=1), normal_vector, dd) %>% proj2pc(fit.pca) %>% {.[1:2]}
#'     to.list[[k]] <- proj.function( matrix(MU[k,]+beta[k,],nrow=1), normal_vector, dd) %>% proj2pc(fit.pca) %>% {.[1:2]}
#'     
#'   }
#'   
#'   
#'   
#'   df.range <- rbind( do.call("rbind", to.list), fit.pca$x[,1:2] )
#'   
#'   if(xylim.identical){
#'     xlim <- ylim <- range(df.range) * xylim.scale
#'   } else {
#'     xlim <- range(df.range[,1]) * xylim.scale
#'     ylim <- range(df.range[,2]) * xylim.scale
#'   }
#'   
#'   
#'   if(add){
#'     points( fit.pca$x[,1:2], pch=18, col=scales::alpha(col1, alpha),
#'             cex=cex, ...)
#'   } else {
#'     plot( fit.pca$x[,1:2], pch=18, col=scales::alpha(col1, alpha),
#'           cex=cex, xlim=xlim, ylim=ylim, ...)
#'   }
#'   
#'   points( proj.function( matrix(MU[k,],nrow=1), normal_vector, dd) %>% proj2pc(fit.pca) %>% {.[1:2]} %>% matrix(.,nrow=1), pch=15, col="black", cex=cex*1.5, ...)
#'   
#'   
#'   
#'   
#'   for( k in 1:(nrow(MU)-1) ){
#'     
#'     # from <- proj.function( matrix(MU[k,],nrow=1), normal_vector, dd) %>% proj2pc(fit.pca)
#'     # to <- proj.function( matrix(MU[k,]+beta[k,],nrow=1), normal_vector, dd) %>% proj2pc(fit.pca)
#'     # 
#'     # lines(rbind(from[1:2], to[1:2]), lwd=lwd, col=scales::alpha("blue", alpha))
#'     
#'     # arrows ----
#'     
#'     # lines(rbind(from.list[[k]], to.list[[k]]), lwd=lwd, col=scales::alpha("blue", alpha))
#'     
#'     
#'     
#'     arrows(from.list[[k]][1], from.list[[k]][2], 
#'            to.list[[k]][1], to.list[[k]][2], 
#'            col=scales::alpha(conf.color, conf.alpha),
#'            length = 0.1,
#'            lwd=lwd)
#'     
#'     
#'     
#'     # library(shape)
#'     # shape::Arrows(from[1],from[2],to[1],to[2],lwd=2, arr.type="triangle", col=scales::alpha("blue", alpha))
#'     
#'     
#'     # plot.conf ----
#'     if(plot.conf){
#'       
#'       if( is.null(fit.summary) ){
#'         fit.summary <- summary(fit)
#'       }
#'       
#'       df <- fit.summary$df
#'       p <- ncol(fit$X)
#'       
#'       
#'       S <- fit.summary$cov.list[[k]] * scale.factor^2
#'       S.proj <- t(fit.pca$rotation[,1:2]) %*% S %*% fit.pca$rotation[,1:2]
#'       
#'       
#'       for(h in c(1,2,3)[1] ){
#'         fit.eig <- eigen(S.proj)
#'         chisq.val <- qchisq(conf.level - c(0.0,0.15,0.45)[h], df = df)
#'         
#'         theta <- seq(0, 2*pi, length.out = 100)
#'         a <- sqrt(chisq.val * fit.eig$values[1])
#'         b <- sqrt(chisq.val * fit.eig$values[2])
#'         
#'         ellipse_data <- cbind(a * cos(theta), b * sin(theta))
#'         
#'         rotated_data <- (t(fit.eig$vectors[,1:2] %*% t(ellipse_data)))
#'         rotated_data2 <- tcrossprod(rep(1,100), to.list[[k]]) + rotated_data[,1:2]
#'         
#'         lines(rotated_data2[,1], rotated_data2[,2], lwd=lwd, col=scales::alpha(conf.color, conf.alpha))
#'       }
#'       
#'     }
#'     
#'   }
#'   
#'   
#'   
#'   
#' }
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
#' value_to_color <- function(value, min_val, max_val) {
#'   scaled <- (value - min_val) / (max_val - min_val)
#'   r <- 1 - scaled
#'   g <- 0
#'   b <- scaled
#'   rgb(r, g, b)
#' }
#' # blue: high; red: low
#' 
#' 
#' 
#' 
#' 
#' sphere2plane <- function(points, normal_vector, d) {
#'   points <- matrix(points, ncol=3)
#'   
#'   # 평면의 방정식: ax + by + cz + d = 0
#'   a <- normal_vector[1]
#'   b <- normal_vector[2]
#'   c <- normal_vector[3]
#'   
#'   # 각 점에서 평면까지의 거리 계산
#'   distances <- (a*points[,1] + b*points[,2] + c*points[,3] + d) / 
#'     sqrt(a^2 + b^2 + c^2)
#'   
#'   # 점들을 평면에 투영
#'   projected_points <- points - outer(distances, normal_vector)
#'   
#'   projected_points
#' }
#' 
#' 
#' 
#' euclidean_project_to_plane <- function(points, normal_vector, dd=NULL) {
#'   # euclidean_project_to_plane
#'   
#'   ProjMat <- qr.Q(qr( cbind(normal_vector,1:3,c(1,3,5)) ))[,-1] %>% {
#'     . %*% solve( crossprod(.) ) %*% t(.)
#'   }
#'   
#'   projected_points <- points %*% ProjMat + tcrossprod(rep(1,nrow(points)), normal_vector)
#'   
#'   projected_points
#' }
#' 
#' 
#' central_project_to_plane <- function(points, normal_vector, dd) {
#'   # central_project_to_plane
#'   
#'   # ax + by + cz + d = 0
#'   a <- normal_vector[1]
#'   b <- normal_vector[2]
#'   c <- normal_vector[3]
#'   
#'   projected_points <- t(apply(points, 1, function(p) {
#'     t <- -dd / sum(p * normal_vector)
#'     t * p
#'   }))
#'   
#'   projected_points
#' }
#' 
