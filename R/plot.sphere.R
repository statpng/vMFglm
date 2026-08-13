#' @export plot.sphere
plot.sphere <- function(df, X=NULL, col="red", cex=0.01, lwd=1, opacity=FALSE, add=FALSE, main=NULL, main.cex=2, main.position="topleft", alpha=0.5, lit=TRUE, useNULL=FALSE){
  
  library(rgl)
  
  
  if(FALSE){
    df=Y; col=col1; opacity=opacity; add=add
    cex=0.01; axis.arrange=FALSE
    main <- "This is the main title"
    main.cex <- 2
  }
  
  if(FALSE){
    df = sim.vMFglm(n=100,s=1)$Y
    
    X=NULL; col="red"; cex=0.01; opacity=FALSE; add=FALSE; main=NULL; main.cex=2; main.position="topleft"
  }
  
  
  
  if(FALSE){
    df <- sim.vMFglm(n=100, s=10)$X
    col="red"; cex=0.01; opacity=TRUE; add=FALSE
  }
  
  
  
  
  
  if( ncol(df) != 3 ) stop("Check the number of columns")
  
  
  if(!is.null(X)){
    
    col1 <- "red"
    col2 <- "blue"
    col3 <- "purple"
    
    value_to_color <- function(value, min_val, max_val) {
      scaled <- (value - min_val) / (max_val - min_val)
      r <- 1 - scaled
      g <- 0
      b <- scaled
      rgb(r, g, b)
    }
    # blue: high; red: low
    
    col <- value_to_color(X, min(X), max(X))
    
  }
  
  
  
  
  # if(axis.arrange){
  #   x <- df[,1]
  #   y <- df[,3]
  #   z <- df[,2]
  # } else {
  x <- df[,1]
  y <- df[,2]
  z <- df[,3]
  # }
  
  
  if(!add){
    
    open3d(useNULL = useNULL)
    par3d(windowRect = c(20, 30, 800, 800))
    
    # getr3dDefaults()
    # rgl.par3d.names
    # par3d()$viewport
    
    
    # par3d(userMatrix = rotationMatrix(45*pi/180, 1, 0, 0),
    #       # userProjection = matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,byrow=T),
    #       # windowRect=c(0,45,800,800)
    #       windowRect=c(0,45,500,500)
    # )
    
    
    # if(!is.null(main) | !is.null(submain)){
    #   
    #   bgplot3d({
    #     plot.new()
    #     if(!is.null(main)) title(main = main, line = 3)
    #     if(!is.null(submain)) mtext(side = 1, submain, line = 4)
    #     # use here any other way you fancy to write your title
    #   }, magnify = 0.5)
    #   
    # }
    
    if(opacity){
      spheres3d(0,0,0,lit=FALSE,color="white")
      spheres3d(0,0,0,radius=1.0,lit=FALSE,color="black",front="lines")
    } else {
      plot.sphere.grid(add=TRUE, lwd=lwd)
    }
    
    
  }
  
  
  if(!is.null(main)) legend3d(main.position, legend = main, adj=0.1, cex=main.cex, bty="n", border = FALSE)
  
  # legend3d("topright", legend = paste('Type', c('A', 'B', 'C')), pch = 16, col = rainbow(3), cex=1, inset=c(0.02))
  
  
  spheres3d(x,y,z,col=col,radius=cex,alpha=alpha, lit=lit, fastTransparency = TRUE)
  
  # rgl::text3d(1.1,0,0, texts="(1,0,0)")
  # rgl::text3d(0,1.1,0, texts="(0,1,0)")
  # rgl::text3d(0,0,1.1, texts="(0,0,1)")
  # 
  # rgl::text3d(-1.1,0,0, texts="(-1,0,0)")
  # rgl::text3d(0,-1.1,0, texts="(0,-1,0)")
  # rgl::text3d(0,0,-1.1, texts="(0,0,-1)")
  
  
  
  # viewpoint ----
  # Rotate viewpoint to the direction with the data points
  {
    angles <- calculate_view_angles(df)
    view3d(theta = angles['theta'], phi = angles['phi'], zoom = 0.75)
    
    rglwidget()
  }
  
}








#' @export dashed3d
dashed3d <- function(p1, p2, n = 20, dash_frac = 0.5, ...) {
  xs <- seq(p1[1], p2[1], length.out = n)
  ys <- seq(p1[2], p2[2], length.out = n)
  zs <- seq(p1[3], p2[3], length.out = n)
  
  # Draw only some segments to mimic dashes
  for (i in seq(1, n-1, by = 2)) {
    segments3d(rbind(c(xs[i], ys[i], zs[i]),
                     c(xs[i+1], ys[i+1], zs[i+1])), ...)
  }
}



#' @export plot.sphere.MuBeta
plot.sphere.MuBeta <- function(df, mu, beta, X, col="red", cex=0.01, opacity=FALSE, add=FALSE, main=NULL, main.cex=2, main.position="topleft", alpha=0.5, lit=TRUE, scale=FALSE, scale.factor=1, plot.CondMean=TRUE, useNULL=FALSE){
  
  
  arrow.type <- "rotation"
  col2 <- "blue"
  barblen <- 0.015
  lwd=0.5
  width=lwd
  arrow.num <- ifelse(arrow.type=="lines", 100, 3)
  n=arrow.num

  
  
  if(FALSE){
    list(df=simdata$Y, mu=simdata$mu, beta=simdata$B, X=simdata$X, col="red", cex=0.01, opacity=FALSE, add=FALSE, main=NULL, main.cex=2, main.position="topleft", alpha=0.5, lit=TRUE) %>% 
      list2env(envir=.GlobalEnv)
  }
  
  library(rgl)
  
  df <- as.matrix(df)
  q <- ncol(df)
  mu <- as.vector(mu)
  beta <- as.matrix(beta)
  X <- as.matrix(X)
  
  
  if(FALSE){
    df=Y; col=col1; opacity=opacity; add=add
    cex=0.01; axis.arrange=FALSE
    main <- "This is the main title"
    main.cex <- 2
  }
  
  if(FALSE){
    df = sim.vMFglm(n=100,s=1)$Y
    
    X=NULL; col="red"; cex=0.01; opacity=FALSE; add=FALSE; main=NULL; main.cex=2; main.position="topleft"
  }
  
  
  if(FALSE){
    devtools::load_all()
    
    set.seed(1)
    
    df <- sim.vMFglm(n=100, mu = c(10,0,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    
    df <- sim.vMFglm(n=100, mu = c(0,-10,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    
    df <- sim.vMFglm(n=100, mu = c(10,10,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    
    df <- sim.vMFglm(n=100, mu = c(-10,10,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    
    df <- sim.vMFglm(n=100, mu = c(-10,-10,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    
    df <- sim.vMFglm(n=100, mu = c(-10,10,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    
    
    
    df <- sim.vMFglm(n=100, mu = c(-10,10,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    df <- sim.vMFglm(n=100, mu = c(10,-10,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    df <- sim.vMFglm(n=100, mu = c(-10,-10,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    df <- sim.vMFglm(n=100, mu = c(10,10,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    
    df <- sim.vMFglm(n=100, mu = c(-10,10,-10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    df <- sim.vMFglm(n=100, mu = c(10,-10,-10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    df <- sim.vMFglm(n=100, mu = c(-10,-10,-10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    df <- sim.vMFglm(n=100, mu = c(10,10,-10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    
    
    
    
    plot.sphere(sim.vMFglm(n=100)$Y, col="blue", opacity=TRUE, add=TRUE)
    
    set.seed(1)
    plot.sphere(sim.vMFglm(n=100, mu=c(1,0,0), s=10)$Y, opacity=T)
    set.seed(1)
    plot.sphere(sim.vMFglm(n=100, mu = c(0,1,0), s=10)$Y, col="blue", add=TRUE)
  }
  
  if(FALSE){
    df <- sim.vMFglm(n=100, s=10)$X
    col="red"; cex=0.01; opacity=TRUE; add=FALSE
  }
  
  
  
  
  
  if( ncol(df) != 3 ) stop("Check the number of columns")
  
  
  if(!is.null(X)){
    
    col1 <- "red"
    col2 <- "blue"
    col3 <- "purple"
    
    value_to_color <- function(value, min_val, max_val) {
      scaled <- (value - min_val) / (max_val - min_val)
      r <- 1 - scaled
      g <- 0
      b <- scaled
      rgb(r, g, b)
    }
    # blue: high; red: low
    
    col <- value_to_color(X, min(X), max(X))
    
  }
  
  
  
  
  # if(axis.arrange){
  #   x <- df[,1]
  #   y <- df[,3]
  #   z <- df[,2]
  # } else {
  x <- df[,1]
  y <- df[,2]
  z <- df[,3]
  # }
  
  
  if(!add){
    
    open3d(useNULL = useNULL)
    par3d(windowRect = c(20, 30, 800, 800))
    
    # getr3dDefaults()
    # rgl.par3d.names
    # par3d()$viewport
    
    
    # par3d(userMatrix = rotationMatrix(45*pi/180, 1, 0, 0),
    #       # userProjection = matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,byrow=T),
    #       # windowRect=c(0,45,800,800)
    #       windowRect=c(0,45,500,500)
    # )
    
    
    # if(!is.null(main) | !is.null(submain)){
    #   
    #   bgplot3d({
    #     plot.new()
    #     if(!is.null(main)) title(main = main, line = 3)
    #     if(!is.null(submain)) mtext(side = 1, submain, line = 4)
    #     # use here any other way you fancy to write your title
    #   }, magnify = 0.5)
    #   
    # }
    
    if(opacity){
      spheres3d(0,0,0,lit=FALSE,color="white")
      spheres3d(0,0,0,radius=1.0,lit=FALSE,color="black",front="lines")
    } else {
      plot.sphere.grid(add=TRUE, lwd=lwd)
    }
    
    
  }
  
  
  if(!is.null(main)) legend3d(main.position, legend = main, adj=0.1, cex=main.cex, bty="n", border = FALSE)
  
  # legend3d("topright", legend = paste('Type', c('A', 'B', 'C')), pch = 16, col = rainbow(3), cex=1, inset=c(0.02))
  
  
  spheres3d(x,y,z,col=col,radius=cex,alpha=alpha, lit=lit, fastTransparency = TRUE)
  
  # rgl::text3d(1.1,0,0, texts="(1,0,0)")
  # rgl::text3d(0,1.1,0, texts="(0,1,0)")
  # rgl::text3d(0,0,1.1, texts="(0,0,1)")
  # 
  # rgl::text3d(-1.1,0,0, texts="(-1,0,0)")
  # rgl::text3d(0,-1.1,0, texts="(0,-1,0)")
  # rgl::text3d(0,0,-1.1, texts="(0,0,-1)")
  
  
  norm.mu <- ifelse( abs(norm(mu, "2")-1) < 0.1, norm(mu, "2"), norm(mu, "2") * scale.factor )
  
  if(scale) mu <- mu / norm.mu
  if(scale) beta <- apply(beta, 1, function(x) x / norm.mu ) %>% t()
  
  
  
  if( !is.null(mu) ){
    
    switch(3,
           `1`=rgl::arrow3d(c(0,0,0), mu, type=c("lines", "rotation")[1], col=col3, s=0.1, barblen=0.02, width=lwd, add=TRUE),
           `2`=rgl::arrow3d(c(0,0,0), mu, type=c("lines", "rotation")[2], col=col3, s=0.1, barblen=0.02, width=lwd, add=TRUE),
           `3`=rgl::arrow3d(c(0,0,0), mu, type=arrow.type, col=col3, s=0.02, barblen=barblen, width=lwd, add=TRUE, n=arrow.num)
    )
    
  }
  
  

  
  if( !is.null(beta) & !is.null(mu) ){
    
    for( k in 1:nrow(beta) ){
      # rgl::arrow3d(MU[k,]-beta[k,], MU[k,]+beta[k,], type=c("lines", "rotation")[1], col=col2, s=0.05, barblen=0.03, width=0.5, add=TRUE, n=100)
      
      rgl::arrow3d(mu, mu+beta[k,], type=arrow.type, col=col2, s=0.05, 
                   # barblen=0.04*sqrt(norm(beta[k,], "2")),
                   barblen=barblen,
                   width=lwd, add=TRUE, n=arrow.num)
      
      rgl::arrow3d(mu, mu-beta[k,], type=arrow.type, col=col2, s=0.05, 
                   # barblen=0.04*sqrt(norm(beta[k,], "2")),
                   barblen=barblen,
                   width=lwd, add=TRUE, n=arrow.num)
      
      dashed3d(rep(0,q), mu-beta[k,], n=100, col=col2, lwd=2)
      dashed3d(rep(0,q), mu+beta[k,], n=100, col=col2, lwd=2)
      
      
      df.MuBeta <- rbind(
        (mu+beta[k,])/norm(mu+beta[k,],"2"),
        (mu-beta[k,])/norm(mu-beta[k,],"2")
      )
      
      
      if(plot.CondMean){
        spheres3d(df.MuBeta[,1], df.MuBeta[,2], df.MuBeta[,3], col="orange",radius=cex*5, alpha=0.9, lit=FALSE, fastTransparency = TRUE)
      }
      
      # 
      # rgl::points3d( (mu+beta[k,])/norm(mu+beta[k,],"2"), col="red", size = 20 )
      # rgl::points3d( (mu-beta[k,])/norm(mu-beta[k,],"2"), col="red", size = 20 )
      
      
      
      # rgl::arrow3d(MU[k,]-beta[k,], MU[k,]+beta[k,], type=arrow.type, col=col2, s=0.05, 
      #              # barblen=0.04*sqrt(norm(beta[k,], "2")),
      #              barblen=barblen,
      #              width=lwd, add=TRUE, n=arrow.num)
      
    }
    
  }
  
  
  
  
  # viewpoint ----
  # Rotate viewpoint to the direction with the data points
  {
    angles <- calculate_view_angles(df)
    view3d(theta = angles['theta'], phi = angles['phi'], zoom = 0.75)
    
    # {
    #   current_matrix <- par3d("userMatrix")
    #   additional_rotation <- get_pca_view_matrix(df)
    #   new_matrix <- additional_rotation %*% current_matrix
    #   par3d(userMatrix = new_matrix)
    # }
    
  }
  
}






#' @export calculate_view_angles
calculate_view_angles <- function(df) {
  
  rotation_matrix_y <- function(angle) {
    angle_rad <- angle * pi / 180
    matrix(c(cos(angle_rad), 0, sin(angle_rad),
             0, 1, 0,
             -sin(angle_rad), 0, cos(angle_rad)), 
           nrow = 3, byrow = TRUE)
  }
  
  rotation_matrix_x <- function(angle) {
    angle_rad <- angle * pi / 180
    matrix(c(1, 0, 0,
             0, cos(angle_rad), sin(angle_rad),
             0, -sin(angle_rad), cos(angle_rad)), 
           nrow = 3, byrow = TRUE)
  }
  
  # 정규화된 평균 벡터 계산
  xbar <- colMeans(df)
  xbar <- xbar / sqrt(sum(xbar^2))
  
  # theta 계산
  v_xz <- c(xbar[1], 0, xbar[3])
  v_xz <- v_xz / sqrt(sum(v_xz^2))
  theta1 <- acos(sum(v_xz * c(0,0,1))) * 180 / pi
  theta2 <- acos(-sum(v_xz * c(0,0,1))) * 180 / pi - 180
  
  # 최적의 theta 선택
  v1 <- rotation_matrix_y(theta1) %*% c(0,0,1)
  v2 <- rotation_matrix_y(theta2) %*% c(0,0,1)
  theta <- ifelse(sum(xbar * v1) > sum(xbar * v2), theta1, theta2)
  
  # phi 계산을 위한 중간 단계
  R_theta <- rotation_matrix_y(theta)
  TT <- R_theta %*% c(0,0,1)
  
  phi1 <- acos(sum(xbar * (R_theta %*% c(0,0,1)))) * 180 / pi
  phi2 <- acos(-sum(xbar * (R_theta %*% c(0,0,1)))) * 180 / pi - 180
  
  # 최적의 phi 선택
  v1 <- rotation_matrix_x(phi1) %*% TT
  v2 <- rotation_matrix_x(phi2) %*% TT
  phi <- ifelse(sum(xbar * v1) > sum(xbar * v2), phi1, phi2)
  
  if(abs(theta) > 90){
    phi <- -phi
  }
  
  return(c(theta = theta, phi = phi))
  
  # angles <- calculate_view_angles(df)
  # view3d(theta = angles['theta'], phi = angles['phi'], zoom = 0.75)
}








#' @export plot.sphere.grid
plot.sphere.grid <- function (radius = 1, lwd=1, col.long = "red", col.lat = "blue", deggap = 15, longtype = "H", add = FALSE, radaxis = TRUE, radlab = "Radius", useNULL=FALSE){
  if(FALSE){
    radius = 1; col.long = "red"; col.lat = "blue"; deggap = 15;
    longtype = "H"; add = FALSE; radaxis = TRUE; radlab = "Radius"
  }
  
  
  # if(FALSE){
  #   library(sphereplot)
  #   sphereplot::pointsphere()
  #   
  #   rgl.sphgrid()
  #   rgl.sphpoints(pointsphere(100,c(0,90),c(0,45),c(0.25,0.8)),deg=T)
  #   
  #   rgl.sphgrid(radaxis=F, radlab=F)
  #   rgl.sphpoints(40,50,0.5,deg=TRUE,col='red',cex=2)
  # }
  
  
  
  if (add == F) {
    open3d(useNULL = useNULL)
    
    par3d(userMatrix = rotationMatrix(-15*pi/180, 1, 0, 0),
          # userProjection = matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,byrow=T),
          # windowRect=c(0,45,800,800)
          windowRect=c(0,45,800,800)
    )
    
  }
  for (lat in seq(-90, 90, by = deggap)) {
    if (lat == 0) {
      col.grid = "grey50"
    }
    else {
      col.grid = "grey"
    }
    plot3d(png.sph2coord(long = seq(0, 360, len = 100), lat = lat,
                         radius = radius, deg = T), col = col.grid, add = T, lwd=lwd,
           type = "l")
  }
  for (long in seq(0, 360 - deggap, by = deggap)) {
    if (long == 0) {
      col.grid = "grey50"
    }
    else {
      col.grid = "grey"
    }
    plot3d(png.sph2coord(long = long, lat = seq(-90, 90, len = 100),
                         radius = radius, deg = T), col = col.grid, add = T, lwd=lwd,
           type = "l")
  }
  if (longtype == "H") {
    scale = 15
  }
  if (longtype == "D") {
    scale = 1
  }
  # rgl.sphtext(long = 0, lat = seq(-90, 90, by = deggap), radius = radius,
  #             text = seq(-90, 90, by = deggap), deg = TRUE, col = col.lat)
  # rgl.sphtext(long = seq(0, 360 - deggap, by = deggap), lat = 0,
  #             radius = radius, text = seq(0, 360 - deggap, by = deggap)/scale,
  #             deg = TRUE, col = col.long)
  # if (radaxis) {
  #   radpretty = pretty(c(0, radius))
  #   radpretty = radpretty[radpretty <= radius]
  #   lines3d(c(0, 0), c(0, max(radpretty)), c(0, 0), col = "grey50")
  #   for (i in 1:length(radpretty)) {
  #     lines3d(c(0, 0), c(radpretty[i], radpretty[i]),
  #             c(0, 0, radius/50), col = "grey50")
  #     text3d(0, radpretty[i], radius/15, radpretty[i],
  #            col = "darkgreen")
  #   }
  #   text3d(0, radius/2, -radius/25, radlab)
  # }
  
}



png.coord2sph <- function(x,y,z){
  R=1
  lat = asin(z / R)
  long = atan2(y, x)
  cbind.data.frame(long=long, lat=lat)
}





png.sph2coord <- function (long, lat, radius = 1, deg = TRUE){
  # if(FALSE){
  #   x = R * cos(lat) * cos(long)
  #   y = R * cos(lat) * sin(long)
  #   z = R * sin(lat)
  #   cbind.data.frame(x=x,y=y,z=z)
  # }
  
  
  if (is.matrix(long) || is.data.frame(long)) {
    if (ncol(long) == 1) {
      long = long[, 1]
    }
    else if (ncol(long) == 2) {
      lat = long[, 2]
      long = long[, 1]
    }
    else if (ncol(long) == 3) {
      radius = long[, 3]
      lat = long[, 2]
      long = long[, 1]
    }
  }
  if (missing(long) | missing(lat)) {
    stop("Missing full spherical 3D input data.")
  }
  if (deg) {
    long = long * pi/180
    lat = lat * pi/180
  }
  return = cbind(x = radius*cos(lat)*cos(long),
                 y = radius*cos(lat)*sin(long),
                 z = radius*sin(lat))
}

























create_plane_mesh <- function(a, d, alpha = 0.05, depth_mask = FALSE, front = "lines", back = "lines", n = 20, col = "blue") {
  # a는 normal vector (길이 3의 벡터)
  # d는 평면의 방정식 ax + by + cz + d = 0에서의 d 값
  
  # 평면 위의 점 생성 함수
  generate_plane_points <- function(a, d, range = c(-1, 1), n = 20) {
    max_idx <- which.max(abs(a))
    vars <- setdiff(1:3, max_idx)
    
    grid <- expand.grid(
      seq(range[1], range[2], length.out = n),
      seq(range[1], range[2], length.out = n)
    )
    names(grid) <- c("v1", "v2")
    
    grid$v3 <- (-d - a[vars[1]] * grid$v1 - a[vars[2]] * grid$v2) / a[max_idx]
    
    result <- matrix(0, nrow = n^2, ncol = 3)
    result[, vars] <- as.matrix(grid[, c("v1", "v2")])
    result[, max_idx] <- grid$v3
    
    return(result)
  }
  
  # 점 생성
  points <- generate_plane_points(a, d, n = n)
  
  # 선분 생성 함수
  create_lines <- function(points, n) {
    lines <- matrix(NA, nrow = n * (n-1) * 4, ncol = 3)
    idx <- 1
    
    # 수평 선분
    for (i in 1:n) {
      for (j in 1:(n-1)) {
        lines[idx,] <- points[(i-1)*n + j,]
        lines[idx+1,] <- points[(i-1)*n + j+1,]
        idx <- idx + 2
      }
    }
    
    # 수직 선분
    for (i in 1:n) {
      for (j in 1:(n-1)) {
        lines[idx,] <- points[(j-1)*n + i,]
        lines[idx+1,] <- points[j*n + i,]
        idx <- idx + 2
      }
    }
    
    return(lines)
  }
  
  # 선분 생성 및 그리기
  lines <- create_lines(points, n)
  segments3d(x = lines[,1], y = lines[,2], z = lines[,3], 
             col = col, alpha = alpha, depth_mask = depth_mask)
}

