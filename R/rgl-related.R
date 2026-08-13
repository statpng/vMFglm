transform3d <- function(obj,matrix,...) rotate3d.mesh3d(obj, matrix = matrix, angle = NA, x = NA, y = NA, z = NA)

translate3d <- function(obj,x,y,z,...) {
  if (is.matrix(obj)) n <- dim(obj)[1]
  else n <- 1
  
  if (length(obj) == 3 || (is.matrix(obj) && dim(obj)[2] == 3))
    return(obj + cbind(rep(x,n), rep(y,n), rep(z,n)))
  else if (length(obj) == 4 || (is.matrix(obj) && dim(obj)[2] == 4))
    return(obj %*% translationMatrix(x,y,z))
  else stop("Unsupported object for translation")
}


scale3d <- function(obj,x,y,z,...) {
  if (is.matrix(obj)) n <- dim(obj)[1]
  else n <- 1
  
  if (length(obj) == 3 || (is.matrix(obj) && dim(obj)[2] == 3))
    return(obj * cbind(rep(x,n), rep(y,n), rep(z,n)))
  else if (length(obj) == 4 || (is.matrix(obj) && dim(obj)[2] == 4))
    return(obj %*% scaleMatrix(x,y,z))
  else stop("Unsupported object for scaling")
}



rotate3d <- function(obj,angle,x,y,z,matrix,...) {
  if (length(obj) == 3 || (is.matrix(obj) && dim(obj)[2] == 3))
    return(asEuclidean(asHomogeneous(obj) %*% rotationMatrix(angle,x,y,z,matrix)))
  else if (length(obj) == 4 || (is.matrix(obj) && dim(obj)[2] == 4))
    return(obj %*% rotationMatrix(angle,x,y,z,matrix))
  else stop("Unsupported object for rotation")
}    


rotate3d.mesh3d <- function( obj,angle,x,y,z,matrix, ... ) {
  obj$vb <- t(rotate3d(t(obj$vb), angle, x, y, z, matrix))
  if ( !is.null(obj$normals) ) {
    normals <- obj$normals  # 4xn matrix
    if ( missing(matrix) ) 
      normals <- rotate3d(t(normals), angle, x, y, z) # nx4
    else {
      if (nrow(matrix) == 4) matrix[4,1:3] <- 0
      if (ncol(matrix) == 4) matrix[1:3,4] <- 0
      normals <- rotate3d(t(normals), angle, x, y, z, t(solve(matrix))) # nx4
      # Reverse the normals if the transformation
      # changes the orientation of the mesh
      if (det(matrix) < 0)
        normals[, 1:3] <- -normals[, 1:3]  # nx4
    }
    normals <- asEuclidean(normals)  # nx3
    normals <- normals/sqrt(rowSums(normals^2))
    obj$normals <- rbind( t(normals), 1 ) # 4xn
  }
  return(obj)                            
}  




rotationMatrix <- function(angle,x,y,z,matrix) {
  if (missing(matrix)) {
    if (angle == 0) return(identityMatrix())
    
    u <- c(x,y,z)/sqrt(x^2+y^2+z^2)
    cosa <- cos(angle)
    sina <- sin(angle)
    matrix <- (1-cosa)*outer(u,u)
    matrix <- matrix + diag(3)*cosa
    matrix[1,2] <- matrix[1,2] - sina*u[3]
    matrix[1,3] <- matrix[1,3] + sina*u[2]
    matrix[2,1] <- matrix[2,1] + sina*u[3]
    matrix[2,3] <- matrix[2,3] - sina*u[1]
    matrix[3,1] <- matrix[3,1] - sina*u[2]
    matrix[3,2] <- matrix[3,2] + sina*u[1]
  }
  if (identical(all.equal(dim(matrix), c(3,3)), TRUE)) 
    matrix <- cbind(rbind(matrix,c(0,0,0)),c(0,0,0,1))
  return(matrix)    
}



translationMatrix <- function(x,y,z) {
  result <- diag(4)
  result[4,1:3] <- c(x,y,z)
  result
}


identityMatrix <- function() diag(nrow=4)

scaleMatrix <- function(x,y,z) diag(c(x,y,z,1))






# Coordinate conversions

asHomogeneous <- function(x) {
  if (is.matrix(x)) {
    if (ncol(x) == 3) return(cbind(x,1))
    if (ncol(x) == 4) return(x)
  } else {
    if (length(x) %% 3 == 0) 
      return(cbind(matrix(x, ncol = 3, byrow = TRUE),1))
    if (length(x) %% 4 == 0)
      return(matrix(x, ncol = 4, byrow = TRUE))
  }
  stop("Don't know how to convert x")
}

asHomogeneous2 <- function(x) {
  if (is.matrix(x)) {
    if (nrow(x) == 3) return(rbind(x,1))
    if (nrow(x) == 4) return(x)
  } else {
    if (length(x) %% 3 == 0) 
      return(rbind(matrix(x, nrow = 3), 1))
    if (length(x) %% 4 == 0)
      return(matrix(x, nrow = 4))
  }
  stop("Don't know how to convert x")
}

asEuclidean <- function(x) {
  if (is.matrix(x)) {
    if (ncol(x) == 4) return(x[, 1:3, drop = FALSE]/x[, 4])
    if (ncol(x) == 3) return(x)
  } else {
    if (length(x) %% 4 == 0) 
      return(asEuclidean(matrix(x, ncol = 4, byrow = TRUE)))
    if (length(x) %% 3 == 0)
      return(matrix(x, ncol = 3, byrow = TRUE))
  }
  stop("Don't know how to convert x")
}

asEuclidean2 <- function(x) {
  if (is.matrix(x)) {
    if (nrow(x) == 4) return(rbind(x[1,]/x[4,], x[2,]/x[4,], x[3,]/x[4,]))
    if (nrow(x) == 3) return(x)
  } else {
    if (length(x) %% 4 == 0) 
      return(asEuclidean2(matrix(x, nrow = 4)))
    if (length(x) %% 3 == 0)
      return(matrix(x, nrow = 3))
  }
  stop("Don't know how to convert x")
}








# shapelist3d -------------------------------------------------------------

scale3d.shapelist3d <- function( obj, x, y, z, ... ) {
  structure(lapply( obj, function(item) scale3d(item, x,y,z,...) ),
            class = class(obj))
}


# mesh3d ------------------------------------------------------------------


scale3d.mesh3d <- function( obj, x, y, z, ... ) {
  obj$vb <- t(scale3d(t(obj$vb), x, y, z))
  if ( !is.null(obj$normals) ) {
    obj$normals <- scale3d(t(obj$normals), 1/x, 1/y, 1/z)
    obj$normals <- t( obj$normals/sqrt(apply(obj$normals[,1:3]^2, 1, sum)) )
    if (nrow(obj$normals) == 4)
      obj$normals[4,] <- 1
  }
  return(obj)
}




translate3d.mesh3d <- function( obj, x, y, z, ... ) {
  obj$vb <- t(translate3d(t(obj$vb), x, y, z))
  return(obj)                            
}  
