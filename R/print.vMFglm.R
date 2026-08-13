#' @method print vMFglm
#' @export 
print.vMFglm <- function(fit, all=FALSE){
  
  
  if(all){
    
    print(fit)
    
  } else{
    
    print(
      list(beta = fit$beta,
         # beta.scaled = fit$beta2,
         params = c(gamma = fit$gamma,
                    maxit = fit$params$maxit,
                    eps = fit$params$eps,
                    orthogonal = fit$params$orthogonal)
      )
    )
  }
  
  
}

