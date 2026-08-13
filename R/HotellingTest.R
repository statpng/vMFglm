#' @export HotellingStat
HotellingStat <- function(Y, group, euclidean=FALSE){
  lev <- levels(as.factor(group))
  p <- dim(Y)[2]
  
  if(euclidean){
    Y1 <- Y[group == lev[1],,drop=F]
    Y2 <- Y[group == lev[2],,drop=F]
    
    
    euc_tmp <- euclideanization( t(Y1), t(Y2) )
    
    Y_euc <- matrix(NA, nrow(Y), ncol(Y)-1)
    Y_euc[group == lev[1],] <- euc_tmp$euclideanG1
    Y_euc[group == lev[2],] <- euc_tmp$euclideanG2
    
    Y <- Y_euc
  }
  
  
  Y1 = Y[group == lev[1], ,drop=F]
  Y2 = Y[group == lev[2], ,drop=F]
  
  n1=nrow(Y1);  n2=nrow(Y2)
  
  D <- colMeans(Y1) - colMeans(Y2)
  COV <- ( (n1-1) * cov(Y1) + (n2-1) * cov(Y2) ) / (n1+n2-2)
  TestStat <- c( t(D) %*% MASS::ginv(COV*(1/n1+1/n2)) %*% D )
  
  n <- n1+n2-1
  F_value<-((n-p)/(p*(n-1)))*TestStat
  df1<-p
  df2<-n-p
  p_value<-1-pf(F_value,df1,df2)
  
  result <- list(Y=Y, group=group, T2=TestStat, Fvalue=F_value, pvalue=p_value)
  result
}


#' @export HotellingPermTest
HotellingPermTest <- function(Y, group, nperm=1000, euclidean=FALSE){
  result <- NULL
  result$Y <- Y
  result$group <- group
  
  fit.Hotelling <- HotellingStat(Y, group, euclidean=euclidean)
  Y <- fit.Hotelling$Y
  group <- fit.Hotelling$group
  
  TestStat <- fit.Hotelling$T2
  result$TestStat <- TestStat
  
  res.perm <- sapply(1:nperm, function(b){
    HotellingStat(Y, sample(group), euclidean=FALSE)$T2
  })
  
  pvalue <- ( sum( TestStat <= res.perm ) + 1 ) / (nperm + 1)
  
  result$nperm <- nperm
  result$perm <- res.perm
  result$pvalue <- pvalue
  
  result
}
