#########################################
###                                   ###
###      MINRES HEYWOOD CASE          ### 
###                                   ###
#########################################
## MINRES EXTRACTION BY Zegers & ten Berge (1983)
# C -> direct, covariance or correlation matrix
# q -> number of components
# conv -> convergence value
# N -> Number of subjects
# A -> Loadings
# Um -> Index of fit
# m <- number of items 
###############################################################################
my_minres_hc <- function(C, q, conv, N){
 
  # r number of factors
  A <- pcaR(C,q)$A
  
  n <- dim(A)[1]
  m <- dim(A)[2]
  
  Co <- C-diag(diag(C))
  
  g1 <- t(Co - (A%*%t(A)) - diag(diag(A%*%t(A))))
  g2 <- Co - (A%*%t(A)) - diag(diag(A%*%t(A)))
  g <- sum(diag(g1 %*% g2))
  gold = g*2
  
  iter <- 0 
  
  while(abs(gold - g) > conv){
    for(i in 1:n){
      Ai <- A
      for (j in 1: m){
        Ai[i,j] <- 0.0
      }

      ci <- as.matrix(c(Co[,i]), ncol=1)
      ai <- (solve(t(Ai)%*% Ai))%*%(t(Ai)%*%(ci))
    
      ai <- (solve(t(Ai)%*% Ai))%*%(t(Ai)%*%(ci))
    
     #HEYWOOD CASE
      condition <- t(ai)%*%ai
      if(condition > 1.0){
        ai <- BERGE(Ai, ci)
      }
    
      for(j in 1:m){
        A[i,j] <- ai[j,]
      }
    }
    
    gold <- g
    g1 <- t(Co - (A%*%t(A)) - diag(diag(A%*%t(A))))
    g2 <- Co - (A%*%t(A)) - diag(diag(A%*%t(A)))
    g <- sum(diag(g1 %*% g2))
    iter <- iter + 1
    if(iter > 1000) {
      gold <- g
    }    
    
  }
  
  Re <- (A%*%t(A))+ diag(diag(C)) - diag(diag(A%*%t(A)))
  Um = (N-1)*(log(det(Re))-log(det(C)))
  if(Um < 0){
    Um <- 0
  }
  return(list("Loadings"=A, "Index of Fit"=Um))
}