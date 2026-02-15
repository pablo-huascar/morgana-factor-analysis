
#############################
###                       ###
###         PCA           ### 
###                       ###
#############################
## PCA - COMPONENT ANALYSIS 
#
# R <- Correlation matrix
# q <- number of components
# B <- weights   F = Z*B
# A <- optimal weights for estiamting  Z from Z = Z*B
# f <- explained variance 
# 
pcaR <- function(R, q){
  tmp0 <- svd(R)
  P <- tmp0$u
  D <- diag(tmp0$d)
  Q <- tmp0$v
  
  
  Kq <- Q[, 1:q]
  Dq <- D[1:q, 1:q]
  N <- diag(q)
  
  B <- Kq %*% solve(Dq)^.5%*%N
  A <- t(solve(t(B)%*%R%*%B)%*%t(B)%*%R)
  f <- sum(diag(t(B)%*%R%*%R%*%B%*%(t(B)%*%R%*%B)^(-1)))

  return(list(B=B, A=A, f=f))
}
