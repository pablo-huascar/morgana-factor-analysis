
# MORGANA EFA
# 
# 
# 
# IN: 
#   X        Matrix of observed data (N x m)
#   r        Number of factors to be extrated
#   conver   Convergence value for MINRES
#   maxiter  Maximun iterations
#   cor      1 -> covariance
#            2 -> pearson
#            3 -> polychoric
# 
# OUT:
#   A        Loading matrix
#   RESCOV   Variance/Covariance Residuals
#   RESCOR   Correlation residuals
#   doublets  Matrix of doublets 


MorganaEFA <- function(X, r, conver, maxiter,cor) {
 
  N <- nrow(X)
  m <- ncol(X)
  
  doublets <- LargestDoubletsEREC(X,r,cor)$doublets
  res <- Morgana(X,r,conver,maxiter, cor, doublets)
  
  A <- res$L
  RESCOV <- res$RESCOV
  RESCOR <- res$RESCOR
  
  results <- list(doublets=doublets, A=A, RESCOV=RESCOV, RESCOR=RESCOR)
  
  return(results)
}

