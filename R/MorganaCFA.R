
# MORGANA CFA
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
#   doblets  Matrix of doblets 


MorganaCFA <- function(X, r, conver, maxiter,cor,doblets) {
  N <- nrow(X)
  m <- ncol(X)
  
  doublets_erec <- DoubletsEREC(X,r,cor,doblets)$doublets
  doublets <- doublets_erec
  
  res <- Morgana(X,r,conver,maxiter, cor, doublets_erec)
  
  A <- res$L
  RESCOV <- res$RESCOV
  RESCOR <- round(res$RESCOR,3)
  
  results <- list(doublets=doublets, A=A, RESCOV=RESCOV, RESCOR=RESCOR)
  
  return(results)
}

