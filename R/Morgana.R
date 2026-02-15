# MORGANA 
# 
# Case 2 based on Mulaik's approach and applying EREC doublets 
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
#   L        Loading matrix
#   RESCOV   Variance/Covariance Residuals
#   RESCOR   Correlation residuals
#   

Morgana <- function(X, r, conver, maxiter, corr, doblets){
   
  if(corr==2){
    R <- cor(X)
  } else{
    if(corr==3){
      R <- psych::polychoric(X)$rho
    } else{R <- cov(X)}
  }
  m <- ncol(R)
  
  # Computes the matrix of partial correlations 
  # AI is Guttman's AntiImage matrix, and P the matrix of partial correltaions 
  
  D  <- diag(diag(solve(R))) # solve = a ^(-1)
  tmp0 <- diag(D^(-1/2))
  AI <- diag(tmp0)%*%solve(R)%*%diag(tmp0)
  P <- (2* eye(m))-AI
    
  Pr <- eye(m)
  
  d <- nrow(doblets)
  tmp <- ncol(doblets)
        
  for(i in 1:d){
    Pr[doblets[i,1],doblets[i,2]] <- doblets[i,3] * sign(P[doblets[i,1], doblets[i,2]])
    Pr[doblets[i,2],doblets[i,1]] <- Pr[doblets[i,1], doblets[i,2]]
  }
 
  # initial estimated solution (no iterations)
  
  U <- diag(diag(D^(-1/2)))
  RESI <- U%*%Pr%*%U
  
  # RESI includes PHI2  in the diagonal and the correlated residuals off-diagonal
  
  RED <- R-RESI
  
  # A is the orthogonal pattern in r factors based on the complete residual matrix
  L <- pcaR(RED, r)$A
  
  f <- sum(diag(t(L)%*%L))
  
  fold <- f*10
  
  iter <- 1
  
  while (abs(f - fold) > conver) {
    W <- R - (L %*% t(L))
    UU <- diag(diag(W))
    U <- sqrt(UU)
    RED <- R - (U %*% Pr %*% U)
    result <- pcaR(RED, r)
    B <- result$B
    L <- result$A
    iter <- iter + 1
    fold <- f
    
    if (iter < maxiter) {
      f <- sum(diag(t(L) %*% L))
    }
  }
  
  
  RESCOV <- R-(L%*%t(L))  
  RESI <- solve(U)%*% RESCOV%*%solve(U)
  RESCOR <- eye(m)
  
  for(i in 1:d){
    RESCOR[doblets[i,1],doblets[i,2]] = RESI[doblets[i,1],doblets[i,2]];
    RESCOR[doblets[i,2],doblets[i,1]] = RESI[doblets[i,2],doblets[i,1]];
  }
  results <- list(L=L, RESCOV=RESCOV, RESCOR=RESCOR)
  return(results)
}

