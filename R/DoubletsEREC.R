
##################################
###                            ###
###      DoubletsEREC          ### 
###                            ###
##################################

# Goodness of bit based on EREC
# 
# IN: 
#   X           Matrix of observed data (N x m)
#   r           Number of factors
#   cor         1 <- covariance
#               2 <- pearson
#               3 <- polychoric
# 
# OUT: 
#   doublets    Pairs of doblets to be free in the exploratory MORGANA
#   nd          Number of pairs


DoubletsEREC <- function(X, r, cor, PairsDoublets){
  
  N <- nrow(X)
  m <- ncol(X)
  conver <- 0.00001
  PairsDoublets <- as.matrix(PairsDoublets)
  
  Ledermann <- floor((2*m+1-sqrt(8*m+1))/2)
  df <- 1/2 * ((m-r) * (m-r+1)) - m
  # max number of doblets
  nd <- Ledermann - r
  if((df - nd) < 1){
    nd <- df - 1
  }
  
  
  R <- switch(cor, "1" = cov(X), "2" = cor(X), "3" = psych::polychoric(X)$rho)
  
  if(cor>1){
    Ac <- my_minres_hc(R, r, conver, N)$Loadings
  } else{
    Ac <- as.matrix(my_minres(R, r, conver, N)$Loadings)
  }
  
  
  PairOut <- zeros(2,1)
  m_Core <- m-2
  Core <- zeros(m_Core, 1)
  RESFREE <- eye(m)
  IPSIs <- c()
  
  
  PPP <- matrix(0, nrow = 0, ncol = nrow(PairOut))
  
  
  d <- 0
  columna <- fila<- c()
  
  for (i in 1:(m-1)) {
    for (j in (i+1):m) {
      PairOut <- matrix(0, nrow = 2, ncol = 1)
      PairOut[1, 1] <- i
      PairOut[2, 1] <- j
      PPP <- rbind(PPP, t(PairOut))
      k <- 1
      Core <- integer(m-2)
      for (h in 1:m) {
        if (h != i && h != j) {
          Core[k] <- h
          k <- k + 1
        }
      }
      
      # inestable 
      RCore <- R[Core, Core]
      
      if(cor>1){
        Ac <- my_minres_hc(RCore, r, conver, N)$Loadings
      } else{
        Ac <- as.matrix(my_minres(RCore, r, conver, N)$Loadings)
      }
      
      # opción estable
      
      # Ac <- A[Core, ]
      
      
      lam <- zeros(2,r)
      for(h in 1:2){
        SelExtension <- c(PairOut[h,1], t(Core))
        R_Ex <- R[SelExtension, SelExtension]
        vecr <- R_Ex[2:(m_Core+1), 1]
        lam[h,] <- (solve(t(Ac)%*%Ac))%*%t(Ac)%*%vecr
      }
      # Residual covariances 
      
      RESFREE[PairOut[1,1], PairOut[2,1]] <- R[PairOut[1,1], PairOut[2,1]]-crossprod(lam[1,], lam[2,])
      RESFREE[PairOut[2,1], PairOut[1,1]] <- RESFREE[PairOut[1,1], PairOut[2,1]]
      d <- d+1
    }
    
    # Collection of IPSIS for correcting EPC1
    
    Li <- rbind(Ac, lam)
    Core <- as.matrix(Core)
    temp <- as.matrix(order(c(Core, PairOut)))
    Li <- Li[temp,]
    
    
    IPSIi <- as.matrix(diag(1, nrow=m)- diag(diag(Li%*%(t(Li)))))^(-0.5)
    diagIPSIi <- t(as.matrix(diag(IPSIi)))
    IPSIs<- rbind(IPSIs, diagIPSIi)
    
  }
  Ips <- matrix(0, nrow=nd, ncol=m)
  for(i in 1:nd){
    for(j in 1:m){
      if(IPSIs[i,j]< 0.001){
        IPSIs[i,j] <- max(IPSIs)
      } else{IPSIs[i,j] <- IPSIs[i,j]}
    }
  }
  
  IPSI <- diag(colmin(IPSIs))
  EREC <- IPSI%*%RESFREE%*%IPSI
  
  
  d <- nrow(PairsDoublets)
  tmp <- ncol(PairsDoublets)
  
  if(d < nd){
    nd <- d
  }
  
  doublets <- zeros(nd, 3)
  for(i in 1:nd){
    doublets[i,1] <- PairsDoublets[i,1]
    doublets[i,2] <- PairsDoublets[i,2]
    doublets[i,3] <- EREC[PairsDoublets[i,1],PairsDoublets[i,2]]
  }
  
  

  result <- list(doublets=doublets, nd=nd)
  return(result)
  
}
