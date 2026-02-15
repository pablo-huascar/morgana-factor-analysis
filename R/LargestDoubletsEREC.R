#########################################
###                                   ###
###      LargestDoubletsEREC          ### 
###                                   ###
#########################################

# Goodness of Fit based on EREC
# 
# IN: 
#   X           Matrix of observed data (N x m)
#   r           Number of factors
#   cor         1 <- covariance
#               2 <- pearson
#               3 <- polychoric
# 
# OUT: 
#   doublets    Pairs of doublets to be free in the exploratory MORGANA
#   nd          Number of pairs


LargestDoubletsEREC <- function(X, r, cor){
  
  N <- nrow(X)
  m <- ncol(X)
  conver <- 0.00001
  
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

      RCore <- R[Core, Core]
      if(cor>1){
        Ac <- my_minres_hc(RCore, r, conver, N)$Loadings
      } else{
        Ac <- as.matrix(my_minres(RCore, r, conver, N)$Loadings)
      }
       
      lam <- zeros(2,r)
      for(h in 1:2){
        SelExtension <- c(PairOut[h,1], t(Core))
        R_Ex <- R[SelExtension, SelExtension]
        vecr <- R_Ex[(2:(m_Core+1)), 1]
        vecr <- as.matrix(vecr, ncol=1)
        lam[h,] <- (solve(t(Ac)%*%Ac))%*%t(Ac)%*%vecr
      }
      
      # Collection of IPSIS for correcting EPC1
      Li <- rbind(Ac, lam)
      Core <- as.matrix(Core)
      temp <- as.matrix(order(c(Core, PairOut)))
      Li <- Li[temp,]
   
      IPSIi <- as.matrix(diag(1, nrow=m)- diag(diag(Li%*%(t(Li)))))^(-0.5)
      diagIPSIi <- t(as.matrix(diag(IPSIi)))
      IPSIs<- rbind(IPSIs, diagIPSIi)      

      # Residual covariances 
      RESFREE[PairOut[1,1], PairOut[2,1]] <- R[PairOut[1,1], PairOut[2,1]]-crossprod(lam[1,], lam[2,])
      RESFREE[PairOut[2,1], PairOut[1,1]] <- RESFREE[PairOut[1,1], PairOut[2,1]]
      d <- d+1
    }
  }

  for(i in 1:nd){
    for(j in 1:m){
      if(IPSIs[i,j]< 0.001){
        IPSIs[i,j] <- max(IPSIs)
      } else{IPSIs[i,j] <- IPSIs[i,j]}
    }
  }
  
  IPSI <- diag(colmin(IPSIs))
  EREC <- IPSI%*%RESFREE%*%IPSI
  
  total_pairs <- m*(m-1)/2
  EREC_row <- c(total_pairs)
  Pair <- zeros(total_pairs, 2)

  h <- 1
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      Pair[h,1] <- i
      Pair[h,2] <- j
      EREC_row[h] <- abs(EREC[i,j])
      h <- h+1
    }
  }
  J <- order(EREC_row, decreasing = TRUE)
  
  doublets <- zeros(nd, 3)
  for(i in 1:nd){
    doublets[i, 1] <- Pair[J[i],1]
    doublets[i, 2] <- Pair[J[i],2]
    doublets[i, 3] <- EREC_row[J[i]]
  }
  
  result <- list(doublets=doublets, nd=nd)
  return(result)
  
}

