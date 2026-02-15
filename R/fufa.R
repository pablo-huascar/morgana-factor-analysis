
##########################
###                    ###
###      FUFA1         ### 
###                    ###
##########################


#Free Unrestricted Factor Analysis (FUFA)
#
# Case 1
#
# IN:
#   X         Matrix of observed data (N x m)
#   r         Number of factors to be extracted
#   doblets   Matrix (p x 2) with p pairs of variables with correlatd errors
#   conver    Convergence value for MINRES
#   cor       1 -> covariance
#             2 -> pearson
#             3 -> polychoric
#
# OUT:
#   L         Loading matrix
#   RESCOV    Variance/Covariance residuals
#   RESCOR    Correlation residuals  


fufa1 <- function(X, r, doblets, conver, cor) {
  doblets <- as.matrix(doblets)
  N <- nrow(X)
  m <- ncol(X)
  
  # Calculate the correlation or covariance matrix based on the specified method
  switch(cor,
         2, {
           R <- cor(X)
         },
         3, {
           R <- psych::polychoric(X)$rho
         },
         {
           R <- cov(X)
         }
  )
  
  p <- nrow(doblets)
  
  # Extract the pairs of variables
  Pairs <- c(t(doblets))
  
  # Identify the core variables that are not in the pairs
  Core <- rep(-1, m) # Inicializar 'Core' con -1
  h <- 1
  
  for (i in 1:m) {
    campana <- 0
    for (j in 1:(p * 2)) {
      if (i == Pairs[j]) {
        campana <- 1
        break  # Salir del bucle anidado si encontramos la variable en los pares
      }
    }
    
    if (campana == 0) {
      Core[h] <- i
      h <- h + 1
    }
  }
  Core <- as.matrix(ClearBuits(Core),ncol=1)
  m_core <- nrow(Core)
  
  
  
  PairsSort <- zeros(p,2)
  rpar <- zeros(p,1)
  J <- rep(0, m)
  # Sort the pairs based on the partial correlations of the first item
  for (h in 1:p) {
    
    J[1] <- doblets[h, 1]
    J[2] <- doblets[h, 2]
    i <- 3
    
    for (j in 1:(p * 2)) {
      if ((Pairs[j] != J[1]) && (Pairs[j] != J[2])) {
        J[i] <- Pairs[j]
        i <- i + 1
      }
    }
    for (j in Core) {
      J[i] <- j
      i <- i + 1
    }
    
    # Calculate the partial correlations
    Rh <- R[J, J]
    Ra <- PartCor(Rh)
    
    # Order the pair based on the average of squared partial correlations
    if (mean(Ra[3:m, 1] * Ra[3:m, 1]) < mean(Ra[3:m, 2] * Ra[3:m, 2])) {
      PairsSort[h, 1] <- J[1]
      PairsSort[h, 2] <- J[2]
      rpar[h] <- mean(Ra[3:m, 1] * Ra[3:m, 1])
    } else {
      PairsSort[h, 1] <- J[2]
      PairsSort[h, 2] <- J[1]
      rpar[h] <- mean(Ra[3:m, 2] * Ra[3:m, 2])
    }
  }
  
  tmp <- order(rpar)
  PairsSort <- PairsSort[tmp, ]
  
  PairIn <- zeros(m,1)-1
  PairProhibit <- zeros(m,1)
  PairIn[1,1]<- PairsSort[1,1]
  PairProhibit[1,1] <- PairsSort[1,2]
  
  h <- 1
  
  # Create the PairIn list by considering the pairs and prohibitions
  for (i in 2:p) {
    campana <- 1
    
    # Check if the first element can enter the PairIn list
    for (j in 1:h) {
      if (PairsSort[i, 1] == PairIn[j,] || PairsSort[i, 2] == PairIn[j,] || PairsSort[i, 1] == PairProhibit[j,]) {
        campana <- 0
      }
    }
    
    # If it can enter, add it to the list
    if (campana == 1 ) {
      h <- h + 1
      PairIn[h,] <- PairsSort[i, 1]
      PairProhibit[h,] <- PairsSort[i, 2]
    } else {
      campana <- 0
      
      # Check if the second element can enter the PairIn list
      for (j in 1:h) {
        if (PairsSort[i, 2] == PairIn[j,] || PairsSort[i, 2] == PairProhibit[j,]) {
          campana <- 1
        }
      }
      
      # If it can enter, add it to the list
      if (campana==0) {
        h <- h + 1
        PairIn[h,1] <- PairsSort[i, 2]
        PairProhibit[h,1] <- PairsSort[i, 1]
      }
    }
  }
  
  PairIn <- ClearBuits(PairIn)
  m_IN <- length(PairIn)
  
  # Create the PairOut list for the variables that are not in PairIn or Core
  PairOut <- zeros(m,1)-1
  m_OUT <- length(PairOut)
  
  h <- 1
  for(i in 1:m){
    campana <- 0
    # item is a core?
    for(j in 1:m_core){
      if(Core[j,] == i){
        campana <- 1
      }
    }
    if(campana ==0){
      # item is in PairIn
      for(j in 1:m_IN){
        if(PairIn[j]==i){
          campana <- 1
        }
      }
    }
    # any list -> PairOut
    if(campana ==0){
      PairOut[h,1]<- i
      h <- h+1
    }
  }
  
  
  PairOut <- ClearBuits(PairOut)
  m_OUT <- length(PairOut)
  
  
  
  # Step 1: Create the nuclear matrix without pairs and obtain the canonical solution
  
  SelCore <- c(Core, PairIn)
  Xcore <- X[, SelCore]
  Rcore <- cor(Xcore)
  CANCORE <- my_minres(Rcore, r, conver, N)$Loadings
  
  
  # Step 2: Estimate the remaining structural parameters by extension
  P <- CANCORE[1:m_core,]
  CANEXT <- zeros(m_OUT,r) 
  
   
  
  
  
  for (i in 1:m_OUT) {
    SelExtension <- c(PairOut[i], Core)
    R_Ex <- R[SelExtension, SelExtension]
    vecr <- R_Ex[2:(m_core + 1), 1]
    lam <- solve(t(P)%*% P)%*%t(P)%*%(vecr)
    CANEXT[i, ] <- lam
  }
  
  # Step 3: Construct the complete pattern and calculate the uniqueness matrix
  L <- rbind(CANCORE, CANEXT)
  J <- order(c(Core, PairIn, PairOut))
  L <- L[J, ]
  
  PSI2 <- diag(diag(L %*% t(L)))
  PSI <- sqrt(PSI2)
  IPSI <- solve(PSI)
  
  # Step 4: Calculate residuals
  RESCOV <- R - (L %*% t(L))
  RESI <- IPSI %*% RESCOV %*% IPSI
  RESCOR <- diag(m)
  
  for (i in 1:p) {
    RESCOR[doblets[i, 1], doblets[i, 2]] <- RESI[doblets[i, 1], doblets[i, 2]]
    RESCOR[doblets[i, 2], doblets[i, 1]] <- RESI[doblets[i, 2], doblets[i, 1]]
  }
  
  return(list(L = L, RESCOV = RESCOV, RESCOR = RESCOR, RESI = RESI))
}
