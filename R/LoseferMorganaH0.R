



###########################################
###                                     ###
###           LoseferMorganaH0          ### 
###                                     ###
###########################################



# empirical distribution for H0
# 
# FINAL TRANSFORMTION of a value xi taking the parameters returned by the function
# 
# STEP0: calculate degrees of freedom
# 
# STEP1: log
# 
# STEP2: z
# 
# STEP3:  normal
# 
# STEP4:  chiq
# 
# 

LoseferMorganaH0 <- function(A, RESCOR, doblets, N, montecarlo) {
  conver <- 0.0001
  
  m <- nrow(A)
  r <- ncol(A)
  nd <- nrow(doblets)
  
  AA <- A %*% t(A)
  PSI2 <- diag(m) - diag(diag(AA))
  PSI <- sqrt(PSI2)
  RESCORZ <- RESCOR - diag(diag(RESCOR))
  Rr <- (AA + PSI %*% RESCORZ %*% PSI)
  Rr <- Rr - diag(diag(Rr)) + diag(m)
  df <- 1/2 * ((m - r) * (m - r + 1)) - m - nd
  
  EV <- eigen(Rr)
  U <- EV$vectors
  V <- diag(EV$values)
  if (min(diag(V)) < 0.001) {
    Rr <- Rr + diag(rep(N, m))
  }
  
  NP <- 100000
  Zi <- matrix(rnorm(NP * m), ncol = m)
  L <- chol(Rr)
  XP <- Zi %*% L
 
   
  DCHI <- zeros(montecarlo,1)
  h <- 1
  while (h < montecarlo) {
    J <- order(runif(NP))
    Xi <- XP[J[1:N], ]
    Ri <- cor(Xi)
    
    result <- fufa1(Xi, r, doblets, conver, 2)
    Ai <- result$L
    RESCOVi <- result$RESCOV
    RESCORi <- result$RESCOR
    
    AAi <- Ai %*% t(Ai)
    PSI2 <- eye(m) - diag(diag(AAi))
    PSI <- sqrt(abs(PSI2))
    RESCORZi <- RESCORi - eye(m)
    RESi <- Ri - (AAi + PSI %*% RESCORZi %*% PSI)
    sum <- 0
    for(i in 1:(m-1)){
      for(j in (1+i):m){
        sum <- sum+RESi[i,j]*RESi[i,j]
      }
    }
   
    
    Fi <- sum
    c2i <- Fi * (N - 1)
    
    if (c2i > 0.01) {
      DCHI[h,] <- c2i
      h <- h + 1
    }
  }
  
  md <- median(DCHI)
  sx <- sd(DCHI)
  sum <- 0
  sumq <- 0
  n <- 0
  h <- 1
  x <- c()
  xi <- 0
  while (h < montecarlo) {
    if (DCHI[h,] < (md + 5.5 * sx)) {
      xi <- DCHI[h,] ^ (1/3)
      sum <- sum+xi
      sumq <- sumq +xi*xi
      x <- c(x, xi)
      n <- n + 1
    }
    h <- h + 1
  }
  res <- Trans2Normal(x)
  v <- res$md
  vv <- res$sx
  a_n <- res$a
  b1_n <- res$b1
  b2_n <- res$b2
  
  skew <- sqrt(8/df)
  kurt <- 12/df
  coeffs <- fit_fleishman_from_sk(skew, kurt)
  a_f <- -coeffs$B
  b_f <- coeffs$A
  c_f <- coeffs$B
  d_f <- coeffs$C
  
  return(list(v = v, vv = vv, a_n = a_n, b1_n = b1_n, b2_n = b2_n, a_f = a_f, b_f = b_f, c_f = c_f, d_f = d_f))
}

