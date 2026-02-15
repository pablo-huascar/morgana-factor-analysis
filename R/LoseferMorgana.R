

#########################################
###                                   ###
###           LoseferMorgana          ### 
###                                   ###
#########################################

LoseferMorgana <- function(S, A, N, RESCOR, nd, mchi, sxchi, a_n, b1_n, b2_n, a_f, b_f, c_f, d_f, display){
  
  # STEP1: c2: uncorrected chi-square; compute d = degrees of freedom
  
  m <- nrow(A)
  r <- ncol(A)
  
  
  # STEP4: Residuals
  AA <- A%*%t(A)
  PSI2 <- eye(m)-diag(diag(AA))
  PSI <- sqrt(PSI2)
  RESCORZ <- RESCOR - eye(m)
  RES <- S-(AA + PSI %*% RESCORZ %*% PSI)
  sum <- 0
  
  
  for(i in 1:(m-1)){
    for(j in (1+i):m){
      sum <- sum + RES[i,j]%*%RES[i,j]
    }
  }
  FF <- sum
  c2 <- FF*(N-1)
  d <- 1/2 * ((m-r)* (m-r+1)) - m - nd
  
  v <- mchi
  vv <- sxchi
  
  logc2 <- c2^(1/3)
  z <- (logc2 - v)/vv 
  zn <- a_n +(b1_n*z)+(b2_n*z*z)
  
  Um <- (a_f + zn*(b_f +zn*(c_f+zn*d_f)))*sqrt(2*d)+d
  
  if(Um > d+300*sqrt(d*2)){
    Um <- d+300*sqrt(d*2)+runif(1,1)*mchi*4-(mchi*2)
  }
  
  sum <- 0
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      sum <- sum + S[i,j]*S[i,j]
    }
  }
  
  FF <- sum
  C0 <- FF*(N-1)
  d0 <- 1/2 * m * (m-1)
  
  CFI <- abs(((C0 - d0) - (Um - d))/(C0-d0))
  if(CFI > .999){
    CFI <- .999
  }else{
    if(CFI < 0){
      CFI <- 0.01
    }
  }
  
  Mo <- C0/d0
  Mk <- Um/d
  TL <- (Mo-Mk)/(Mo-1)
  
  if(TL > .999){
    TL <- .999
  } else{
    if(TL < 0){
      TL <- CFI
    }
  }
  
  Fo <- Um/(N-1)
  FFo <- Fo - (d/(N-1))
  if(FFo < 0.0){
    FFo <- 0.0
  }
  RMSEA <- sqrt(FFo/d)
  if(RMSEA < 0){
    RMSEA <- 0.001
  }
  
  if(display==1){
    buff <- sprintf("Chi Square adjusted by Losefer with %i degrees of freedom = %.3f", d, Um)
    cat(buff)
    cat("\n")
    buff <- sprintf("Goodness of Fit (Null Model) = %.2f", C0)
    cat(buff)
    cat("\n")
    buff <- sprintf("Degrees of freedom =%.0f", d0)
    cat(buff)
    cat("\n")
    buff <- sprintf("RMSEA (Root Mean Square Error of Approximation)= %.3f", RMSEA)
    cat(buff)
    cat("\n")
    buff <- sprintf("CFI (Comparative Fit Index)= %.3f", CFI)
    cat(buff)
    cat("\n")
    buff <- sprintf("NNFI (Coeficiente de Tucker y Lewis)= %.3f", TL)
    cat(buff)
    cat("\n")
  }
  
  return(list(Um=Um, d=d, d0=d0, C0=C0, RMSEA=RMSEA, CFI=CFI, TL=TL))
}
#########################################
