# fit_fleishman_from_sk
fit_fleishman_from_sk <- function(skew, kurt){

max_iter<-100

  A <- NaN
  B <- NaN
  C <- NaN
  
  op <- (-1.13168 + 1.58837 * (skew^2))
  
  if(kurt < op){
    A <- NaN
    B <- NaN
    C <- NaN
  } else{
    b <- 0.95357 - 0.05679 * kurt + 0.03520 * skew^2 + 0.00133 * kurt^2
    c <- 0.10007 * skew + 0.00844 * skew^3
    d <- 0.30978 - 0.31655 * b
    b2 <- b * b
    c2 <- c * c
    d2 <- d * d
    bd <- b * d
    varf <- b2 + 6*bd + 2*c2 + 15*d2
    skewf <- 2 * c * (b2 + 24*bd + 105*d2 + 2)
    kurtf <- 24 * (bd + c2 * (1 + b2 + 28*bd) + d2 * (12 + 48*bd + 141*c2 + 225*d2))
    f <- matrix(0,3,1) 
    J <- matrix(0,3,3)         
    f[1,1] <- varf - 1
    f[2,1] <- skewf - skew
    f[3,1] <- kurtf - kurt
    for (i in 1:max_iter){
        b2 <- b * b
        c2 <- c * c
        d2 <- d * d
        bd <- b * d
        J[1,1] <- 2*b + 6*d
        J[1,2] <- 4*c
        J[1,3] <- 6*b + 30*d
        J[2,1] <- 4*c * (b + 12*d)
        J[2,2] <- 2 * (b2 + 24*bd + 105*d2 + 2)
        J[2,3] <- 4 * c * (12*b + 105*d)
        J[3,1] <- 24 * (d + c2 * (2*b + 28*d) + 48 * d^3)
        J[3,2] <- 48 * c * (1 + b2 + 28*bd + 141*d2)
        J[3,3] <- 24 * (b + 28*b * c2 + 2 * d * (12 + 48*bd + 141*c2 + 225*d2) + d2 * (48*b + 450*d))
        delta <- -(solve(t(J) %*% J) %*% t(J) %*% f)  
        b <- delta[1,1] + b
        c <- delta[2,1] + c
        d <- delta[3,1] + d
        b2 <- b * b
        c2 <- c * c
        d2 <- d * d
        bd <- b * d
        varf <- b2 + 6*bd + 2*c2 + 15*d2
        skewf <- 2 * c * (b2 + 24*bd + 105*d2 + 2)
        kurtf <- 24 * (bd + c2 * (1 + b2 + 28*bd) + d2 * (12 + 48*bd + 141*c2 + 225*d2))
        f[1,1] <- varf - 1
        f[2,1] <- skewf - skew
        f[3,1] <- kurtf - kurt
    }
    A <- b
    B <- c
    C <- d
  }
  lista <- list(A=A, B=B, C=C)
  return(lista)
}
