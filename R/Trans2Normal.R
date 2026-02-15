#TRANSNORMAL     Transform an initial variable to have the first four moments
#                of a standard normal distribution
#
#                The initial variable has only positive values and 
#                its distribution is rightly skewed 
#   x is expected to be previously transformed as x = log(Y)
Trans2Normal <- function(x) {
max_iter <- 100
  
n <- nrow(as.matrix(x))
md <- mean(x)
sx <- sqrt(var(x))
z <- (x-md)/sx
z2 <- z*z
z3 <- z2*z
z4 <- z3*z
zp1 <- (z2-matrix(1,n,1))
zp2 <- zp1*zp1
zp3 <- zp2*zp1
zp4 <- zp2*zp2
t3 <- mean(z3)
t4 <- mean(z4) 
pt2 <- mean(zp2)
pt3 <- mean(zp3) 
pt4 <- mean(zp4)
c1 <- cov(z,zp1) 
c2 <- cov(z2,zp1) 
c3 <- cov(z,zp2) 
c4 <- cov(z3,zp1) 
c5 <- cov(z,zp3) 
c6 <- cov(z2,zp2) 
b1 <- 0.95
b2 <- 0.15

B <- matrix(0,2,1)
F <- matrix(0,2,1)
J <- matrix(0,2,2)
B[1,1] <- b1
B[2,1] <- b2
for (i in 1:max_iter){
  F[1,1] <- (b1^2 * (1-(b1*t3))) + (b2^2*(pt2-(b2*pt3))) + b1*b2*((2*c1)-(3*b1*c2)-(3*b2*c3)) - 1
  F[2,1] <- ((t4*b1^4)+(pt4*b2^4)) + (4*(b1^3)*b2*c4) + (4*(b2^3)*b1*c5) + (6*(b1^2)*(b2^2)*c6) - 3
  b11 <- b1+0.00001
  b21 <- b2+0.00001
  J[1,1] <- ((b11^2 * (1-(b11*t3))) + (b2^2*(pt2-(b2*pt3))) + b11*b2*((2*c1)-(3*b11*c2)-(3*b2*c3)) - 1 - F[1,1]) / 0.00001
  J[1,2] <- ((b1^2 * (1-(b1*t3))) + (b21^2*(pt2-(b21*pt3))) + b1*b21*((2*c1)-(3*b1*c2)-(3*b21*c3)) - 1 - F[1,1]) / 0.00001
  J[2,1] <- (((t4*b11^4)+(pt4*b2^4)) + (4*(b11^3)*b2*c4) + (4*(b2^3)*b11*c5) + (6*(b11^2)*(b2^2)*c6) - 3 - F[2,1]) / 0.00001
  J[2,2] <- (((t4*b1^4)+(pt4*b21^4)) + (4*(b1^3)*b21*c4) + (4*(b21^3)*b1*c5) + (6*(b1^2)*(b21^2)*c6) - 3 - F[2,1]) / 0.00001
  B <- B - (solve(J) %*% F)
  b1 <- B[1,1]
  b2 <- B[2,1]
}
a  <- -b2
  
  return(list(md = md, sx = sx, a = a, b1 = b1, b2 = b2))
}

