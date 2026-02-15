
##########################
###                    ###
###      BERGE         ### 
###                    ###
##########################


## A general solution to Mosier's oblique procrustes problem
## Ten Berge and Nevels (1977). Psychometrika

# f -> Factor solution to be rotated
# PHI -> Target matrix
# t -> transformation matrix

###############################################################################





BERGE <- function(f, PHI){
  
  n <- dim(f)[1]
  m <- dim(f)[2]
  n2 <- dim(PHI)[1]
  m2 <- dim(PHI)[2]
  
  phi <- PHI[, 1]
  
  ########################################################
  # First step in the paper
  ########################################################
  
  # EigenVector descomposition of A'A = UCU'
  
  tmp0 <- svd(f)
  
  V <- tmp0$u
  C <- diag(tmp0$d)
  U <- tmp0$v
  
  # Evaluation of q
  
  q <- 0
  
  for(i in 1:m){
    if(C[i,i]==C[m,m]){
      q <- q+1
    }
  }

  ########################################################
  # Second step in the paper
  ########################################################
  
  # compute x = U'F'phi
  
  x <- t(U)%*%t(f)%*%phi
  
  
  s <- ssq(x[(m-q+1):m])
  
  ########################################################
  # Third Step in the paper 
  ########################################################
  
  
  # s > 0
  # case 1
  # phim'(0)
  
  
  if(s>0){
    tmp <- 0
    
    for(i in 1:m){
      tmp <- tmp + ( (x[i]*x[i]) / (C[i,i]%*%C[i,i]))
    }
    
    phim <- 1- tmp
    
    # b1 <- 38 or 39
    
    if(m>1){
      C_x <- min(diag(C[1:m-q, 1:m-q])) - eye(m-q) * abs(x[1:m-q])
      bas <- C[m,m]-(ssq(x[m-q+1:m])^0.5)
    }
    else{
      C_x <- C-abs(x)
      bas <- C-abs(x)
    }
    
    if(phim <= 0){
      start <- c(0, C_x, bas)
    }
    else{start <-c(C_x, bas)}
    
    
    b1 <- min(start)
    
    # bo = (33)
    
    bo <- -(t(x)%*%solve(C)%*%x)
    
    dis <- b1 + abs(bo)
    
    # iteration Phil
  
    l <- m
    p <- matrix()
    
    for(i in 1:20){
      D <- (C[1:l, 1:l] - b1 %*% eye(l))^2
      fb <- 1- sum(t(x)%*%solve(D)%*%x)
      
      if(fb<0){
        b1 <- b1 - dis
      }
      else{
        b1 <- b1 + dis
      }
      
      dis <- dis /2
      p <- c(p, fb)
    }
    
    # Wo <- (C - boI)^1 *x
    
    D <- C - eye(m)%*% b1
    Wo <- solve(D) %*% x
    
  }
  else{
    tmp <- 0
    for(i in 1:(m-q)){
      tmp <- tmp + ( (x[i]*x[i]) / (C[i,i] - C[m,m]))
    }
    phim_q <- 1-tmp
    
    # CASE 3
    
    if(phim_q >= 0){
      # Wo =13
      Wo <- zeros(m, 1)
      tmp <- solve(C[1:m-q, 1:m-q] - C[m,m] %*% eye(m-q)) %*% x[1:m-q]
      l <- m-q
      D <- (C[1:l, 1:l] - C[m,m] %*% eye(l))^2
      fb <- (1 - sum(t(x[1:l])%*% solve(D)%*%x[1:l]))^0.5
      
      for(i in 1:m-q){
        Wo[i] <- x[i]/(C[i,i] - C[m,m])
      }
      Wo[m-q+1] <- fb
    } else{
      # CASE 2
      b1 <- C[m,m]
      # bo <- 33
      bo <- -(t(x)%*%pinv(C)%*%x)
      dis <- b1 + abs(bo)
      
      # iteration Phil (b)
      
      
      l <- m-q
      p <- c()
      
      for(i in 1:20){
        D <- (C[1:l, 1:l]- b1 %*%x[1:l])
        fb <- 1-sum(t(x[1:l])%*%pinv(D)%*%x[1:l])
        if(fb<0){
          b1 <- b1-dis
        } else{
          b1 <- b1+dis
        }
        dis <- dis/2
        p <- c(p, fb)
      }
      D <- C- eye(m)%*%b1
      Wo <- solve(D)%*%x
    }
    
  
  }
    t0 <- U%*%Wo
    t_mat <- matrix(0, nrow = nrow(f), ncol = 1)
    t_mat[,1] <- t0
    return(t_mat)



}
