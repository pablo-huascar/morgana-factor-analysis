
posit <- function(A){
  n <- nrow(A)
  m <- ncol(A)  
  tmp <- zeros(n,m)
  id <- ones(n,1)
  
  for(i in 1: m){

    
    lugmin <- which.min(A[,i])
    lugmax <- which.max(A[,i])
    
    a <- abs(min(A[,i])) # min 
    b <- abs(max(A[,i])) # max
    
    
    if(a<b & lugmin<lugmax){
      tmp[,i] <- A[,i] * (-id)
    }else{
      tmp[,i] <- A[,i]
    }
    
    if(all(A[,i] < 0)){
      tmp[,i] <- A[,i] * (-id)
    }else{
      tmp[,i] <- A[,i]
    }
    

  }
  return(tmp)
}

