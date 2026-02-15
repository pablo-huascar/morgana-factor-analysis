



colmin <- function(matrix){
  mat <- c()
  for(i in 1: ncol(matrix)){
    mat[i] <- min(matrix[,i])
  }
  return(mat)
}
