
#############################
###                       ###
###      PARTCOR          ### 
###                       ###
#############################




PartCor <- function(R) {
  # Compute partial correlation
  m <- nrow(R)
  ir <- as.matrix(inverse(R)$result)
  D <- diag(diag(ir))
  Dd <- diag(diag(D^(-1/2)))
  
  AI <- Dd %*% ir %*% Dd
  # AI es la matriz Anti-Imagen de Guttman, y Ra es la matriz de correlaciones parciales
  Ra <- 2 * diag(m) - AI
  
  Ra <- Nanas(Ra)
  return(Ra)
}

