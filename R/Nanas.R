

Nanas <- function(matriz) {
  # Obtener las posiciones de los NaN
  posiciones_nan <- which(is.nan(matriz), arr.ind = TRUE)
  
  # Reemplazar NaN por 0 en esas posiciones
  matriz[posiciones_nan] <- 0
  
  # Poner la diagonal a 1
  diag(matriz) <- 1
  
  return(matriz)
}

