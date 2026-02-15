

# Cargar el paquete MASS para usar la funci?n ginv


# Funci?n que calcula la inversa de una matriz y, si es singular, calcula la pseudo-inversa
inverse <- function(mat) {
  # Intentar calcular la inversa de la matriz
  inv_mat <- tryCatch(solve(mat), error = function(e) NULL)
  
  if (is.null(inv_mat)) {
    # Si la matriz es singular y no se pudo calcular la inversa, calcular la pseudo-inversa
    pseudo_inv_mat <- ginv(mat)
    return(list(result = pseudo_inv_mat, method = "pseudo-inverse"))
  } else {
    # Si se pudo calcular la inversa, devolverla
    return(list(result = inv_mat, method = "inverse"))
  }
}


