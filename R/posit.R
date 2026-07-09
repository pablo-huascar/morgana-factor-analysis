#' Reflect factor loading signs
#'
#' Reflects (multiplies by -1) every column of a loading matrix whose sum of
#' loadings is negative, so that factors point in the "positive" direction.
#' The sign of a factor is arbitrary in factor analysis; this only affects
#' presentation, never fit.
#'
#' @param A A numeric loading matrix (items in rows, factors in columns).
#'
#' @return A matrix with the same dimensions as `A`, with reflected columns
#'   where needed.
#'
#' @export
#' @examples
#' A <- cbind(c(-.7, -.8, -.6), c(.5, .6, .4))
#' posit(A)
posit <- function(A) {
  A <- as.matrix(A)
  signs <- ifelse(colSums(A) < 0, -1, 1)
  sweep(A, 2, signs, `*`)
}
