\name{posit}
\alias{posit}
\title{Orient Factor Loadings to Positive Direction}
\description{
Applies a sign-orientation rule to columns of a loading matrix to improve interpretability,
flipping columns when needed so dominant loadings follow a positive direction.
}
\usage{
posit(A)
}
\arguments{
  \item{A}{Numeric loading matrix.}
}
\value{
A numeric matrix with the same dimensions as \code{A}, potentially with sign-flipped columns.
}
\examples{
A <- matrix(c(-0.6, -0.4,
               0.5,  0.1,
               0.7, -0.2), ncol = 2, byrow = TRUE)
posit(A)
}
