\name{LargestDoubletsEREC}
\alias{LargestDoubletsEREC}
\alias{DoubletsEREC}
\title{EREC-Based Detection and Evaluation of Doublets}
\description{
Functions for identifying correlated-residual pairs (doublets) using the EREC index.
\code{LargestDoubletsEREC()} ranks all candidate pairs and returns the largest values under identifiability constraints.
\code{DoubletsEREC()} evaluates EREC only for user-provided pairs.
}
\usage{
LargestDoubletsEREC(X, r, cor)

DoubletsEREC(X, r, cor, PairsDoublets)
}
\arguments{
  \item{X}{Numeric matrix/data frame of observed scores.}
  \item{r}{Number of factors.}
  \item{cor}{Dispersion type (1 covariance, 2 Pearson, 3 polychoric).}
  \item{PairsDoublets}{Two-column matrix/data frame with candidate pairs to evaluate.}
}
\value{
A list with:
\itemize{
  \item \code{doublets}: matrix with columns pair-1, pair-2, and EREC value,
  \item \code{nd}: number of returned pairs.
}
}
\examples{
\donttest{
db <- as.matrix(read.table(system.file("extdata", "database.dat", package = "morganafa")))
LargestDoubletsEREC(db[1:200, ], r = 2, cor = 2)$doublets

pairs <- matrix(c(1, 2,
                  6, 7), ncol = 2, byrow = TRUE)
DoubletsEREC(db[1:200, ], r = 2, cor = 2, PairsDoublets = pairs)$doublets
}
}
