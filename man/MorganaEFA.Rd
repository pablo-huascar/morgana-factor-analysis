\name{MorganaEFA}
\alias{MorganaEFA}
\alias{MorganaCFA}
\alias{Morgana}
\title{Core MORGANA Estimation Routines}
\description{
Core routines for MORGANA factor extraction with correlated residuals.
\code{MorganaEFA()} performs exploratory extraction with automatic doublet selection.
\code{MorganaCFA()} evaluates user-specified pairs using EREC before extraction.
\code{Morgana()} executes the core unrestricted estimation step given free doublets.
}
\usage{
MorganaEFA(X, r, conver, maxiter, cor)

MorganaCFA(X, r, conver, maxiter, cor, doblets)

Morgana(X, r, conver, maxiter, corr, doblets)
}
\arguments{
  \item{X}{Numeric matrix/data frame of observed scores.}
  \item{r}{Number of factors.}
  \item{conver}{Convergence criterion.}
  \item{maxiter}{Maximum number of iterations.}
  \item{cor}{Dispersion type (1 covariance, 2 Pearson, 3 polychoric).}
  \item{corr}{Same as \code{cor}; parameter name used in \code{Morgana()}.}
  \item{doblets}{Matrix of item-pair indices (and EREC values when passed to \code{Morgana()}).}
}
\value{
\code{MorganaEFA()} and \code{MorganaCFA()} return a list with:
\itemize{
  \item \code{doublets}: selected/evaluated doublets,
  \item \code{A}: loading matrix,
  \item \code{RESCOV}: residual covariance matrix,
  \item \code{RESCOR}: residual correlation matrix.
}
\code{Morgana()} returns \code{L}, \code{RESCOV}, and \code{RESCOR}.
}
\examples{
\donttest{
db <- as.matrix(read.table(system.file("extdata", "database.dat", package = "morganafa")))
out <- MorganaEFA(db[1:200, ], r = 2, conver = 1e-4, maxiter = 100, cor = 2)
dim(out$A)
}
}
\seealso{
\code{\link{LargestDoubletsEREC}}, \code{\link{DoubletsEREC}}
}
