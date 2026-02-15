\name{solutionMEFA}
\alias{solutionMEFA}
\alias{solutionMCFA}
\title{High-Level MORGANA Workflows for EFA and CFA}
\description{
Runs end-to-end MORGANA analyses with correlated residuals.
\code{solutionMEFA()} performs exploratory analysis with automatic EREC-based doublet detection.
\code{solutionMCFA()} performs confirmatory analysis using user-provided doublet pairs.
}
\usage{
solutionMEFA(
  database,
  r,
  cor = 2,
  conver = 1e-04,
  maxiter = 500,
  display = 1,
  montecarlo = 5000
)

solutionMCFA(
  database,
  r,
  doublets,
  cor = 2,
  conver = 1e-04,
  maxiter = 500,
  display = 1,
  montecarlo = 5000
)
}
\arguments{
  \item{database}{Numeric matrix or data frame of observed scores (rows are subjects, columns are items).}
  \item{r}{Number of common factors.}
  \item{cor}{Dispersion type: \code{1} covariance, \code{2} Pearson correlation, \code{3} polychoric correlation.}
  \item{conver}{Convergence criterion used in iterative routines.}
  \item{maxiter}{Maximum number of iterations.}
  \item{display}{If \code{1}, prints fit summary from LOSEFER; if \code{0}, suppresses printing.}
  \item{montecarlo}{Number of Monte Carlo replications used in \code{LoseferMorganaH0}.}
  \item{doublets}{Two-column matrix/data frame of item index pairs to evaluate in confirmatory mode (only for \code{solutionMCFA}).}
}
\value{
A list with:
\itemize{
  \item \code{loadings}: estimated loading matrix,
  \item \code{residual_covariance}: residual covariance matrix,
  \item \code{residual_correlation}: residual correlation matrix (selected/free residuals),
  \item \code{doublets}: matrix with item pairs and EREC values,
  \item \code{fit}: fit statistics from \code{LoseferMorgana},
  \item \code{h0}: transformation parameters estimated by \code{LoseferMorganaH0}.
}
}
\examples{
\donttest{
library(morganafa)

db <- as.matrix(read.table(system.file("extdata", "database.dat", package = "morganafa")))

res_efa <- solutionMEFA(db[1:250, ], r = 2, cor = 2, display = 0, montecarlo = 60)
res_efa$fit$CFI

pairs <- matrix(c(1, 2,
                  6, 7), ncol = 2, byrow = TRUE)
res_cfa <- solutionMCFA(db[1:250, ], r = 2, doublets = pairs, cor = 2, display = 0, montecarlo = 60)
res_cfa$fit$RMSEA
}
}
\seealso{
\code{\link{MorganaEFA}}, \code{\link{MorganaCFA}},
\code{\link{LoseferMorgana}}, \code{\link{LoseferMorganaH0}}
}
