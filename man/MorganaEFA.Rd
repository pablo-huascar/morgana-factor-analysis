\name{MorganaEFA}
\alias{MorganaEFA}
\alias{MorganaCFA}
\alias{Morgana}
\title{Core MORGANA Estimation Routines}
\description{
Core routines for MORGANA factor extraction with correlated residuals, as
described in Ferrando et al. (2022).

\code{MorganaEFA()} is the exploratory entry point: it calls
\code{\link{LargestDoubletsEREC}} to identify the most substantial doublets
and then runs \code{Morgana()} to produce the final unrestricted solution.

\code{MorganaCFA()} is the confirmatory entry point used by the mild-restricted
CFA workflow (Lorenzo-Seva 2025): it calls \code{\link{DoubletsEREC}} to
compute EREC values for user-specified pairs and then runs \code{Morgana()}.

\code{Morgana()} is the central estimation function. It implements the
two-stage Mulaik-style approach:
\enumerate{
  \item Initialise correlated-residual estimates from the EREC values and the
    partial-correlation matrix (Guttman anti-image approach, Eq. 4 of
    Ferrando et al. 2022).
  \item Iteratively refine the loading matrix on the reduced correlation matrix
    \eqn{R_\text{red} = R - U \hat\Psi_r U}, where U and \eqn{\hat\Psi_r}
    are updated at each step, until convergence.
  \item Compute the final residual covariance (\code{RESCOV}) and residual
    correlation (\code{RESCOR}) matrices.
}
}
\usage{
MorganaEFA(X, r, conver, maxiter, cor)

MorganaCFA(X, r, conver, maxiter, cor, doblets)

Morgana(X, r, conver, maxiter, corr, doblets)
}
\arguments{
  \item{X}{Numeric matrix or data frame of observed scores
    (rows = subjects, columns = items).}
  \item{r}{Number of common factors.}
  \item{conver}{Convergence criterion for the iterative loading refinement
    loop.}
  \item{maxiter}{Maximum number of iterations in \code{Morgana()}.}
  \item{cor}{Dispersion type: \code{1} covariance, \code{2} Pearson,
    \code{3} polychoric.}
  \item{corr}{Same as \code{cor}; the parameter name used in \code{Morgana()}
    to avoid masking base R's \code{cor()}.}
  \item{doblets}{Matrix of item-pair doublets (nd x 2 or nd x 3). When passed
    to \code{Morgana()}, column 3 contains the EREC value used to initialise
    the partial-correlation adjustment.}
}
\value{
\code{MorganaEFA()} and \code{MorganaCFA()} return a list with:
\itemize{
  \item \code{doublets}: selected/evaluated doublets (nd x 3: P1, P2, EREC).
  \item \code{A}: loading matrix (m x r).
  \item \code{RESCOV}: residual covariance matrix (m x m).
  \item \code{RESCOR}: residual correlation matrix (m x m).
}

\code{Morgana()} returns a list with \code{L} (loading matrix),
\code{RESCOV}, and \code{RESCOR}.
}
\details{
The choice of ULS (MINRES) as the estimation criterion follows the rationale
in Ferrando et al. (2022, p. 632): ULS is computationally robust and
defensible with tetrachoric/polychoric correlation matrices.

The same routines are implemented in the program \strong{FACTOR}
(\url{https://psico.fcep.urv.cat/utilitats/factor/}). This R implementation is
independent and has not been formally validated against FACTOR.
}
\references{
Ferrando, P. J., Hernandez-Dorado, A., & Lorenzo-Seva, U. (2022).
Detecting correlated residuals in exploratory factor analysis: New proposals
and a comparison of procedures.
\emph{Structural Equation Modeling: A Multidisciplinary Journal}, 29(4), 630--638.
\doi{10.1080/10705511.2021.2004543}

Lorenzo-Seva, U. (2025). Mild-restricted confirmatory factor analysis: a new
approach to assessing population factor models.
\emph{Communications in Statistics -- Simulation and Computation}.
\doi{10.1080/03610918.2025.2453473}
}
\examples{
\donttest{
db <- as.matrix(read.table(system.file("extdata", "database.dat",
                                        package = "morganafa")))
out <- MorganaEFA(db[1:200, ], r = 2, conver = 1e-4, maxiter = 100, cor = 2)
dim(out$A)          # loading matrix
out$doublets        # doublets with EREC values
}
}
\seealso{
\code{\link{LargestDoubletsEREC}},
\code{\link{DoubletsEREC}},
\code{\link{solutionMEFA}},
\code{\link{solutionMCFA}}
}
