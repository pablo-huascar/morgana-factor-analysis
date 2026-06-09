\name{LoseferMorgana}
\alias{LoseferMorgana}
\alias{LoseferMorganaH0}
\title{LOSEFER Empirical Fit-Index Correction for MORGANA Solutions}
\description{
LOSEFER (\emph{LOw-Sampling-Error Fit Empirical-Reference}) is a simulation-based
procedure that corrects the chi-square test statistic so that its distribution
approximates the theoretical chi-square distribution more closely under
small-to-moderate samples. The adjusted statistic is then used to compute
RMSEA, CFI, and NNFI (Tucker-Lewis index).

\code{LoseferMorganaH0()} estimates the empirical null distribution of the
fit function value by generating Monte Carlo samples from the model reproduced
under H0 (i.e. with loadings and freed residual correlations as population
values but no additional misspecification). It returns parameters for a
nonlinear transformation -- a combination of a cubic-root normalisation and a
Fleishman polynomial -- that maps the simulated distribution to the theoretical
chi-square shape.

\code{LoseferMorgana()} applies those transformation parameters to the observed
fit function value and computes RMSEA, CFI, and NNFI/TL from the adjusted
chi-square statistic.

These functions implement the LOSEFER methodology described in
Lorenzo-Seva & Ferrando (2023) and are called internally by
\code{\link{solutionMEFA}} and \code{\link{solutionMCFA}}.
}
\usage{
LoseferMorgana(
  S,
  A,
  N,
  RESCOR,
  nd,
  mchi,
  sxchi,
  a_n,
  b1_n,
  b2_n,
  a_f,
  b_f,
  c_f,
  d_f,
  display
)

LoseferMorganaH0(A, RESCOR, doblets, N, montecarlo)
}
\arguments{
  \item{S}{Observed covariance or correlation matrix (m x m).}
  \item{A}{Estimated loading matrix (m x r).}
  \item{N}{Sample size.}
  \item{RESCOR}{Residual correlation matrix (m x m), with non-zero
    off-diagonal elements only at freed doublet positions.}
  \item{nd}{Number of freely estimated residual pairs (doublets).}
  \item{mchi, sxchi}{Mean and standard deviation of the cubic-root-transformed
    H0 distribution, returned by \code{LoseferMorganaH0} as \code{v} and
    \code{vv} respectively.}
  \item{a_n, b1_n, b2_n}{Quadratic normalisation parameters from
    \code{LoseferMorganaH0}.}
  \item{a_f, b_f, c_f, d_f}{Fleishman polynomial parameters from
    \code{LoseferMorganaH0}.}
  \item{display}{If \code{1}, prints a fit summary table; if \code{0}, silent.}
  \item{doblets}{Two-column integer matrix of freed doublet pairs passed to
    \code{fufa1} during H0 simulation.}
  \item{montecarlo}{Number of Monte Carlo replications for H0 estimation
    (default \code{5000} when called via \code{solutionMEFA}).}
}
\value{
\code{LoseferMorgana()} returns a named list:
\itemize{
  \item \code{Um}: adjusted chi-square statistic.
  \item \code{d}: degrees of freedom of the current model.
  \item \code{d0}: degrees of freedom of the null (independence) model.
  \item \code{C0}: chi-square of the null model.
  \item \code{RMSEA}: Root Mean Square Error of Approximation.
  \item \code{CFI}: Comparative Fit Index.
  \item \code{TL}: Tucker-Lewis (NNFI) index.
}

\code{LoseferMorganaH0()} returns a named list of transformation parameters:
\code{v}, \code{vv}, \code{a_n}, \code{b1_n}, \code{b2_n}, \code{a_f},
\code{b_f}, \code{c_f}, \code{d_f}.
}
\details{
LOSEFER is based on a bootstrapped approach analogous to Bollen & Stine (1992)
but uses nonlinear transformations to match the four moments of the theoretical
chi-square distribution, making it more accurate than the plain bootstrap.

Typical threshold values: RMSEA < 0.07, CFI > 0.95, chi-square/df < 2
(see Lorenzo-Seva 2025, Section 4.4).

The same procedure is implemented in the program \strong{FACTOR}
(\url{https://psico.fcep.urv.cat/utilitats/factor/}). This R implementation is
independent and has not been formally validated against FACTOR.
}
\references{
Lorenzo-Seva, U., & Ferrando, P. J. (2023). A simulation-based scaled test
statistic for assessing model-data fit in least-squares unrestricted
factor-analysis solutions.
\emph{Methodology}, 19(2), 96--115.
\doi{10.5964/meth.9839}

Lorenzo-Seva, U. (2025). Mild-restricted confirmatory factor analysis: a new
approach to assessing population factor models.
\emph{Communications in Statistics -- Simulation and Computation}.
\doi{10.1080/03610918.2025.2453473}

Bollen, K. A., & Stine, R. A. (1992). Bootstrapping goodness-of-fit measures
in structural equation models.
\emph{Sociological Methods & Research}, 21(2), 205--229.
\doi{10.1177/0049124192021002004}
}
\seealso{
\code{\link{solutionMEFA}}, \code{\link{solutionMCFA}}
}
