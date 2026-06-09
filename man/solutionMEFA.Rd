\name{solutionMEFA}
\alias{solutionMEFA}
\alias{solutionMCFA}
\title{High-Level MORGANA Workflows for EFA and CFA}
\description{
Runs end-to-end MORGANA analyses with correlated residuals.

\code{solutionMEFA()} performs \strong{exploratory} analysis with automatic
EREC-based doublet detection. It orchestrates the full pipeline:
\enumerate{
  \item \code{\link{LargestDoubletsEREC}} -- identifies the largest doublets;
  \item \code{\link{MorganaEFA}} -- runs the three-stage ULS extraction;
  \item \code{\link{LoseferMorganaH0}} -- estimates the H0 distribution;
  \item \code{\link{LoseferMorgana}} -- computes RMSEA, CFI, and NNFI/TL;
  \item \code{\link{ParallelAnalysisEREC}} -- (optional) computes PA reference
    thresholds for EREC significance assessment.
}

\code{solutionMCFA()} performs \strong{confirmatory} analysis using
user-provided doublet pairs, following the mild-restricted CFA approach
described in Lorenzo-Seva (2025). The same steps are run except doublets are
pre-specified and PA is not computed.

Both functions require a positive-definite correlation/covariance matrix; see
Lorenzo-Seva & Ferrando (2021) for handling non-positive-definite cases.
}
\usage{
solutionMEFA(
  database,
  r,
  cor = 2,
  conver = 1e-04,
  maxiter = 500,
  display = 1,
  montecarlo = 5000,
  pa_reps = 500
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
  \item{database}{Numeric matrix or data frame of observed scores
    (rows = subjects, columns = items/variables).}
  \item{r}{Number of common factors to extract.}
  \item{cor}{Type of dispersion matrix:
    \code{1} = covariance matrix,
    \code{2} = Pearson correlation (default),
    \code{3} = polychoric correlation (requires \pkg{psych}).}
  \item{conver}{Convergence criterion passed to iterative routines
    (default \code{1e-4}).}
  \item{maxiter}{Maximum number of iterations in \code{Morgana()}
    (default \code{500}).}
  \item{display}{If \code{1} (default), prints fit summary from
    \code{LoseferMorgana()}; if \code{0}, suppresses printing.}
  \item{montecarlo}{Number of Monte Carlo replications for \code{LoseferMorganaH0()}
    (default \code{5000}). Reduce for faster exploratory runs.}
  \item{pa_reps}{Number of Monte Carlo replications for Parallel Analysis of
    EREC thresholds via \code{\link{ParallelAnalysisEREC}} (default \code{500}).
    Set to \code{0} to skip Parallel Analysis and omit \code{$pa} from the
    output. Computationally intensive for large item sets.}
  \item{doublets}{(\code{solutionMCFA} only) Two-column integer matrix or
    data frame whose rows are the item-index pairs to evaluate as doublets.}
}
\value{
A named list. Elements common to both functions:
\itemize{
  \item \code{loadings}: estimated loading matrix (m x r), sign-corrected by
    \code{\link{posit}}.
  \item \code{residual_covariance}: residual covariance matrix (m x m).
  \item \code{residual_correlation}: residual correlation matrix (m x m);
    off-diagonal elements are non-zero only at the freed doublet positions.
  \item \code{doublets}: integer/numeric matrix (nd x 3) with columns
    \code{P1}, \code{P2} (item indices), and \code{EREC} (estimated residual
    correlation for that pair).
  \item \code{fit}: list returned by \code{\link{LoseferMorgana}} containing
    \code{Um} (adjusted chi-square), \code{d} (df), \code{d0} (null df),
    \code{C0} (null chi-square), \code{RMSEA}, \code{CFI}, \code{TL}.
  \item \code{h0}: transformation parameters from \code{\link{LoseferMorganaH0}}.
}

\code{solutionMEFA} additionally returns:
\itemize{
  \item \code{pa}: list from \code{\link{ParallelAnalysisEREC}}, or \code{NULL}
    when \code{pa_reps = 0}. Contains \code{doublets} (observed EREC values),
    \code{pa_thresholds} (PA cutoff per rank position), and \code{significant}
    (logical vector flagging doublets that exceed the PA threshold).
}
}
\details{
The EREC-based doublet detection and the MORGANA estimation procedure are
described in Ferrando et al. (2022). The mild-restricted CFA workflow (used
by \code{solutionMCFA}) is described in Lorenzo-Seva (2025). LOSEFER fit
correction is described in Lorenzo-Seva & Ferrando (2023).

The same workflows are implemented in the program \strong{FACTOR}
(\url{https://psico.fcep.urv.cat/utilitats/factor/}). This R implementation is
independent and has not been formally validated against FACTOR; numerical
differences may exist.
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

Lorenzo-Seva, U., & Ferrando, P. J. (2021). Not positive definite correlation
matrices in exploratory item factor analysis: Causes, consequences and a
proposed solution.
\emph{Structural Equation Modeling: A Multidisciplinary Journal}, 28(1), 138--147.
\doi{10.1080/10705511.2020.1735393}
}
\examples{
\donttest{
library(morganafa)

db <- as.matrix(read.table(system.file("extdata", "database.dat",
                                        package = "morganafa")))

# EFA com Parallel Analysis (pa_reps reduzido para o exemplo)
res_efa <- solutionMEFA(db[1:250, ], r = 2, cor = 2,
                        display = 0, montecarlo = 60, pa_reps = 50)
res_efa$doublets
res_efa$pa$significant
res_efa$fit$CFI

# CFA com pares de doublets fornecidos pelo usuario
pairs <- matrix(c(1, 2,
                  6, 7), ncol = 2, byrow = TRUE)
res_cfa <- solutionMCFA(db[1:250, ], r = 2, doublets = pairs,
                        cor = 2, display = 0, montecarlo = 60)
res_cfa$fit$RMSEA
}
}
\seealso{
\code{\link{LargestDoubletsEREC}},
\code{\link{ParallelAnalysisEREC}},
\code{\link{MorganaEFA}},
\code{\link{MorganaCFA}},
\code{\link{LoseferMorgana}},
\code{\link{LoseferMorganaH0}}
}
