\name{morganafa-package}
\alias{morganafa-package}
\alias{morganafa}
\title{MORGANA Factor Analysis with Correlated Residuals}
\description{
\code{morganafa} implements the MORGANA (\emph{Minimum expected bias Of
Residual and loadinG vAlues iN sAmple estimates}) methodology for exploratory
and confirmatory factor analysis with correlated residuals (doublets).

The statistical methodology is described in:

\itemize{
  \item Ferrando, P. J., Hernandez-Dorado, A., & Lorenzo-Seva, U. (2022).
    Detecting correlated residuals in exploratory factor analysis: New proposals
    and a comparison of procedures.
    \emph{Structural Equation Modeling: A Multidisciplinary Journal},
    29(4), 630--638. \doi{10.1080/10705511.2021.2004543}

  \item Lorenzo-Seva, U. (2025). Mild-restricted confirmatory factor analysis:
    a new approach to assessing population factor models.
    \emph{Communications in Statistics -- Simulation and Computation}.
    \doi{10.1080/03610918.2025.2453473}
}

The same methods are also available in the program \strong{FACTOR}
(\url{https://psico.fcep.urv.cat/utilitats/factor/}), developed and maintained
by the authors of the articles at the Universitat Rovira i Virgili.

\strong{Important disclaimer:} this R package was developed independently of
the FACTOR program. No formal comparison between the numerical outputs of the
two implementations has been performed, and results may differ due to
implementation choices (convergence tolerances, random-number generation,
operation order). Users are advised to validate results carefully, and to
report which implementation was used when publishing.
}
\details{
Main capabilities:
\itemize{
  \item EREC-based doublet detection with Ledermann-bound restriction:
    \code{\link{LargestDoubletsEREC}}, \code{\link{DoubletsEREC}}.
  \item Parallel Analysis reference thresholds for EREC:
    \code{\link{ParallelAnalysisEREC}}.
  \item Three-stage ULS unrestricted estimation with correlated residuals:
    \code{\link{Morgana}}, \code{\link{MorganaEFA}}, \code{\link{MorganaCFA}}.
  \item LOSEFER empirical correction for fit indices (RMSEA, CFI, NNFI/TL):
    \code{\link{LoseferMorgana}}, \code{\link{LoseferMorganaH0}}.
  \item High-level workflows:
    \code{\link{solutionMEFA}} (exploratory) and \code{\link{solutionMCFA}} (confirmatory).
}
}
\author{MORGANA Developers}
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
\seealso{
\code{\link{solutionMEFA}},
\code{\link{solutionMCFA}},
\code{\link{LargestDoubletsEREC}},
\code{\link{ParallelAnalysisEREC}},
\code{\link{MorganaEFA}},
\code{\link{MorganaCFA}}
}
\keyword{package}
