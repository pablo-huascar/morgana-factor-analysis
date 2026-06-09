\name{LargestDoubletsEREC}
\alias{LargestDoubletsEREC}
\alias{DoubletsEREC}
\title{EREC-Based Detection and Evaluation of Doublets}
\description{
Functions for identifying correlated-residual pairs (doublets) using the
\strong{EREC} (\emph{Expected REsidual Correlation direct Change}) index,
also called EPC1 in Ferrando et al. (2022, Equation 5).

\code{LargestDoubletsEREC()} scans all m*(m-1)/2 candidate pairs, estimates
the EREC for each, and returns the \code{nd} pairs with the largest absolute
EREC values. The number of pairs is bounded by the Ledermann (1937) constraint
\eqn{n_d \le g - r}, where \emph{g} is the Ledermann bound on the total number
of factors (common + minor) and \emph{r} is the number of common factors.
This restriction ensures the extended model remains identifiable.

\code{DoubletsEREC()} evaluates EREC only for the pairs explicitly provided
by the user (confirmatory use).

For each candidate pair (j, k), both functions implement the three-stage ULS
estimation of Ferrando et al. (2022):
\enumerate{
  \item Fit MINRES to the \emph{core set} of m-2 variables (omitting j and k).
  \item Estimate loadings of j and k by \emph{extension analysis}
    (Eq. 7: \eqn{\hat\lambda_j = (\hat\Lambda_\text{core}'\hat\Lambda_\text{core})^{-1}
    \hat\Lambda_\text{core}' r_{j,\text{core}}}).
  \item Estimate the residual correlation for the pair (Eq. 8:
    \eqn{\hat{R}_{uu} = \hat\Psi^{-1}(R - \hat\Lambda\hat\Lambda')\hat\Psi^{-1}}).
}
The off-diagonal element of \eqn{\hat{R}_{uu}} at position (j, k) is the EREC
for that pair. An IPSI correction is applied across all pairs to standardise
the estimates to a correlation scale before ranking.

To assess whether an observed EREC value exceeds chance, use
\code{\link{ParallelAnalysisEREC}}, which provides Monte Carlo reference values
as described in the section \emph{"Reference values and upper bounds"} of
Ferrando et al. (2022).
}
\usage{
LargestDoubletsEREC(X, r, cor)

DoubletsEREC(X, r, cor, PairsDoublets)
}
\arguments{
  \item{X}{Numeric matrix or data frame of observed scores
    (rows = subjects, columns = items).}
  \item{r}{Number of common factors.}
  \item{cor}{Type of dispersion matrix:
    \code{1} = covariance, \code{2} = Pearson correlation,
    \code{3} = polychoric correlation.}
  \item{PairsDoublets}{(\code{DoubletsEREC} only) Two-column integer matrix
    or data frame specifying the candidate pairs to evaluate.}
}
\value{
A list with:
\itemize{
  \item \code{doublets}: numeric matrix (nd x 3) with columns for the two
    item indices (P1, P2) and the estimated EREC value.
  \item \code{nd}: integer, number of doublet positions evaluated.
}
}
\details{
The estimation procedure is described as MORGANA in Ferrando et al. (2022)
and is the recommended method for doublet detection, outperforming fitted
residuals and partial correlations in both sensitivity and specificity across
simulation conditions (see Table 6 in that article).

The same procedure is implemented in the program \strong{FACTOR}
(\url{https://psico.fcep.urv.cat/utilitats/factor/}). This R implementation is
independent and has not been formally validated against FACTOR.
}
\references{
Ferrando, P. J., Hernandez-Dorado, A., & Lorenzo-Seva, U. (2022).
Detecting correlated residuals in exploratory factor analysis: New proposals
and a comparison of procedures.
\emph{Structural Equation Modeling: A Multidisciplinary Journal}, 29(4), 630--638.
\doi{10.1080/10705511.2021.2004543}

Ledermann, W. (1937). On the rank of the reduced correlational matrix in
multiple-factor analysis. \emph{Psychometrika}, 2, 85--93.
\doi{10.1007/BF02288062}
}
\examples{
\donttest{
db <- as.matrix(read.table(system.file("extdata", "database.dat",
                                        package = "morganafa")))
# Maiores doublets por EREC (EFA exploratorio)
LargestDoubletsEREC(db[1:200, ], r = 2, cor = 2)$doublets

# EREC para pares fornecidos pelo usuario (uso confirmatorio)
pairs <- matrix(c(1, 2,
                  6, 7), ncol = 2, byrow = TRUE)
DoubletsEREC(db[1:200, ], r = 2, cor = 2, PairsDoublets = pairs)$doublets
}
}
\seealso{
\code{\link{ParallelAnalysisEREC}},
\code{\link{MorganaEFA}},
\code{\link{solutionMEFA}}
}
