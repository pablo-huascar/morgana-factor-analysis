\name{ParallelAnalysisEREC}
\alias{ParallelAnalysisEREC}
\title{Parallel Analysis Reference Thresholds for the EREC Index}
\description{
Computes Monte Carlo reference (cutoff) values for the EREC doublet-detection
index, following the Parallel Analysis procedure described in the section
\emph{"Reference values and upper bounds"} of Ferrando et al. (2022).

The rationale is that observed EREC values should be compared not to an
arbitrary fixed cutoff, but to the distribution of EREC values that would
arise by sampling error alone (i.e. when no doublets exist in the population).
\code{ParallelAnalysisEREC} generates this null distribution, estimates
rank-specific thresholds, and returns a significance flag for each doublet.
}
\usage{
ParallelAnalysisEREC(X, r, cor = 2, pa_reps = 500, percentile = 95)
}
\arguments{
  \item{X}{Numeric matrix or data frame of observed scores
    (rows = subjects, columns = items).}
  \item{r}{Number of common factors.}
  \item{cor}{Dispersion type used for the observed data:
    \code{1} = covariance, \code{2} = Pearson correlation (default),
    \code{3} = polychoric correlation.}
  \item{pa_reps}{Number of Monte Carlo replications (default \code{500},
    as recommended in Ferrando et al. 2022). Each replication requires a full
    \code{\link{LargestDoubletsEREC}} run on a simulated dataset; reduce for
    speed in exploratory analyses.}
  \item{percentile}{Percentile used as the PA threshold (default \code{95},
    i.e. the 95th centile of the null distribution at each rank position).
    Ferrando et al. (2022) report that both the mean and the 95th centile
    produce similar results, with the mean being slightly more accurate; the
    95th centile is the more conservative choice.}
}
\value{
A named list with:
\itemize{
  \item \code{doublets}: numeric matrix (nd x 3) with columns P1, P2, and the
    observed EREC value for each doublet position.
  \item \code{pa_thresholds}: numeric vector (length nd) with the PA cutoff
    value at each rank position (percentile of the null distribution).
  \item \code{significant}: logical vector (length nd); \code{TRUE} when the
    observed EREC exceeds the PA threshold at that rank position.
  \item \code{nd}: integer, number of doublet positions evaluated.
  \item \code{pa_reps}: number of replications used.
  \item \code{percentile}: percentile used as the cutoff.
}
}
\details{
\strong{Algorithm:}
\enumerate{
  \item Compute EREC for all m*(m-1)/2 pairs in the observed data using
    \code{\link{LargestDoubletsEREC}}.
  \item Fit a baseline MINRES model (no doublets) to obtain loadings and
    uniqueness estimates. Construct the null population correlation matrix
    \eqn{R_0 = \hat\Lambda\hat\Lambda' + \hat\Psi}.
  \item Generate a pool of 100,000 null observations via Cholesky decomposition
    of \eqn{R_0} (same approach used in \code{\link{LoseferMorganaH0}}).
  \item For each of \code{pa_reps} replications: draw N rows at random and
    run \code{\link{LargestDoubletsEREC}} with Pearson correlation (appropriate
    for continuous normally distributed null data regardless of the \code{cor}
    setting of the observed data). Store the nd ranked EREC values.
  \item Compute the \code{percentile}-th quantile across replications for each
    rank position to obtain \code{pa_thresholds}.
  \item A doublet at rank position k is flagged as \code{significant} when its
    observed EREC > \code{pa_thresholds[k]}.
}

\strong{Note on computational cost:} each replication requires m*(m-1)/2
MINRES fits. For m = 20 items this means 190 fits per replication, and 95,000
fits for \code{pa_reps = 500}. For large item sets, consider starting with
\code{pa_reps = 100}.

\strong{Note on the null model for polychoric data (\code{cor = 3}):} the
null samples are generated as continuous normal data and EREC is computed
with Pearson correlation, regardless of the \code{cor} argument. The original
MORGANA-3 reference implementation includes a separate graded variant
(\code{LoseferMorganaH0Graded}) for ordinal data; that variant is not yet
implemented in this package.

The same Parallel Analysis procedure is implemented in the program
\strong{FACTOR} (\url{https://psico.fcep.urv.cat/utilitats/factor/}).
This R implementation is independent and has not been formally validated
against FACTOR.
}
\references{
Ferrando, P. J., Hernandez-Dorado, A., & Lorenzo-Seva, U. (2022).
Detecting correlated residuals in exploratory factor analysis: New proposals
and a comparison of procedures.
\emph{Structural Equation Modeling: A Multidisciplinary Journal}, 29(4), 630--638.
\doi{10.1080/10705511.2021.2004543}

Timmerman, M. E., & Lorenzo-Seva, U. (2011). Dimensionality assessment of
ordered polytomous items with parallel analysis.
\emph{Psychological Methods}, 16(2), 209--220.
\doi{10.1037/a0023353}
}
\examples{
\donttest{
db <- as.matrix(read.table(system.file("extdata", "database.dat",
                                        package = "morganafa")))

# PA com 50 replicas (reduzido para o exemplo; use 500 na pratica)
pa <- ParallelAnalysisEREC(db[1:200, ], r = 2, cor = 2,
                           pa_reps = 50, percentile = 95)

pa$doublets          # doublets observados com EREC
pa$pa_thresholds     # limiares por posicao de ranque
pa$significant       # doublets que superam o limiar de chance
}
}
\seealso{
\code{\link{LargestDoubletsEREC}},
\code{\link{solutionMEFA}}
}
