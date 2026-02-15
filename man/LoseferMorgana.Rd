\name{LoseferMorgana}
\alias{LoseferMorgana}
\alias{LoseferMorganaH0}
\title{LOSEFER Fit Correction for MORGANA Solutions}
\description{
\code{LoseferMorganaH0()} estimates empirical transformation parameters under the null model.
\code{LoseferMorgana()} applies those parameters to compute adjusted fit statistics (including RMSEA, CFI, and Tucker-Lewis index).
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
  \item{S}{Observed covariance/correlation matrix.}
  \item{A}{Loading matrix.}
  \item{N}{Sample size.}
  \item{RESCOR}{Residual correlation matrix.}
  \item{nd}{Number of freely estimated residual pairs.}
  \item{mchi, sxchi, a_n, b1_n, b2_n, a_f, b_f, c_f, d_f}{Transformation parameters from \code{LoseferMorganaH0()}.}
  \item{display}{If \code{1}, prints fit output.}
  \item{doblets}{Two-column matrix of free residual pairs used in H0 simulation.}
  \item{montecarlo}{Number of simulated replications for H0 estimation.}
}
\value{
\code{LoseferMorgana()} returns a list with adjusted fit values:
\code{Um}, \code{d}, \code{d0}, \code{C0}, \code{RMSEA}, \code{CFI}, and \code{TL}.

\code{LoseferMorganaH0()} returns transformation parameters used to adjust the fit statistic.
}
