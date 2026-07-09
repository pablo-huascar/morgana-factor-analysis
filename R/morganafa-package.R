#' morganafa: MORGANA Factor Analysis with Correlated Residuals
#'
#' Unofficial R implementation of the MORGANA methodology (Minimum expected
#' bias Of Residual and loadinG vAlues iN sAmple estimates) for exploratory
#' and confirmatory factor analysis with correlated residuals (doublets),
#' including EREC-based doublet detection, Parallel Analysis reference
#' values for EREC, and LOSEFER-corrected fit indices.
#'
#' @section Important disclaimer:
#' The MORGANA, EREC, and LOSEFER methods were proposed and developed by
#' Pere J. Ferrando, Ana Hernandez-Dorado, and Urbano Lorenzo-Seva
#' (Universitat Rovira i Virgili). **This package is an independent effort
#' to make their proposal more accessible as an R package. The original
#' authors are not involved in this package and bear no responsibility for
#' it.**
#'
#' The methods are officially implemented in the FACTOR program
#' (<https://psico.fcep.urv.cat/utilitats/factor/>), and the authors also
#' released reference R code. **No formal comparison between the results of
#' this package and the official implementations has been conducted yet.**
#' The package is provided as is, without any warranty; validate important
#' results against the FACTOR program whenever possible.
#'
#' @section Main functions:
#' * [solutionMEFA()] - complete exploratory workflow (doublet detection +
#'   MORGANA + LOSEFER + Parallel Analysis).
#' * [solutionMCFA()] - complete confirmatory (mild-restricted) workflow
#'   with user-specified doublets.
#' * [LargestDoubletsEREC()], [DoubletsEREC()], [ParallelAnalysisEREC()],
#'   [Morgana()], [MorganaEFA()], [MorganaCFA()], [LoseferMorgana()],
#'   [LoseferMorganaH0()] - lower-level building blocks kept close to the
#'   authors' reference code.
#'
#' @references
#' Ferrando, P. J., Hernandez-Dorado, A., & Lorenzo-Seva, U. (2022).
#' Detecting correlated residuals in exploratory factor analysis: New
#' proposals and a comparison of procedures. *Structural Equation Modeling:
#' A Multidisciplinary Journal*, 29(4), 630-638.
#' \doi{10.1080/10705511.2021.2004543}
#'
#' Lorenzo-Seva, U. (2025). Mild-restricted confirmatory factor analysis: a
#' new approach to assessing population factor models. *Communications in
#' Statistics - Simulation and Computation*.
#' \doi{10.1080/03610918.2025.2453473}
#'
#' @keywords internal
#' @aliases morganafa
"_PACKAGE"

#' @importFrom stats cor cov median quantile rnorm runif sd setNames uniroot var
#' @importFrom tibble tibble
NULL
