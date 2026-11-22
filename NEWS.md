# morganafa 0.1.1

## Bug fixes

* Polychoric correlations (`cor = 3` / `"polychoric"`) no longer fail in
  Positron with "Can't fork the R session". `psych::polychoric()`
  parallelises via `parallel::mcmapply()`, which forks the R session, and
  Positron forbids forking. When running inside Positron, morganafa now
  computes the polychoric matrix with `mc.cores = 1` (restored afterwards).
  Other front-ends keep `psych`'s default parallel behaviour. Affects every
  function that estimates the polychoric matrix internally
  (`LargestDoubletsEREC()`, `DoubletsEREC()`, `MorganaEFA()`,
  `solutionMEFA()`, `solutionMCFA()`, etc.).

# morganafa 0.1.0

First public (GitHub) release. Unofficial R implementation of the MORGANA /
EREC / LOSEFER methodology by Ferrando, Hernández-Dorado & Lorenzo-Seva.

## User-facing interface

* `solutionMEFA()` and `solutionMCFA()` were redesigned: they now accept
  data frames or matrices, `cor` as a string (`"pearson"`, `"polychoric"`,
  `"covariance"`) or the original numeric codes, optional factor rotation
  (`"oblimin"` default for r > 1, via GPArotation; also `"varimax"` and
  `"none"`), and optional `item_labels` with the text of each item.
* Results are classed objects (`morgana_efa` / `morgana_cfa`) with tidy
  tibble components (`$loadings`, `$doublets`, `$fit`), item names carried
  throughout, and a readable `print()` method. `ParallelAnalysisEREC()`
  gained a `morgana_pa` class and print method.
* `solutionMCFA()` accepts doublet pairs as item names or indices, and
  validates them.

## Bug fixes relative to the authors' released R code

These change results relative to the reference code (see README):

* `fufa1()`: numeric `switch()` only worked for Pearson correlations;
  residuals were standardised by the communality instead of the uniqueness.
* `DoubletsEREC()`: the uniqueness-collection block sat outside the inner
  pair loop; EREC values now match `LargestDoubletsEREC()` for the same
  pairs. The EREC of user pairs is now reported in absolute value, as in
  the exploratory flow.
* `BERGE()`: `1:m-q`-style indexing and scalar `%*%` matrix errors made the
  Heywood-case correction fail for r > 1; rewritten via the stationary
  equations of Ten Berge & Nevels (1977).
* `LoseferMorganaH0()`: off-by-one `while` loops left one Monte Carlo
  replicate empty; the near-singularity repair added N to the diagonal
  instead of raising the smallest eigenvalue to 0.001.
* `posit()`: sign-reflection rule contained dead code; now reflects columns
  whose loading sum is negative.

## Performance

* The dispersion matrix (correlation/covariance) is computed once per
  analysis and passed down (`R` argument on the low-level functions) —
  a large saving for polychoric correlations.
* Monte Carlo subsampling uses `sample.int()` instead of sorting 100,000
  uniform draws per replication; various loops were vectorised.

## Infrastructure

* Full roxygen2 documentation with method references and disclaimers;
  testthat test suite; GitHub Actions R CMD check.
