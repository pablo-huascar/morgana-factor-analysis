# MINRES extraction with Heywood-case protection.
# See .minres_engine() in my_minres.R.
#' @noRd
my_minres_hc <- function(C, q, conv, N) {
  .minres_engine(C, q, conv, N, heywood = TRUE)
}
