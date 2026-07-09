morgana_example_data <- function() {
  as.matrix(read.table(
    system.file("extdata", "database.dat", package = "morganafa")
  ))
}
