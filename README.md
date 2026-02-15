# morganafa

Pacote R para Análise Fatorial MORGANA (EFA/CFA) com resíduos correlacionados.

## Instalação

```r
install.packages("remotes")
remotes::install_github("SEU_USUARIO/morgana-factor-analysis")
```

## Uso rápido

```r
library(morganafa)

database <- as.matrix(
  read.table(system.file("extdata", "database.dat", package = "morganafa"))
)

# EFA (detecção automática de doublets com EREC)
res_efa <- solutionMEFA(
  database = database,
  r = 2,
  cor = 2,
  conver = 0.0001,
  maxiter = 500,
  display = 1
)

# CFA (doublets fornecidos pelo usuário)
doublets <- matrix(c(1, 2, 6, 7), ncol = 2, byrow = TRUE)
res_cfa <- solutionMCFA(
  database = database,
  r = 2,
  doublets = doublets,
  cor = 2,
  conver = 0.0001,
  maxiter = 500,
  display = 1
)
```
