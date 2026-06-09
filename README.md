# morganafa

Pacote R que implementa o método **MORGANA** (*Minimum expected bias Of Residual and loadinG vAlues iN sAmple estimates*) para análise fatorial exploratória e confirmatória com resíduos correlacionados.

---

## Base metodológica

O pacote é uma implementação independente em R das propostas descritas nos seguintes artigos:

- **Ferrando, P. J., Hernandez-Dorado, A., & Lorenzo-Seva, U. (2022).** Detecting correlated residuals in exploratory factor analysis: New proposals and a comparison of procedures. *Structural Equation Modeling: A Multidisciplinary Journal*, 29(4), 630–638. <https://doi.org/10.1080/10705511.2021.2004543>

- **Lorenzo-Seva, U. (2025).** Mild-restricted confirmatory factor analysis: a new approach to assessing population factor models. *Communications in Statistics – Simulation and Computation*. <https://doi.org/10.1080/03610918.2025.2453473>

Os mesmos métodos estão implementados no programa **FACTOR** (disponível gratuitamente em <https://psico.fcep.urv.cat/utilitats/factor/>), desenvolvido e mantido pelos autores dos artigos na Universitat Rovira i Virgili.

> **Aviso importante:** Este pacote R foi desenvolvido de forma **independente** do programa FACTOR e **não foi submetido a uma comparação formal** entre os resultados produzidos por ambas as implementações. Diferenças numéricas podem existir em virtude de escolhas de implementação distintas (tolerâncias de convergência, geração de números aleatórios, ordem de operações). Recomenda-se cautela ao comparar os resultados e, sempre que possível, validar os achados com o programa FACTOR ou com dados simulados de referência.

---

## Visão geral do método

Na Análise Fatorial Exploratória (EFA) clássica, os resíduos são assumidos não-correlacionados. Na prática — especialmente em análises de itens — podem existir pares de variáveis com resíduo correlacionado (*doublets*), o que, se ignorado, produz cargas fatoriais e resíduos enviesados.

O MORGANA propõe um procedimento ULS em três estágios para estimar os doublets mais substanciais e corrigi-los:

1. **Ajuste do conjunto núcleo:** para cada par candidato (i, j), o modelo é ajustado nas m−2 variáveis restantes (*core set*) via MINRES.  
2. **Análise de extensão:** as cargas de i e j são estimadas por extensão a partir do núcleo ajustado.  
3. **Estimação dos resíduos:** a correlação residual do par é estimada via Ψ̂⁻¹(R − Λ̂Λ̂ʹ)Ψ̂⁻¹.

O índice resultante **EREC** (*Expected REsidual Correlation direct Change*, EPC1 no artigo) quantifica a correlação residual estimada para cada par. Os pares com maiores valores de EREC, limitados pelo *Ledermann bound*, são liberados como doublets antes do ajuste final do modelo com a função `Morgana()`.

Os índices de ajuste (RMSEA, CFI, NNFI/TL) são corrigidos pelo procedimento **LOSEFER** (*LOw-Sampling-Error Fit Empirical-Reference*), que usa simulação de Monte Carlo para aproximar a distribuição empírica da estatística qui-quadrado.

---

## Instalação

```r
# A partir do código-fonte local (diretório do projeto):
install.packages(".", repos = NULL, type = "source")

# Ou via R CMD no terminal:
# R CMD INSTALL .
```

---

## Funções principais

| Função | Descrição |
|---|---|
| `solutionMEFA()` | Fluxo completo de EFA (detecção automática de doublets + MORGANA + LOSEFER + PA) |
| `solutionMCFA()` | Fluxo completo de CFA (doublets fornecidos pelo usuário + MORGANA + LOSEFER) |
| `LargestDoubletsEREC()` | Identifica e ordena os maiores doublets por EREC; aplica restrição de Ledermann |
| `DoubletsEREC()` | Calcula EREC para pares explicitamente especificados pelo usuário |
| `ParallelAnalysisEREC()` | Gera limiares de referência para EREC via Análise Paralela (Monte Carlo) |
| `MorganaEFA()` | Extração MORGANA em modo exploratório |
| `MorganaCFA()` | Extração MORGANA em modo confirmatório |
| `Morgana()` | Rotina central de estimação (procedimento ULS em três estágios) |
| `LoseferMorgana()` | Calcula RMSEA, CFI e NNFI corrigidos pelo LOSEFER |
| `LoseferMorganaH0()` | Estima parâmetros da distribuição H0 por Monte Carlo |

---

## Uso rápido

### EFA com detecção automática de doublets e Análise Paralela

```r
library(morganafa)

database <- as.matrix(
  read.table(system.file("extdata", "database.dat", package = "morganafa"))
)

# EFA completa: EREC + MORGANA + LOSEFER + Parallel Analysis (500 réplicas)
res_efa <- solutionMEFA(
  database    = database,
  r           = 2,
  cor         = 2,       # 1 = covariância, 2 = Pearson, 3 = policórica
  montecarlo  = 5000,    # réplicas para LOSEFER H0
  pa_reps     = 500      # réplicas para PA dos limiares EREC; use 0 para desativar
)

# Doublets detectados e seus valores EREC
res_efa$doublets

# Limiares da Análise Paralela por posição de ranque (percentil 95)
res_efa$pa$pa_thresholds

# Doublets que superam o limiar de chance
res_efa$pa$significant

# Índices de ajuste
res_efa$fit$RMSEA
res_efa$fit$CFI
res_efa$fit$TL
```

### CFA com doublets fornecidos pelo usuário

```r
doublets <- matrix(c(1, 2,
                     6, 7), ncol = 2, byrow = TRUE)

res_cfa <- solutionMCFA(
  database   = database,
  r          = 2,
  doublets   = doublets,
  cor        = 2,
  montecarlo = 5000
)

res_cfa$fit$RMSEA
```

### Análise Paralela para EREC de forma independente

```r
pa <- ParallelAnalysisEREC(
  X          = database,
  r          = 2,
  cor        = 2,
  pa_reps    = 500,
  percentile = 95        # percentil de referência (padrão: 95)
)

# Doublets observados com limiares e flags de significância
pa$doublets
pa$pa_thresholds
pa$significant
```

---

## Notas sobre a Análise Paralela (ParallelAnalysisEREC)

O procedimento de Análise Paralela para o índice EREC é descrito na seção *"Reference values and upper bounds"* de Ferrando et al. (2022). Sua lógica é:

1. Ajustar o modelo sem doublets (apenas fatores comuns) para obter a matriz de correlações nula.  
2. Gerar `pa_reps` amostras de tamanho N a partir dessa população nula.  
3. Para cada réplica, calcular os `nd` maiores valores de EREC.  
4. O limiar para a posição de ranque k é o percentil escolhido (padrão: 95) da distribuição simulada de EREC na posição k.

Um doublet é considerado *significativo* quando seu EREC observado supera o limiar da Análise Paralela na mesma posição de ranque.

> **Nota de desempenho:** cada réplica exige m*(m−1)/2 ajustes de modelo. Para matrizes grandes (m > 20), `pa_reps = 500` pode ser computacionalmente intensivo. Considere usar `pa_reps = 100` em explorações iniciais.

---

## Referências

Ferrando, P. J., Hernandez-Dorado, A., & Lorenzo-Seva, U. (2022). Detecting correlated residuals in exploratory factor analysis: New proposals and a comparison of procedures. *Structural Equation Modeling: A Multidisciplinary Journal*, 29(4), 630–638. <https://doi.org/10.1080/10705511.2021.2004543>

Lorenzo-Seva, U. (2025). Mild-restricted confirmatory factor analysis: a new approach to assessing population factor models. *Communications in Statistics – Simulation and Computation*. <https://doi.org/10.1080/03610918.2025.2453473>
