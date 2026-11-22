# morganafa

<!-- badges: start -->
[![R-CMD-check](https://github.com/pablo-huascar/morgana-factor-analysis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pablo-huascar/morgana-factor-analysis/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

Pacote R que implementa o método **MORGANA** (*Minimum expected bias Of Residual and loadinG vAlues iN sAmple estimates*) para análise fatorial exploratória e confirmatória com resíduos correlacionados (*doublets*), incluindo a detecção de doublets pelo índice **EREC**, limiares de referência por Análise Paralela e índices de ajuste corrigidos pelo procedimento **LOSEFER**.

*An unofficial R implementation of the MORGANA method for factor analysis with correlated residuals (doublets), proposed by Ferrando, Hernández-Dorado & Lorenzo-Seva. See the disclaimer below.*

---

## ⚠️ Aviso importante / Important disclaimer

**PT:** A ideia e o desenvolvimento originais dos métodos MORGANA, EREC e LOSEFER são de **Pere J. Ferrando, Ana Hernández-Dorado e Urbano Lorenzo-Seva** (Universitat Rovira i Virgili). Este pacote é apenas uma tentativa independente de tornar a proposta mais acessível como pacote de R, **sem garantias de qualquer natureza**. Os autores dos artigos **não estão envolvidos neste projeto e não têm qualquer responsabilidade sobre ele**.

Os métodos estão oficialmente implementados no programa [FACTOR](https://psico.fcep.urv.cat/utilitats/factor/), desenvolvido e mantido pelos autores, que também disponibilizaram um código R de referência. **Ainda não foram realizados testes formais comparando os resultados deste pacote com a implementação oficial do FACTOR nem com o código R disponibilizado pelos autores.** Diferenças numéricas podem existir por escolhas de implementação (tolerâncias de convergência, geração de números aleatórios, correções aplicadas — veja [Diferenças em relação ao código de referência](#diferenças-em-relação-ao-código-de-referência-dos-autores)). Sempre que possível, valide resultados importantes com o programa FACTOR.

**EN:** The MORGANA, EREC and LOSEFER methods were conceived and developed by **Pere J. Ferrando, Ana Hernández-Dorado and Urbano Lorenzo-Seva** (Universitat Rovira i Virgili). This package is an independent attempt to make their proposal more accessible as an R package, **provided as is, with no warranty of any kind**. The original authors **are not involved in this project and bear no responsibility for it**. The methods are officially implemented in the [FACTOR](https://psico.fcep.urv.cat/utilitats/factor/) program; **no formal comparison between this package and the official implementations has been conducted yet.**

---

## Instalação

O pacote pode ser instalado diretamente do GitHub:

```r
# install.packages("pak")
pak::pak("pablo-huascar/morgana-factor-analysis")

# ou, com remotes:
# install.packages("remotes")
remotes::install_github("pablo-huascar/morgana-factor-analysis")
```

## Base metodológica

- **Ferrando, P. J., Hernández-Dorado, A., & Lorenzo-Seva, U. (2022).** Detecting correlated residuals in exploratory factor analysis: New proposals and a comparison of procedures. *Structural Equation Modeling: A Multidisciplinary Journal, 29*(4), 630–638. <https://doi.org/10.1080/10705511.2021.2004543>
- **Lorenzo-Seva, U. (2025).** Mild-restricted confirmatory factor analysis: A new approach to assessing population factor models. *Communications in Statistics – Simulation and Computation*. <https://doi.org/10.1080/03610918.2025.2453473>

### Visão geral do método

Na análise fatorial clássica, assume-se que os resíduos são não correlacionados. Em análises de itens, porém, é comum haver pares de itens com resíduos correlacionados (*doublets*) — por redundância de conteúdo, por exemplo — o que, se ignorado, enviesa cargas e resíduos.

O MORGANA estima os doublets mais substanciais em três estágios: (1) para cada par candidato, o modelo é ajustado nas m−2 variáveis restantes (*core set*) via MINRES; (2) as cargas do par são estimadas por análise de extensão; (3) a correlação residual do par é estimada a partir da matriz residual reescalonada. O índice resultante, **EREC** (*Expected REsidual Correlation direct change*; EPC1 no artigo), quantifica a correlação residual esperada de cada par. Os maiores EREC — limitados pelo *Ledermann bound* — são liberados no ajuste final. Os índices de ajuste (RMSEA, CFI, NNFI/TLI) são corrigidos pelo procedimento **LOSEFER**, que aproxima por Monte Carlo a distribuição empírica do qui-quadrado sob H0.

## Uso rápido

### EFA com detecção automática de doublets

```r
library(morganafa)

database <- as.matrix(read.table(
  system.file("extdata", "database.dat", package = "morganafa")
))

set.seed(123)
efa <- solutionMEFA(
  database,
  r          = 2,          # número de fatores
  cor        = "pearson",  # "pearson", "polychoric" ou "covariance"
  rotate     = "oblimin",  # rotação para interpretação ("varimax", "none")
  montecarlo = 5000,       # réplicas para o H0 do LOSEFER
  pa_reps    = 500         # réplicas da Análise Paralela (0 desativa)
)

efa            # impressão amigável de toda a solução
efa$loadings   # tibble: cargas rotacionadas, h2 e u2 por item
efa$doublets   # tibble: doublets, EREC, correlação residual, limiar da PA
efa$fit        # tibble: qui-quadrado LOSEFER, RMSEA, CFI, NNFI/TLI
efa$pa         # resultado completo da Análise Paralela
```

### Vinculando itens ao seu texto

Os nomes das colunas dos dados são usados como nomes dos itens em todas as
saídas. Opcionalmente, é possível fornecer o texto dos itens:

```r
rotulos <- c(
  V1 = "Me sinto nervoso ultimamente",
  V2 = "Tenho dificuldade para relaxar",
  V6 = "Sinto-me triste sem motivo",
  V7 = "Perdi o interesse pelas coisas"
)

efa <- solutionMEFA(database, r = 2, item_labels = rotulos, pa_reps = 0)
efa$doublets
#> # A tibble: 4 × 9
#>    rank item_a item_b label_a            label_b            erec residual_r ...
#> 1     1 V1     V2     Me sinto nervoso…  Tenho dificuldade… 0.501      0.511
#> 2     2 V6     V7     Sinto-me triste…   Perdi o interesse… 0.493      0.499
```

### CFA (mild-restricted) com doublets definidos pelo usuário

```r
set.seed(123)
cfa <- solutionMCFA(
  database,
  r        = 2,
  doublets = rbind(c("V1", "V2"),   # aceita nomes de itens…
                   c(6, 7)),        # …ou índices
  montecarlo = 5000
)
cfa$fit
```

### Análise Paralela para o EREC de forma independente

```r
set.seed(123)
pa <- ParallelAnalysisEREC(database, r = 2, pa_reps = 500, percentile = 95)
pa
```

Um doublet é considerado *significativo* quando seu EREC observado supera o
limiar da Análise Paralela na mesma posição de ranque (seção *"Reference
values and upper bounds"* de Ferrando et al., 2022).

## Tempo de processamento

Quase todo o tempo de `solutionMEFA()` vai para duas simulações Monte Carlo:

* **Análise Paralela** (`pa_reps`, padrão 500): cada réplica exige
  m·(m−1)/2 ajustes de modelo, então o tempo cresce com o **quadrado do
  número de itens**.
* **Distribuição H0 do LOSEFER** (`montecarlo`, padrão 5000).

Referência aproximada (notebook Apple silicon, 952 casos, 58 itens, r = 3):
cerca de 40 min para a PA e 8 min para o LOSEFER com os valores padrão. Com
20 itens, cerca de 1,5 min no total. Barras de progresso são exibidas em
sessões interativas.

Recomendações:

* **Use vários cores** com o argumento `cores` (ou defina uma vez por
  sessão com `options(morganafa.cores = 8)`). Os resultados são
  **idênticos** aos da execução sequencial com a mesma `set.seed()`,
  qualquer que seja o número de cores. Com 58 itens e 9 cores, a PA ficou
  cerca de 5× mais rápida (≈ 7 min em vez de ≈ 37 min). Funciona no
  Positron, no RStudio e no Windows. Iniciar os processos leva alguns
  segundos, então só compensa em execuções longas.

* Em explorações, use valores menores (`pa_reps = 100`,
  `montecarlo = 1000`, ou `pa_reps = 0` para pular a PA); mantenha os
  valores padrão na análise final.
* A matriz policórica é calculada só uma vez por chamada, mas cada função
  chamada separadamente a recalcula. Ao usar várias funções nos mesmos
  dados, calcule-a uma vez e passe-a pelo argumento `R`:

```r
R_poly <- psych::polychoric(dados)$rho
LargestDoubletsEREC(dados, r = 3, cor = 3, R = R_poly)
ParallelAnalysisEREC(dados, r = 3, cor = 3, R = R_poly)
```

## Funções

| Função | Descrição |
|---|---|
| `solutionMEFA()` | Fluxo completo de EFA (EREC + MORGANA + LOSEFER + PA), com saídas em tibbles |
| `solutionMCFA()` | Fluxo completo de CFA *mild-restricted* (doublets do usuário) |
| `LargestDoubletsEREC()` | Identifica e ordena os maiores doublets pelo EREC |
| `DoubletsEREC()` | Calcula o EREC para pares especificados pelo usuário |
| `ParallelAnalysisEREC()` | Limiares de referência para o EREC via Análise Paralela |
| `MorganaEFA()` / `MorganaCFA()` / `Morgana()` | Rotinas de extração de nível mais baixo |
| `LoseferMorgana()` / `LoseferMorganaH0()` | Correção LOSEFER dos índices de ajuste |
| `posit()` | Reflexão de sinais das colunas de uma matriz de cargas |

As funções de nível mais baixo mantêm a estrutura do código de referência dos
autores e são exportadas para permitir fluxos customizados; para uso comum,
prefira `solutionMEFA()`/`solutionMCFA()`.

## Diferenças em relação ao código de referência dos autores

Este pacote partiu do código R de referência disponibilizado pelos autores.
Na revisão, alguns problemas aparentes de tradução do MATLAB original foram
corrigidos — o que significa que **os resultados deste pacote podem diferir
dos do código R de referência** (mais um motivo para a validação formal
pendente contra o programa FACTOR):

1. **`fufa1()` (simulação H0 do LOSEFER):** o `switch()` numérico só
   funcionava para correlações de Pearson; e os resíduos eram padronizados
   por √h² (comunalidade) em vez de √(1−h²) (unicidade). A padronização pela
   unicidade é a consistente com o estimador usado nos dados reais
   (`Morgana()`), e altera levemente o qui-quadrado ajustado.
2. **`DoubletsEREC()`:** o bloco que acumula as unicidades para o
   reescalonamento do EREC ficava fora do laço interno de pares (diferente de
   `LargestDoubletsEREC()`); com a correção, o EREC de um mesmo par é idêntico
   nos fluxos exploratório e confirmatório.
3. **`BERGE()` (caso Heywood no MINRES):** erros de indexação
   (`1:m-q` em vez de `1:(m-q)`) e de produto matricial faziam a função
   falhar com mais de um fator; foi reescrita resolvendo as mesmas equações
   estacionárias de Ten Berge & Nevels (1977).
4. **`LoseferMorganaH0()`:** laços `while` descartavam a última réplica de
   Monte Carlo; matriz quase singular era "corrigida" somando N (o tamanho
   amostral!) à diagonal — agora o menor autovalor é elevado a 0,001.
5. **`posit()`:** a regra de reflexão de sinais continha código morto; agora
   reflete colunas com soma de cargas negativa.

As partes determinísticas do algoritmo (detecção EREC, extração MORGANA,
cargas) reproduzem o código de referência com precisão de máquina (≈1e−16)
nos dados de exemplo.

## Licença

GPL (>= 3). Veja o arquivo [LICENSE.md](LICENSE.md).

## Como citar

Se você usar este pacote, cite **os artigos originais dos autores do método**
(referências acima). Para citar o pacote em si:

> Aragão Pinheiro, F. P. H. (2026). *morganafa: MORGANA factor analysis with
> correlated residuals* (implementação não oficial em R). R package version 0.1.3.
> <https://github.com/pablo-huascar/morgana-factor-analysis>
