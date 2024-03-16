## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include = FALSE---------------------------------------------------
library(AlphaSimR)
library(FieldSimR)
library(ggplot2)

## -----------------------------------------------------------------------------
n_traits <- 2 # Number of traits.
n_envs <- 3 # Number of environments (locations).
n_reps <- c(2, 3, 2) # Number of replicates tested within environments 1, 2 and 3.


n_ind <- 100 # Number of founder genotypes in the population.
n_chr <- 10 # Number of chromosomes.
n_seg_sites <- 200 # Number of QTN per chromosome.

## -----------------------------------------------------------------------------
mean <- c(4.9, 5.4, 5.1, 235.2, 228.5, 239.1) # c(Yld:E1, Yld:E2, Yld:E3, Pht:E1, Pht:E2, Pht:E3)

## -----------------------------------------------------------------------------
var <- c(0.08, 13) # c(grain yield, plant height)

## -----------------------------------------------------------------------------
rel_main_eff_A <- c(0.4, 0.6) # c(grain yield, plant height)

## -----------------------------------------------------------------------------
cor_A <- matrix( # Matrix of additive genetic correlations grain yield and plant height.
  c(
    1.0, 0.5,
    0.5, 1.0
  ),
  ncol = 2
)

## ----echo=FALSE---------------------------------------------------------------
cor_A

## -----------------------------------------------------------------------------
mean_DD <- c(0.4, 0.4, 0.4, 0.1, 0.1, 0.1) # c(Yld:E1, Yld:E2, Yld:E3, Pht:E1, Pht:E2, Pht:E3)

var_DD <- c(0.2, 0.2) # c(grain yield, plant height)

rel_main_eff_DD <- 0.4 # Same value set for traits 1 and 2.

cor_DD <- diag(2)

## ----echo=FALSE---------------------------------------------------------------
cor_DD

## -----------------------------------------------------------------------------
input_asr <- compsym_asr_input(
  n_envs = n_envs,
  n_traits = n_traits,
  mean = mean,
  var = var,
  rel_main_eff_A = rel_main_eff_A,
  cor_A = cor_A,
  mean_DD = mean_DD,
  var_DD = var_DD,
  rel_main_eff_DD = rel_main_eff_DD,
  cor_DD = cor_DD
)

## -----------------------------------------------------------------------------
founders <- runMacs( # Simulation of founder genotypes using AlphaSimR's "MAIZE" presets
  nInd = n_ind, # to mimic the species' evolutionary history.
  nChr = n_chr,
  segSites = n_seg_sites,
  species = "MAIZE",
  nThreads = 2
)

SP <- SimParam$new(founders)

## -----------------------------------------------------------------------------
SP$addTraitAD( # Additive + dominance trait simulation.
  nQtlPerChr = n_seg_sites,
  mean = input_asr$mean,
  var = input_asr$var,
  meanDD = input_asr$mean_DD,
  varDD = input_asr$var_DD,
  corA = input_asr$cor_A,
  corDD = input_asr$cor_DD,
  useVarA = FALSE
)

founders <- newPop(founders)

## -----------------------------------------------------------------------------
pool_A <- makeDH(founders[1:10], nDH = 1) # Pool A: 1 DH line from founders 1 to 10, respectively.
pool_B <- makeDH(founders[11:20], nDH = 1) # Pool B: 1 DH line from founders 11 to 20, respectively.

dh_lines <- mergePops(list(pool_A, pool_B))

factorial_plan <- as.matrix(expand.grid(A = pool_A@id, B = pool_B@id)) # Factorial crossing plan.

hybrid_pop <- makeCross(pop = dh_lines, crossPlan = factorial_plan, nProgeny = 1) # Hybrid genotypes.

## -----------------------------------------------------------------------------
gv_df <- compsym_asr_output(
  pop = hybrid_pop,
  n_envs = n_envs,
  n_reps = n_reps,
  n_traits = n_traits
)

## ----echo=FALSE---------------------------------------------------------------
head(gv_df)

## ----echo=TRUE, fig.height = 4, fig.width = 9, fig.align = "center"-----------
ggplot(gv_df, aes(x = gv.Trait.1, fill = factor(env))) +
  geom_histogram(color = "#e9ecef", alpha = 0.8, position = "identity", bins = 50) +
  scale_fill_manual(values = c("violetred3", "goldenrod3", "skyblue2")) +
  labs(x = "Genetic values for grain yield (t/ha)", y = "Count", fill = "Environment")

