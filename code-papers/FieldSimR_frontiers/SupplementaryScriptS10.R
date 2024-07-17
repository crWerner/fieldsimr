

###################################################################
#
# Supplementary Script S10
#
# Data simulation example in Werner, Gemenet, and Tolhurst (2024)
#
# Script author: C.R. Werner
#
###################################################################


library(FieldSimR)
library(AlphaSimR)


### 1. Simulation of genetic values in 'AlphaSimR' for two additive + dominance traits tested ####
### in three environments based on an unstructured model for GxE interaction. ####################

### simulation Parameters
ntraits <- 2 # Number of traits.
nenvs <- 3 # Number of environments (locations).
nreps <- nblocks <- c(2, 2, 3) # Number of complete replicates/blocks within environments 1, 2 and 3.
nind <- 20 # Number of founder genotypes in the population.
nchr <- 10 # Number of chromosomes.
nsegsites <- 300 # Number of QTN per chromosome.

# 1.1 Define the genetic architecture of the simulated traits.
# Mean genetic values and mean dominance degrees for trait 1 in all 3 environments and trait 2
# in all 3 environments.
meanA <- c(4.9, 5.4, 5.1, 235.2, 228.5, 239.1) # Trait 1 x 3 environments, trait 2 x 3 environments.
meanDD <- c(0.4, 0.4, 0.4, 0.1, 0.1, 0.1)

# Additive genetic variances (useVarA = TRUE) and dominance degree variances for traits 1 and 2,
# assuming a separable structure between traits and environments.
varA <- c(0.085, 0.12, 0.06, 15.1, 8.5, 11.7)

# Dominance degree variances for trait 1 in 3 environments and for trait 2 in 3 environments,
# assuming a non-separable structure between traits and environments.
varDD <- rep(0.2, 6)

# Additive genetic correlations between the two simulated traits.
TcorA <- matrix(
  c(
    1.0, 0.6,
    0.6, 1.0
  ),
  ncol = 2
)

# Additive genetic correlations between the three simulated environments.
EcorA <- matrix(
  c(
    1.0, 0.4, 0.6,
    0.4, 1.0, 0.5,
    0.6, 0.5, 1.0
  ),
  ncol = 3
)

corDD <- diag(6)

input_asr <- unstr_asr_input(
  nenvs = nenvs,
  ntraits = ntraits,
  mean = meanA,
  var = varA,
  TcorA = TcorA,
  EcorA = EcorA,
  meanDD = meanDD,
  varDD = varDD,
  corDD = corDD
)


# 1.2 Use input_asr to simulate genetic values in 'AlphaSimR' based on an unstructured model for
# GxE interaction.

founders <- runMacs( # Simulation of founder genotypes using AlphaSimR's "MAIZE" presets
  nInd = nind, # to mimic the species' evolutionary history.
  nChr = nchr,
  segSites = nsegsites,
  inbred = FALSE,
  species = "MAIZE"
)

SP <- SimParam$new(founders)

SP$addTraitAD( # Additive + dominance trait simulation.
  nQtlPerChr = nsegsites,
  mean = input_asr$mean,
  var = input_asr$var,
  meanDD = input_asr$meanDD,
  varDD = input_asr$varDD,
  corA = input_asr$corA,
  corDD = input_asr$corDD,
  useVarA = FALSE
)

founders <- newPop(founders)

pool_A <- makeDH(founders[1:10], nDH = 1)  # Pool A: 1 DH line from founders 1 to 10, respectively.
pool_B <- makeDH(founders[11:20], nDH = 1) # Pool B: 1 DH line from founders 11 to 20, respectively.

dh_lines <- mergePops(list(pool_A, pool_B))

factorial_plan <- as.matrix(expand.grid(A = pool_A@id, B = pool_B@id)) # Factorial crossing plan.

hybrid_pop <- makeCross(pop = dh_lines, crossPlan = factorial_plan, nProgeny = 1) # Hybrid genotypes.

# 1.3 Create a data frame containing the simulated genetic values for the two traits
# in the three environments.

gv_df <- unstr_asr_output(
  pop = hybrid_pop,
  ntraits = ntraits,
  nenvs = nenvs,
  nreps = nreps
)

# A similar data frame is provided in FieldSimR's object 'gv_df_unstr'


### 2. Simulation of plot errors using FieldSimR #################################################

error_ls <- field_trial_error(
  ntraits = 2,
  nenvs = 3,
  nblocks = c(2, 2, 3),
  ncols = c(10, 10, 15),
  nrows = 20,
  block.dir = "col",
  varR = c(0.20, 0.28, 0.14, 15.1, 8.5, 11.7),
  spatial.model = "Bivariate",
  complexity = 10,
  plot.length = 8,
  plot.width = 2,
  prop.spatial = 0.4,
  ext.ord = "zig-zag",
  ext.dir = "row",
  prop.ext = 0.2,
  return.effects = TRUE
)

# A similar data frame is provided in FieldSimR's object 'error_df_bivar'


### 3. Generation of phenotypes through combination of simulated genetic values and plot errors ##

pheno_df <- make_phenotypes(
  gv.df = gv_df,
  error.df = error_ls$error.df,
  randomise = TRUE
)


### Extract and plot the total error, the three error components, and the phenotypic values in
### Environment 1.

total_error_1 <- error_ls$error.df[error_ls$error.df$env == 1, ]
error_components_1 <- error_ls$Trait1[error_ls$Trait1$env == 1, ]
pheno_env_1 <- pheno_df[pheno_df$env == 1, ]

plot_effects(df = total_error_1, effect = "e.Trait1")
plot_effects(df = error_components_1, effect = "e.spat")
plot_effects(df = error_components_1, effect = "e.rand")
plot_effects(df = error_components_1, effect = "e.ext.row")
plot_effects(df = pheno_env_1, effect = "y.Trait1")


save.image()

# end of script

