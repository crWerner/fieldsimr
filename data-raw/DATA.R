# Simulation of genetic values in 'AlphaSimR' for two additive + dominance traits in
# three environments based on an unstructured model for GxE interaction.

set.seed(1312)
library(FieldSimR)

### simulation Parameters
ntraits <- 2 # Number of traits.
nenvs <- 3 # Number of environments.
nreps <- nblocks <- c(2, 2, 3) # Number of replicates/blocks in Environments 1, 2 and 3.
nind <- 20 # Number of founder genotypes in the population.
nchr <- 10 # Number of chromosomes.
nseg_sites <- 300 # Number of QTN per chromosome.

# 1. Define the genetic architecture of the simulated traits.
# Mean genetic values and mean dominance degrees
mean <- c(4.9, 5.4, 5.1, 235.2, 228.5, 239.1) # Trait 1 x 3 environments, Trait 2 x 3 environments
meanDD <- c(0.4, 0.4, 0.4, 0.1, 0.1, 0.1) # Trait 1 and 2, same values for all environments

# Additive genetic variances (useVarA = TRUE) and dominance degree variances for Traits 1 and 2,
# assuming a non-separable structure between traits and environments.
var <- c(0.086, 0.12, 0.06, 15.1, 8.5, 11.7)

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
  mean = mean,
  var = var,
  TcorA = TcorA,
  EcorA = EcorA,
  meanDD = meanDD,
  varDD = varDD,
  corDD = corDD
)


# 2. Use input_asr to simulate genetic values with 'AlphaSimR' based on an unstructured model.

library("AlphaSimR")
founders <- runMacs( # Simulation of founder genotypes using AlphaSimR's "MAIZE" presets
  nInd = nind, # to mimic the species' evolutionary history.
  nChr = nchr,
  segSites = nseg_sites,
  inbred = FALSE,
  species = "MAIZE"
)

SP <- SimParam$new(founders)

SP$addTraitAD( # Additive + dominance trait simulation.
  nQtlPerChr = nseg_sites,
  mean = input_asr$mean,
  var = input_asr$var,
  corA = input_asr$corA,
  meanDD = input_asr$meanDD,
  varDD = input_asr$varDD,
  corDD = input_asr$corDD,
  useVarA = FALSE
)

founders <- newPop(founders)

pool_A <- makeDH(founders[1:10], nDH = 1) # Pool A: 1 DH line from founders 1 to 10, respectively.
pool_B <- makeDH(founders[11:20], nDH = 1) # Pool B: 1 DH line from founders 11 to 20, respectively.

dh_lines <- mergePops(list(pool_A, pool_B))

factorial_plan <- as.matrix(expand.grid(A = pool_A@id, B = pool_B@id)) # Factorial crossing plan.

hybrid_pop <- makeCross(pop = dh_lines, crossPlan = factorial_plan, nProgeny = 1) # Hybrid genotypes.

# 3. Create a data frame containing the simulated genetic values for the two traits
#  in the three environments.

gv_df_unstr <- unstr_asr_output(
  pop = hybrid_pop,
  ntraits = ntraits,
  nenvs = nenvs,
  nreps = nreps
)



# Simulation of plot errors for two traits in three environments using a bivariate
# interpolation model for spatial variation.

block_dir <- "col" # Layout of replicates (“side-by-side”).
ncols <- c(10, 10, 15) # Total umber of columns per location.
nrows <- 20 # Total number of rows per location.
plot_length <- 8 # Plot length; here in meters (column direction).
plot_width <- 2 # Plot width; here in meters (row direction).


H2 <- c(0.3, 0.3, 0.3, 0.5, 0.5, 0.5) # c(Yld:E1, Yld:E2, Yld:E3, Pht:E1, Pht:E2, Pht:E3)

# Error variances for traits 1 and 2.
calc_varR <- function(var, H2) {
  varR <- (var / H2) - var
  return(varR)
}

varR <- calc_varR(var, H2)

# Spatial error correlations between traits 1 and 2.
spatial_model <- "Bivariate" # Spatial error model.
complexity <- 10
prop_spatial <- 0.4 # Proportion of spatial to total error variance.
ScorR <- rand_cor_mat(ntraits, min.cor = 0, max.cor = 0.5, pos.def = TRUE)

prop_ext <- 0.2
ext_dir <- "row"
ext_ord <- "zig-zag"
EcorR <- rand_cor_mat(ntraits, min.cor = 0, max.cor = 0.5, pos.def = TRUE)


# Simulate field error using bivariate interpolation.
error_df_bivar <- field_trial_error(
  ntraits = ntraits,
  nenvs = nenvs,
  nblocks = nblocks,
  block.dir = block_dir,
  ncols = ncols,
  nrows = nrows,
  plot.length = plot_length,
  plot.width = plot_width,
  varR = varR,
  RcorR = NULL,
  spatial.model = spatial_model,
  complexity = complexity,
  prop.spatial = prop_spatial,
  ScorR = ScorR,
  prop.ext = prop_ext,
  ext.dir = ext_dir,
  ext.ord = ext_ord,
  EcorR = EcorR,
  return.effects = FALSE
)

# Combine genetic values and field error to obtain phenotypes.
pheno_df <- make_phenotypes(gv_df_unstr, error_df_bivar, randomise = TRUE)

# str(gv_df_unstr)
# str(error_df_bivar)
# head(pheno_df)
#
# plot_effects(pheno_df,
#   env = 2,
#   effect = "phe.Trait2"
# )



##### IMPLEMENT GENETIC VALUES USED IN PAPER AND UPDATE IDs #####

load("genetic_values_T1_paper.RData")
pheno_df_spat <- pheno_df_spat[, c(5, 6)]
pheno_df_spat <- pheno_df_spat[order(pheno_df_spat$id), ]
gv_replace <- rep(unique(pheno_df_spat$gv), 2)
gv_df_unstr$gv.Trait1[1:200] <- gv_replace
id_new <- rep(1:100, 7)
gv_df_unstr$id <- id_new
gv_df_unstr$id <- as.factor(gv_df_unstr$id)
# reorder columns (ids followed by reps)
gv_df_unstr <- gv_df_unstr[, c(1, 3, 2, 4:5)]

##### REORDER ENVS 2 & 3 AND INSERT PAPER ERRORS FOR ENV1 #####
error_df_bivar$e.Trait1[error_df_bivar$env == 1] <- pheno_df_spat$e_total
error_df_bivar$env <- trimws(error_df_bivar$env)
error_df_bivar$env[error_df_bivar$env == 3] <- 4
error_df_bivar$env[error_df_bivar$env == 2] <- 3
error_df_bivar$env[error_df_bivar$env == 4] <- 2
error_df_bivar$env <- factor(error_df_bivar$env)
error_df_bivar <- error_df_bivar[order(error_df_bivar$env, error_df_bivar$col, error_df_bivar$row), ]
rownames(error_df_bivar) <- NULL

gv_df_unstr$env <- trimws(gv_df_unstr$env)
gv_df_unstr$env[gv_df_unstr$env == 3] <- 4
gv_df_unstr$env[gv_df_unstr$env == 2] <- 3
gv_df_unstr$env[gv_df_unstr$env == 4] <- 2
gv_df_unstr$env <- factor(gv_df_unstr$env)
gv_df_unstr <- gv_df_unstr[order(gv_df_unstr$env, gv_df_unstr$rep, gv_df_unstr$id), ]
rownames(gv_df_unstr) <- NULL


###########################################################



# usethis::use_data(gv_df_unstr, overwrite = TRUE)
# usethis::use_data(error_df_bivar, overwrite = TRUE)
