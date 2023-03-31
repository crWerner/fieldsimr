# Simulation of genetic values in 'AlphaSimR' for two additive + dominance traits tested in
# three environments based on an unstructured model for GxE interaction.

set.seed(1312)
library(FieldSimR)

### simulation Parameters
n_traits <- 2 # Number of traits.
n_envs <- 3 # Number of environments (locations).
n_reps <- c(3, 3, 2) # Number of full replicates within environments 1, 2 and 3.
n_ind <- 20 # Number of founder genotypes in the population.
n_chr <- 10 # Number of chromosomes.
n_seg_sites <- 300 # Number of QTN per chromosome.

# 1. Define the genetic architecture of the simulated traits.
# Mean genetic values and mean dominance degrees for trait 1 in all 3 environments and trait 2
# in all 3 environments.
mean <- c(4.9, 5.4, 5.1, 235.2, 228.5, 239.1) # Trait 1 x 3 environments, trait 2 x 3 environments.
mean_DD <- c(0.4, 0.1) # == c(0.1, 0.1, 0.1, 0.4, 0.4, 0.4)

# Additive genetic variances (useVarA = TRUE) and dominance degree variances for traits 1 and 2,
# assuming a separable structure between traits and environments.
var <- c(0.085, 0.12, 0.06, 15.1, 8.5, 11.7)

# Dominance degree variances for trait 1 in 3 environments and for trait 2 in 3 environments,
# assuming a non-separable structure between traits and environments.
var_DD <- c(0.1, 0.1) # == c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)

# Additive genetic correlations between the two simulated traits.
T_cor_A <- matrix(
  c(
    1.0, 0.6,
    0.6, 1.0
  ),
  ncol = 2
)

# Additive genetic correlations between the three simulated environments.
E_cor_A <- matrix(
  c(
    1.0, 0.4, 0.6,
    0.4, 1.0, 0.5,
    0.6, 0.5, 1.0
  ),
  ncol = 3
)

cor_DD <- diag(6)

input_asr <- unstr_asr_input(
  n_envs = n_envs,
  n_traits = n_traits,
  mean = mean,
  var = var,
  T_cor_A = T_cor_A,
  E_cor_A = E_cor_A,
  mean_DD = mean_DD,
  var_DD = var_DD,
  cor_DD = cor_DD
)


# 2. Use input_asr to simulate genetic values in 'AlphaSimR' based on an unstructured model for
# GxE interaction.

library("AlphaSimR")
founders <- runMacs( # Simulation of founder genotypes using AlphaSimR's "MAIZE" presets
  nInd = n_ind, # to mimic the species' evolutionary history.
  nChr = n_chr,
  segSites = n_seg_sites,
  inbred = FALSE,
  species = "MAIZE"
)

SP <- SimParam$new(founders)

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

pool_A <- makeDH(founders[1:10], nDH = 1) # Pool A: 1 DH line from founders 1 to 10, respectively.
pool_B <- makeDH(founders[11:20], nDH = 1) # Pool B: 1 DH line from founders 11 to 20, respectively.

dh_lines <- mergePops(list(pool_A, pool_B))

factorial_plan <- as.matrix(expand.grid(A = pool_A@id, B = pool_B@id)) # Factorial crossing plan.

hybrid_pop <- makeCross(pop = dh_lines, crossPlan = factorial_plan, nProgeny = 1) # Hybrid genotypes.

# 3. Create a data frame containing the simulated genetic values for the two traits
#  in the three environments.

df_gv_unstr <- unstr_asr_output(
  pop = hybrid_pop,
  n_envs = n_envs,
  n_reps = n_reps,
  n_traits = n_traits
)



# Simulation of plot-level errors for two traits in three environments using a bivariate
# interpolation model for spatial variation.

rep_dir <- "row" # Layout of replicates (“above-and-below”).
n_cols <- 10 # Total umber of columns per location.
n_rows <- c(30, 30, 20) # Total number of rows per location.
plot_length <- 5 # Plot length; here in meters (column direction).
plot_width <- 2 # Plot width; here in meters (row direction).


H2 <- c(0.1, 0.1, 0.1, 0.3, 0.3, 0.3) # c(Yld:E1, Yld:E2, Yld:E3, Pht:E1, Pht:E2, Pht:E3)

# Error variances for traits 1 and 2.
calc_var_R <- function(var, H2) {
  varR <- (var / H2) - var
  return(varR)
}

var_R <- calc_var_R(var, H2)

# Spatial error correlations between traits 1 and 2.
spatial_model <- "Bivariate" # Spatial error model.
prop_spatial <- 0.4 # Proportion of spatial to total error variance.
S_cor_R <- rand_cor_mat(n_traits, min_cor = 0, max_cor = 0.5, pos_def = TRUE)

prop_ext <- 0.2
ext_dir <- "row"
ext_row_cor <- runif(1, min = -0.9, max = -0.6)
E_cor_R <- rand_cor_mat(n_traits, min_cor = 0, max_cor = 0.5, pos_def = TRUE)


# Simulate field error using bivariate interpolation.
error_df <- field_trial_error(
  n_envs = n_envs,
  n_traits = n_traits,
  n_reps = n_reps,
  rep_dir = rep_dir,
  n_cols = n_cols,
  n_rows = n_rows,
  plot_length = plot_length,
  plot_width = plot_width,
  var_R = var_R,
  R_cor_R = NULL,
  spatial_model = spatial_model,
  prop_spatial = prop_spatial,
  S_cor_R = S_cor_R,
  prop_ext = prop_ext,
  ext_dir = ext_dir,
  ext_row_cor = ext_row_cor,
  E_cor_R = E_cor_R,
  return_effects = TRUE
)

# Combine genetic values and field error to obtain phenotypes.
pheno_df <- make_phenotypes(df_gv_unstr, df_error_bivar, randomise = TRUE)

# str(df_gv_unstr)
# str(df_error_bivar)
# head(pheno_df)
#
# plot_effects(pheno_df,
#   env = 2,
#   effect = "phe.Trait.2"
# )

usethis::use_data(df_gv_unstr, overwrite = TRUE)
usethis::use_data(df_error_bivar, overwrite = TRUE)
