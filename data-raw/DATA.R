# Simulation of genetic values in 'AlphaSimR' for two additive + dominance traits tested in
# three environments based on an unstructured model for GxE interaction.

set.seed(1312)
library(FieldSimR)

# 1. Define the genetic architecture of the simulated traits.
# Mean genetic values and mean dominance degrees for trait 1 in all 3 environments and trait 2
# in all 3 environments.
mean <- c(1, 3, 2, 80, 70, 100) # Trait 1 x 3 environments, trait 2 x 3 environments.
mean_DD <- c(0.1, 0.4) # Trait 1 and 2, same values set in all 3 environments for each trait.

# Additive genetic variances (useVarA = TRUE) and dominance degree variances for traits 1 and 2,
# assuming a separable structure between traits and environments.
T_var <- c(0.2, 10) # Genetic variances defined for the two traits.
E_var <- c(0.5, 1, 1.5) # Genetic variances defined for the three environments.

# Dominance degree variances for trait 1 in 3 environments and for trait 2 in 3 environments,
# assuming a non-separable structure between traits and environments.
var_DD <- c(0.1, 0.15, 0.2, 0.2, 0.3, 0.2)

# Additive genetic correlations between the two simulated traits.
T_cor_A <- matrix(
  c(
    1.0, 0.3,
    0.3, 1.0
  ),
  ncol = 2
)

# Additive genetic correlations between the three simulated environments.
E_cor_A <- matrix(
  c(
    1, 0.2, 0.45,
    0.2, 1, -0.15,
    0.45, -0.15, 1
  ),
  ncol = 3
)

# Dominance degree correlation between all six trait-by-environment combinations.
cor_DD <- diag(6) # Assuming independence between traits

input_asr <- unstr_asr_input(
  n_envs = 3,
  n_traits = 2,
  mean = mean,
  T_var = T_var,
  E_var = E_var,
  T_cor_A = T_cor_A,
  E_cor_A = E_cor_A,
  mean_DD = mean_DD,
  var_DD = var_DD,
  cor_DD = cor_DD
)


# 2. Use input_asr to simulate genetic values in 'AlphaSimR' based on an unstructured model for
# GxE interaction.

library("AlphaSimR")
FOUNDERPOP <- quickHaplo(
  nInd = 100,
  nChr = 7,
  segSites = 100
)

SP <- SimParam$new(FOUNDERPOP)

SP$addTraitAD(
  nQtlPerChr = 100,
  mean = input_asr$mean,
  var = input_asr$var,
  meanDD = input_asr$mean_DD,
  varDD = input_asr$var_DD,
  corA = input_asr$cor_A,
  corDD = input_asr$cor_DD,
  useVarA = TRUE
)

# By default, the value provided in 'var' represents the additive variance.
# If useVarA=FALSE, 'var' represents the total genetic variance.

pop <- newPop(FOUNDERPOP)


# 3. Create a data frame containing the simulated genetic values for the two traits
#  in the three environments.

n_envs <- 3
n_traits <- 2
n_reps <- c(3, 3, 2) # Vector containing the number of complete replicates in each
# environment.

df_gv_unstr <- unstr_asr_output(
  pop = pop,
  n_envs = n_envs,
  n_reps = n_reps,
  n_traits = n_traits
)


# Simulation of plot-level errors for two traits in three environments using a bivariate
# interpolation model for spatial variation.

# Field layout
n_cols <- 10 # Total number of columns in each environment.
n_rows <- c(30, 30, 20) # Total number of rows in each environment.
plot_length <- 5 # Plot length set to 5 meters in each environment.
plot_width <- 2 # Plot width set to 2 meters in each environment.

# Error variances for traits 1 and 2.
var_R <- c(0.4, 15)

# Spatial error correlations between traits 1 and 2.
S_cor_R <- matrix(
  c(
    1.0, 0.2,
    0.2, 1.0
  ),
  ncol = 2
)

# Simulate field error using bivariate interpolation.
df_error_bivar <- field_trial_error(
  n_envs = n_envs,
  n_traits = n_traits,
  n_reps = n_reps,
  n_cols = n_cols,
  n_rows = n_rows,
  plot_length = plot_length,
  plot_width = plot_width,
  rep_dir = "row",
  var_R = var_R,
  S_cor_R = S_cor_R,
  spatial_model = "bivariate",
  complexity = NULL,
  prop_spatial = 0.5,
  prop_ext = 0.2,
  ext_dir = "row",
  return_effects = FALSE
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
