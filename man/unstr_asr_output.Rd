% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unstructured_gxe.R
\name{unstr_asr_output}
\alias{unstr_asr_output}
\title{Genetic values based on an unstructured model for GxE interaction using 'AlphaSimR' - Simulated
genetic values}
\usage{
unstr_asr_output(pop, n_envs, n_reps, n_traits)
}
\arguments{
\item{pop}{An 'AlphaSimR' population object (\code{\link[AlphaSimR]{Pop-class}} or
\code{\link[AlphaSimR]{HybridPop-class}}) generated using \link[FieldSimR]{unstr_asr_input}.}

\item{n_envs}{Number of simulated environments (same as in \link[FieldSimR]{unstr_asr_input}).}

\item{n_reps}{A vector defining the number of complete replicates in each environment. If only
one value is provided and \code{n_traits > 1}, all environments will be assigned the same
number of replicates.}

\item{n_traits}{Number of simulated traits (same as in \link[FieldSimR]{unstr_asr_input}).}
}
\value{
A data frame containing the environment id, replicate number, genotype id, and the
  simulated genetic values for each trait.
}
\description{
Creates a data frame of correlated genetic values for multiple traits in multiple environments
based on an unstructured model for genotype-by-environment (GxE) interaction. This function
requires an 'AlphaSimR' population object generated using the \link[FieldSimR]{unstr_asr_input}
function.
}
\examples{
# Simulation of genetic values in 'AlphaSimR' for two additive + dominance traits tested in
# three environments based on an unstructured model for GxE interaction.

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
E_cor_A <- stats::cov2cor(matrix(
  c(
    0.5, 0.4, 0.6,
    0.4, 1.0, 0.5,
    0.6, 0.5, 1.5
  ),
  ncol = 3
))

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
  nChr = 6,
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


# 3. Create a data frame containing the simulated genetic values for the two traits in the
# three environments.

n_reps <- c(2, 3, 2) # Vector containing the number of complete replicates in each environment.

gv_df <- unstr_asr_output(
  pop = pop,
  n_envs = 3,
  n_reps = n_reps,
  n_traits = 2
)
}
