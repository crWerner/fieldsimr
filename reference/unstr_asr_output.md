# Simulate genetic values based on an unstructured model for GxE interaction - Simulation with \`AlphaSimR\`

Creates a data frame of simulated genetic values in multiple
environments for one or more traits based on an unstructured model for
genotype-by-environment (GxE) interaction. The wrapper function
`unstr_asr_output` requires an
[\`AlphaSimR\`](https://CRAN.R-project.org/package=AlphaSimR) population
object generated with
[unstr_asr_input](https://crwerner.github.io/fieldsimr/reference/unstr_asr_input.md).

## Usage

``` r
unstr_asr_output(pop, ntraits = 1, nenvs, nreps = 1)
```

## Arguments

- pop:

  An [\`AlphaSimR\`](https://CRAN.R-project.org/package=AlphaSimR)
  population object
  ([Pop-class](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.html)
  or
  [HybridPop-class](https://gaynorr.github.io/AlphaSimR/reference/HybridPop-class.html))
  generated with
  [unstr_asr_input](https://crwerner.github.io/fieldsimr/reference/unstr_asr_input.md).

- ntraits:

  Number of simulated traits specified in
  [unstr_asr_input](https://crwerner.github.io/fieldsimr/reference/unstr_asr_input.md).

- nenvs:

  Number of simulated environments specified in
  [unstr_asr_input](https://crwerner.github.io/fieldsimr/reference/unstr_asr_input.md).

- nreps:

  A vector defining the number of replicates in each environment. If
  only one value is specified, all environments will be assigned the
  same number.

## Value

A data frame with columns 'env', genotype 'id', and 'rep', followed by
the simulated genetic values for each trait.

## Examples

``` r
# Simulate genetic values with 'AlphaSimR' for two additive + dominance traits
# in two environments based on an unstructured model.

# 1. Define the genetic architecture of the simulated traits.
# Mean genetic values and mean dominance degrees.
mean <- c(4.9, 5.4, 235.2, 228.5) # Trait 1 x 2 environments, Trait 2 x 2 environments
meanDD <- c(0.4, 0.4, 0.1, 0.1) # Trait 1 and 2, same value for both environments

# Additive genetic variances and dominance degree variances.
var <- c(0.086, 0.12, 15.1, 8.5) # Trait 1 x 2 environments, Trait 2 x 2 environments
varDD <- 0.2 # Same value for all environment-within-trait combinations

# Additive genetic correlations between the two simulated traits.
TcorA <- matrix(c(
  1.0, 0.6,
  0.6, 1.0
), ncol = 2)

# Additive genetic correlations between the two simulated environments.
EcorA <- matrix(c(
  1.0, 0.2,
  0.2, 1.0
), ncol = 2)

# Dominance degree correlations between the four environment-within-trait combinations.
corDD <- diag(4) # Assuming independence

input_asr <- unstr_asr_input(
  ntraits = 2,
  nenvs = 2,
  mean = mean,
  var = var,
  TcorA = TcorA,
  EcorA = EcorA,
  meanDD = meanDD,
  varDD = varDD,
  corDD = corDD
)


# 2. Use input_asr to simulate genetic values with 'AlphaSimR' based on an
# unstructured model.

library("AlphaSimR")
FOUNDERPOP <- quickHaplo(
  nInd = 10,
  nChr = 1,
  segSites = 20
)

SP <- SimParam$new(FOUNDERPOP)


SP$addTraitAD(
  nQtlPerChr = 20,
  mean = input_asr$mean,
  var = input_asr$var,
  corA = input_asr$corA,
  meanDD = input_asr$meanDD,
  varDD = input_asr$varDD,
  corDD = input_asr$corDD,
  useVarA = TRUE
)

# By default, the variances in 'var' represent additive genetic variances.
# When useVarA = FALSE, the values represent total genetic variances.

pop <- newPop(FOUNDERPOP)
#> Error in get("SP", envir = .GlobalEnv): object 'SP' not found


# 3. Create a data frame with simulated genetic values for the two traits in
# the two environments, with two replicates of each genotype.

gv_df <- unstr_asr_output(
  pop = pop,
  ntraits = 2,
  nenvs = 2,
  nreps = 2
)
#> Error: object 'pop' not found
```
