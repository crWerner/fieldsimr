# Simulate genetic values based on a compound symmetry model for GxE interaction - Simulation with \`AlphaSimR\`

Creates a data frame of simulated genetic values in multiple
environments for one or more traits based on a compound symmetry model
for genotype-by-environment (GxE) interaction. The wrapper function
`compsym_asr_output()` requires an
[\`AlphaSimR\`](https://CRAN.R-project.org/package=AlphaSimR) population
object generated with
[compsym_asr_input](https://crwerner.github.io/fieldsimr/reference/compsym_asr_input.md).

## Usage

``` r
compsym_asr_output(pop, ntraits = 1, nenvs, nreps = 1, return.effects = FALSE)
```

## Arguments

- pop:

  An [\`AlphaSimR\`](https://CRAN.R-project.org/package=AlphaSimR)
  population object
  ([Pop-class](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.html)
  or
  [HybridPop-class](https://gaynorr.github.io/AlphaSimR/reference/HybridPop-class.html))
  generated with
  [compsym_asr_input](https://crwerner.github.io/fieldsimr/reference/compsym_asr_input.md).

- ntraits:

  Number of traits specified in
  [compsym_asr_input](https://crwerner.github.io/fieldsimr/reference/compsym_asr_input.md).

- nenvs:

  Number of environments specified in
  [compsym_asr_input](https://crwerner.github.io/fieldsimr/reference/compsym_asr_input.md).

- nreps:

  A vector defining the number of replicates in each environment. If
  only one value is specified, all environments will be assigned the
  same number.

- return.effects:

  When `TRUE` (default is `FALSE`), a list is returned with additional
  entries containing the genotype main effects and GxE interaction
  effects for each trait.

## Value

A data frame with columns 'env', genotype 'id', and 'rep', followed by
the simulated genetic values for each trait. When
`return.effects = TRUE`, a list is returned with additional entries
containing the genotype main effects and GxE interaction effects for
each trait.

## Examples

``` r
# Simulate genetic values with 'AlphaSimR' for two additive + dominance traits
# in two environments based on a compound symmetry model.

# 1. Define the genetic architecture of the simulated traits.
# Mean genetic values and mean dominance degrees.
mean <- c(4.9, 5.4, 235.2, 228.5) # Trait 1 x 2 environments, Trait 2 x 2 environments
meanDD <- c(0.4, 0.4, 0.1, 0.1) # Trait 1 and 2, same value for both environments

# Additive genetic variances and dominance degree variances.
var <- c(0.08, 13) # Different values for Traits 1 and 2
varDD <- 0.2 # Same value for Traits 1 and 2

# Proportion of additive and dominance degree main effect variances.
prop.main <- c(0.4, 0.6) # Different values for Traits 1 and 2
prop.mainDD <- 0.4 # Same value for Traits 1 and 2

# Additive and dominance degree correlations between the two simulated traits.
corA <- matrix(c(
  1.0, 0.5,
  0.5, 1.0
), ncol = 2)
corDD <- diag(2) # Assuming independence

input_asr <- compsym_asr_input(
  ntraits = 2,
  nenvs = 2,
  mean = mean,
  var = var,
  prop.main = prop.main,
  corA = corA,
  meanDD = meanDD,
  varDD = varDD,
  prop.mainDD = prop.mainDD,
  corDD = corDD
)


# 2. Use input_asr to simulate genetic values with 'AlphaSimR' based on a
# compound symmetry model.

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

gv_ls <- compsym_asr_output(
  pop = pop,
  ntraits = 2,
  nenvs = 2,
  nreps = 2,
  return.effects = TRUE
)
#> Error: object 'pop' not found
```
