# Simulate genetic values based on a multiplicative model for GxE interaction - Simulation with \`AlphaSimR\`

Creates a data frame of simulated genetic values in multiple
environments for one or more traits based on a (reduced rank)
multiplicative model for genotype-by-environment (GxE) interaction. This
function requires an
[\`AlphaSimR\`](https://CRAN.R-project.org/package=AlphaSimR) population
object generated with
[multi_asr_input](https://crwerner.github.io/fieldsimr/reference/multi_asr_input.md).

## Usage

``` r
multi_asr_output(
  pop,
  ntraits = 1,
  nenvs,
  nreps = 1,
  cov.mat,
  return.effects = FALSE
)
```

## Arguments

- pop:

  An [\`AlphaSimR\`](https://CRAN.R-project.org/package=AlphaSimR)
  population object
  ([Pop-class](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.html)
  or
  [HybridPop-class](https://gaynorr.github.io/AlphaSimR/reference/HybridPop-class.html))
  generated with
  [multi_asr_input](https://crwerner.github.io/fieldsimr/reference/multi_asr_input.md).

- ntraits:

  Number of traits specified in
  [multi_asr_input](https://crwerner.github.io/fieldsimr/reference/multi_asr_input.md).

- nenvs:

  Number of environments specified in
  [multi_asr_input](https://crwerner.github.io/fieldsimr/reference/multi_asr_input.md).

- nreps:

  A vector defining the number of replicates in each environment. If
  only one value is specified, all environments will be assigned the
  same number.

- cov.mat:

  A matrix of covariates that will be used to construct the genetic
  values, typically generated with
  [multi_asr_input](https://crwerner.github.io/fieldsimr/reference/multi_asr_input.md).

- return.effects:

  When `TRUE` (default is `FALSE`), a list is returned with additional
  entries containing the genotype slopes for each multiplicative term.

## Value

A data frame with columns 'env', genotype 'id', and 'rep', followed by
the simulated genetic values for each trait. When
`return.effects = TRUE`, a list is returned with additional entries
containing the genotype slopes for each multiplicative term.

## Examples

``` r
# Simulate genetic values with 'AlphaSimR' for two additive traits in two
# environments based on a multiplicative model with three terms.

# 1. Define the genetic architecture of the simulated traits.
# Mean genetic values.
mean <- c(5, 240) # Trait 1, Trait 2

# Additive genetic variances.
var <- c(0.086, 0.12, 15.1, 8.5) # Trait 1 x 2 environments, Trait 2 x 2 environments

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

# Construct separable additive genetic correlation matrix
corA <- kronecker(TcorA, EcorA)

input_asr <- multi_asr_input(
  ntraits = 2,
  nenvs = 2,
  mean = mean,
  var = var,
  corA = corA,
  nterms = 3
)
#> Warning message: 
#>  'nterms' is less than rank of 'corA', 99.79% of variation captured with 3 terms


# 2. Use input_asr to simulate genetic values in 'AlphaSimR' based on a
# multiplicative model with three terms.

library("AlphaSimR")
FOUNDERPOP <- quickHaplo(
  nInd = 10,
  nChr = 1,
  segSites = 20
)

SP <- SimParam$new(FOUNDERPOP)


SP$addTraitA(
  nQtlPerChr = 20,
  mean = input_asr$mean,
  var = input_asr$var,
  corA = input_asr$corA
)

pop <- newPop(FOUNDERPOP)
#> Error in get("SP", envir = .GlobalEnv): object 'SP' not found


# 3. Create a data frame with simulated genetic values for the two traits in the two
# environments, with two replicates of each genotype.

# The covariates are obtained from input_asr.

gv_ls <- multi_asr_output(
  pop = pop,
  ntraits = 2,
  nenvs = 2,
  nreps = 2,
  cov.mat = input_asr$cov.mat,
  return.effects = TRUE
)
#> Error: object 'pop' not found
```
