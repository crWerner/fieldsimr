# Simulate genetic values based on an unstructured model for GxE interaction - \`AlphaSimR\` input parameters

Creates a list of input parameters for
[\`AlphaSimR\`](https://CRAN.R-project.org/package=AlphaSimR) to
simulate genetic values in multiple environments for one or more traits
based on an unstructured model for genotype-by-environment (GxE)
interaction.  
This function utilises the ability of \`AlphaSimR\` to simulate
correlated traits. The wrapper function `unstr_asr_input()` is used to
specify the input parameters required in \`AlphaSimR\`, and can handle
separable and non-separable structures between traits and environments
(see below). After simulating the genetic values, the wrapper function
[unstr_asr_output](https://crwerner.github.io/fieldsimr/reference/unstr_asr_output.md)
can be used to generate a data frame with output values.

## Usage

``` r
unstr_asr_input(
  ntraits = 1,
  nenvs = 2,
  mean = 0,
  var = NULL,
  Tvar = NULL,
  Evar = NULL,
  corA = NULL,
  TcorA = NULL,
  EcorA = NULL,
  meanDD = NULL,
  varDD = NULL,
  TvarDD = NULL,
  EvarDD = NULL,
  corDD = NULL,
  TcorDD = NULL,
  EcorDD = NULL,
  relAA = NULL,
  corAA = NULL,
  TcorAA = NULL,
  EcorAA = NULL
)
```

## Arguments

- ntraits:

  Number of traits to be simulated.

- nenvs:

  Number of environments to be simulated (minimum of two).

- mean:

  A vector of mean genetic values for each environment-within-trait
  combination. If only one value is specified, all combinations will be
  assigned the same mean.

- var:

  A vector of genetic variances for each environment-within-trait
  combination. If only one value is specified, all combinations will be
  assigned the same variance.  
  **Alternatively**, if a separable structure between traits and
  environments is desired, `Tvar` and `Evar` can be specified.

- Tvar:

  A vector of genetic variances for each trait. Must be provided in
  combination with `Evar`.  
  **Alternatively**, `var` can be specified.

- Evar:

  A vector of genetic variances for each environment. Must be provided
  in combination with `Tvar`.  
  **Alternatively**, `var` can be specified.

- corA:

  A matrix of additive genetic correlations between
  environment-within-trait combinations. By default, a diagonal matrix
  is constructed.  
  **Alternatively**, `TcorA` and `EcorA` can be specified.

- TcorA:

  A matrix of additive genetic correlations between traits. Must be
  provided in combination with `EcorA`.  
  **Alternatively**, `corA` can be specified.

- EcorA:

  A matrix of additive genetic correlations between environments. Must
  be provided in combination with `TcorA`.  
  **Alternatively**, `corA` can be specified.

- meanDD:

  A vector of mean dominance degrees for each environment-within-trait
  combination (similar to `mean`). If only one value is specified, all
  combinations will be assigned the same mean. By default,
  `meanDD = NULL` and dominance is not simulated.

- varDD:

  A vector of dominance degree variances for each
  environment-within-trait combination (similar to `var`). If only one
  value is specified, all combinations will be assigned the same
  variance.  
  **Alternatively**, if a separable structure between traits and
  environments is desired, `TvarDD` and `EvarDD` can be specified.

- TvarDD:

  A vector of dominance degree variances for each trait (similar to
  `Tvar`). Must be provided in combination with `EvarDD`.  
  **Alternatively**, `varDD` can be specified.

- EvarDD:

  A vector of dominance degree variances for each environment (similar
  to `Evar`). Must be provided in combination with `TvarDD`.  
  **Alternatively**, `varDD` can be specified.

- corDD:

  A matrix of dominance degree correlations between
  environment-within-trait combinations (similar to `corA`). If not
  specified and dominance is simulated, a diagonal matrix is
  constructed.  
  **Alternatively**, `TcorDD` and `EcorDD` can be specified.

- TcorDD:

  A matrix of dominance degree correlations between traits (similar to
  `TcorA`). Must be provided in combination with `EcorDD`.  
  **Alternatively**, `corDD` can be specified.

- EcorDD:

  A matrix of dominance degree correlations between environments
  (similar to `EcorA`). Must be provided in combination with `TcorDD`.  
  **Alternatively**, `corDD` can be specified.

- relAA:

  A vector defining the relative magnitude of additive-by-additive
  (epistatic) variance to additive genetic variance for each
  environment-within-trait combination, that is in a diploid organism
  with allele frequency of 0.5. If only one value is specified, all
  environment-within-trait combinations will be assigned the same value.
  By default, `relAA = NULL` and epistasis is not simulated.  
  **Alternatively**, if a separable structure between traits and
  environments is desired, `TrelAA` and `ErelAA` can be specified.

- corAA:

  A matrix of epistatic correlations between environment-within-trait
  combinations (similar to `corA`). If not specified and epistasis is
  simulated, a diagonal matrix is constructed.  
  **Alternatively**, `TcorAA` and `EcorAA` can be specified.

- TcorAA:

  A matrix of epistatic correlations between traits (similar to
  `TcorA`). Must be provided in combination with `EcorAA`.  
  **Alternatively**, `corAA` can be specified.

- EcorAA:

  A matrix of epistatic correlations between environments (similar to
  `EcorA`). Must be provided in combination with `TcorAA`.  
  **Alternatively**, `corAA` can be specified.

## Value

A list with input parameters for \`AlphaSimR\`, which are used to
simulate correlated genetic values based on an unstructured model for
GxE interaction.

## Details

`unstr_asr_input` can handle separable and non-separable structures
between traits and environments.

- For separable structures, provide (1) `Tvar` & `Evar`, and (2) `TcorA`
  & `EcorA`.

- For non-separable structures, provide (1) `var`, and (2) `corA`.  

**Note:** \`AlphaSimR\` can simulate different biological effects (see:
[SimParam](https://gaynorr.github.io/AlphaSimR/reference/SimParam.html)).

- For additive traits use `addTraitA()`.

- For additive + dominance traits use `addTraitAD()`.

- For additive + epistatic traits use `addTraitAE()`.

- For additive + dominance + epistatic traits use `addTraitADE()`.

Check the `useVarA` argument of these functions when simulating
non-additive traits.

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
```
