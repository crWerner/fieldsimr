# Simulate genetic values based on a multiplicative model for GxE interaction - \`AlphaSimR\` input parameters

Creates a list of input parameters for
[\`AlphaSimR\`](https://CRAN.R-project.org/package=AlphaSimR) to
simulate genetic values in multiple environments for one or more traits
based on a (reduced rank) multiplicative model for
genotype-by-environment (GxE) interaction.  
This function utilises the ability of \`AlphaSimR\` to simulate
correlated traits. The wrapper function `multi_asr_input()` is used to
specify the input parameters required in \`AlphaSimR\`. After simulating
the genetic values, the wrapper function
[multi_asr_output](https://crwerner.github.io/fieldsimr/reference/multi_asr_output.md)
can be used to generate a data frame with output values.

## Usage

``` r
multi_asr_input(
  ntraits = 1,
  nenvs = 2,
  mean = 0,
  var = 1,
  corA = NULL,
  nterms = NULL
)
```

## Arguments

- ntraits:

  Number of traits to be simulated.

- nenvs:

  Number of environments to be simulated (minimum of two).

- mean:

  A vector of mean genetic values for each trait or each
  environment-within-trait combination. If only one value is specified,
  all combinations will be assigned the same mean.

- var:

  A vector of additive genetic variances for each trait or each
  environment-within-trait combination. If only one value is specified,
  all combinations will be assigned the same variance.

- corA:

  A matrix of additive genetic correlations between
  environment-within-trait combinations. By default, a diagonal matrix
  is constructed.

- nterms:

  A scalar defining the number of multiplicative terms to be simulated.
  By default, the number of terms is set to the number of
  environment-within-trait combinations. **Note:** when `nterms` is less
  than the number of environment-within-trait combinations, the values
  in `mean` will be approximated.

## Value

A list with input parameters for \`AlphaSimR\`, which are used to
simulate correlated genetic values based on a multiplicative model for
GxE interaction. Covariates are also supplied for use in
[multi_asr_output](https://crwerner.github.io/fieldsimr/reference/multi_asr_output.md).

## Details

Currently supports additive traits only, but other (non-additive) traits
are being implemented.

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

# Construct separable additive genetic correlation matrix.
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
```
