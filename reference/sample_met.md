# Sample environments from a target population

Creates a list of environments sampled from a population with
user-defined sample size.

## Usage

``` r
sample_met(
  ntraits = 1,
  nenvs = 1000,
  nsamples = 10,
  sample.size = 20,
  replace = TRUE,
  cov.mat = NULL
)
```

## Arguments

- ntraits:

  A scalar defining the number of traits.

- nenvs:

  A scalar defining the number of environments in the target population.

- nsamples:

  A scalar defining the number of samples to be taken.

- sample.size:

  A vector defining the number of environments in each sample. When only
  one value is specified, all samples will be assigned the same number.

- replace:

  When `TRUE` (default), samples are taken with replacement. Ignored
  when `nsamples = 1`.

- cov.mat:

  An optional matrix of environmental covariates for one or more traits.
  When supplied, the covariates are sampled and printed.

## Value

A list with elements given by the sample of environments taken from the
target population. When `cov.mat` is supplied, additional entries are
given containing the sampled environmental covariates for each trait.

## Examples

``` r
# Sample environments from a target population of 1000, with each sample containing 20 environments.
cov_ls <- sample_met(
  nenvs = 1000,
  nsamples = 10,
  sample.size = 20,
  replace = TRUE
)
```
