# Genetic values - Example data frame

An example data frame with simulated genetic values of 100 genotypes for
two traits in three environments. Environments 1 and 2 comprise two
replicates of each genotype, while Environment 3 comprises three
replicates. The data frame was generated using the wrapper functions
[unstr_asr_input](https://crwerner.github.io/fieldsimr/reference/unstr_asr_input.md)
and
[unstr_asr_output](https://crwerner.github.io/fieldsimr/reference/unstr_asr_output.md),
which simulate correlated genetic values with
[\`AlphaSimR\`](https://CRAN.R-project.org/package=AlphaSimR). The
simulation is demonstrated in the vignette [Simulation of genetic values
based on an unstructured model for GxE
interaction](https://crwerner.github.io/fieldsimr/articles/unstructured_GxE_demo.html).

## Usage

``` r
gv_df_unstr
```

## Format

A data frame with 700 rows and 5 columns:

- env:

  Environment number

- id:

  Genotype identifier

- rep:

  Replicate number

- gv.Trait1:

  Simulated genetic values for Trait 1

- gv.Trait2:

  Simulated genetic values for Trait 2
