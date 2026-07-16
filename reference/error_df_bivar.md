# Plot errors - Example data frame

An example data frame with simulated plot errors for two traits in three
environments. Environments 1 and 2 comprise two blocks, while
Environment 3 comprises three blocks. The blocks are aligned in the
column direction (side-by-side) and comprise 5 columns and 20 rows. The
data frame was generated using the function
[field_trial_error](https://crwerner.github.io/fieldsimr/reference/field_trial_error.md)
with bivariate interpolation. The simulation is demonstrated in the
vignette [Simulation of plot errors and phenotypes in a plant breeding
field
trial](https://crwerner.github.io/fieldsimr/articles/spatial_error_demo.html).

## Usage

``` r
error_df_bivar
```

## Format

A data frame with 700 rows and 6 columns:

- env:

  Environment number

- block:

  Block number

- col:

  Column number

- row:

  Row number

- e.Trait1:

  Simulated plot errors for Trait 1

- e.Trait2:

  Simulated plot errors for Trait 2
