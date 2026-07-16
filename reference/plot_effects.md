# Graphics for plot effects

Creates a graphical field array for a set of plot effects (e.g.,
phenotypes, genetic values, or plot errors). Requires a data frame
generated with the functions
[field_trial_error](https://crwerner.github.io/fieldsimr/reference/field_trial_error.md)
or
[make_phenotypes](https://crwerner.github.io/fieldsimr/reference/make_phenotypes.md),
or any data frame matching the description below.

## Usage

``` r
plot_effects(df, effect, blocks = TRUE, labels = TRUE, dim.names = NULL)
```

## Arguments

- df:

  A data frame with columns containing the column and row dimensions,
  and the effects to be plotted.

- effect:

  A character defining the effects to be plotted.

- blocks:

  When `TRUE` (default), the field array is split into blocks. This
  requires an additional column in the data frame, as specified by
  `dim.names`.

- labels:

  When `TRUE` (default), column and row labels are displayed.

- dim.names:

  An optional vector defining the column, row and block dimensions
  ('col', 'row' and 'block' by default).

## Value

A graphical field array with x- and y-axes displaying the column and row
numbers, and colour gradient ranging from red (low value) to green (high
value).

## Examples

``` r
# Display the simulated plot errors in the example data frame 'error_df_bivar'
# for Trait 1 in Environment 1.

error_df <- error_df_bivar[error_df_bivar$env == 1, ]

plot_effects(
  df = error_df,
  effect = "e.Trait1",
  labels = TRUE,
)

```
