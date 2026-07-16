# Simulate a random diagonal variance matrix

Creates a diagonal `n x n` variance matrix with user-defined minimum and
maximum variances based on a continuous uniform distribution.

## Usage

``` r
rand_diag_mat(n = 5, min.var = 0, max.var = 1)
```

## Arguments

- n:

  A scalar defining the dimensions of the variance matrix.

- min.var:

  A scalar defining the minimum variance.

- max.var:

  A scalar defining the maximum variance.  
  **Note:** `0 < min.var < max.var`.

## Value

A diagonal `n x n` variance matrix.

## Examples

``` r
# Simulate a random diagonal matrix with 10 columns and rows.
diag_mat <- rand_diag_mat(
  n = 10,
  min.var = 0,
  max.var = 0.2
)
```
