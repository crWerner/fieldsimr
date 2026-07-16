# Simulate a random correlation matrix

Creates a symmetric `n x n` correlation matrix with user-defined minimum
and maximum correlations based on a continuous uniform distribution.

## Usage

``` r
rand_cor_mat(
  n = 5,
  min.cor = -1,
  max.cor = 1,
  pos.def = FALSE,
  small.positive = NULL
)
```

## Arguments

- n:

  A scalar defining the dimensions of the correlation matrix.

- min.cor:

  A scalar defining the minimum correlation.

- max.cor:

  A scalar defining the maximum correlation.  
  **Note:** `-1 < min.cor < max.cor < 1`.

- pos.def:

  When `TRUE` (default is `FALSE`), the function `bend` of the R package
  [\`mbend\`](https://cran.r-project.org/package=mbend) is used to bend
  a non-positive (semi)-definite matrix to be positive (semi)-definite.

- small.positive:

  Argument passed to `bend` when `pos.def = TRUE` (default is 1e-8).
  Eigenvalues smaller than `small.positive` are replaced by this.  
  **Note:** `0 < small.positive < 0.1`.

## Value

A symmetric `n x n` correlation matrix. When `pos.def = TRUE`, the
correlation matrix is guaranteed to be positive (semi)-definite.

## Examples

``` r
# Simulate and visualise a random correlation matrix with 10 columns and rows.
cor_mat <- rand_cor_mat(
  n = 10,
  min.cor = -0.2,
  max.cor = 0.8,
  pos.def = TRUE
)
#> Unweighted bending
#> max.iter = 10000
#> small.positive = 1e-08
#> method = hj
#> Found a correlation matrix.
#> Convergence met after 17 iterations.

plot_matrix(
  mat = cor_mat,
  order = TRUE
)

```
