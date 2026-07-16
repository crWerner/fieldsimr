# Simulate a structured correlation matrix with reduced rank

Creates a symmetric `n x n` correlation matrix with user-defined
structure and rank.

## Usage

``` r
struc_cor_mat(
  n = 5,
  base.cor = 0.5,
  range = NULL,
  rank = 3,
  skew = 0,
  base.mat = NULL,
  pos.def = FALSE,
  small.positive = NULL
)
```

## Arguments

- n:

  A scalar defining the dimensions of the correlation matrix.

- base.cor:

  A scalar defining the baseline correlation.  
  **Note:** `-1 < base.cor < 1`.

- range:

  A scalar defining the range of correlations around the baseline. By
  default, `range = 1 - base.cor` which ensures the matrix is positive
  semi-definite with defined rank. **Note:** `base.cor + range <= 1`.

- rank:

  A scalar defining the rank of the correlation matrix.

- skew:

  A scalar defining the skewness imposed on the correlations. **Note:**
  `-1 < skew < 1`.

- base.mat:

  An optional `n x n` base correlation matrix. When supplied, `base.cor`
  and `skew` are ignored and noise is simulated based on `rank`.

- pos.def:

  When `TRUE` (default is `FALSE`), the function `bend` of the R package
  [\`mbend\`](https://cran.r-project.org/package=mbend) is used to bend
  a non-positive (semi)-definite matrix to be positive (semi)-definite.

- small.positive:

  Argument passed to `bend` when `pos.def = TRUE` (default is 1e-8).
  Eigenvalues smaller than `small.positive` are replaced by this.  
  **Note:** `0 < small.positive < 0.1`.

## Value

A symmetric `n x n` correlation matrix with defined rank. When
`pos.def = TRUE`, the correlation matrix is guaranteed to be positive
(semi)-definite.

## Examples

``` r
# Simulate and visualise a correlation matrix with 10 columns and rows, rank equal to 4 and
# negatively skewed correlations.
cor_mat <- struc_cor_mat(
  n = 10,
  base.cor = 0.3,
  range = 0.7,
  rank = 4,
  skew = -0.5
)

plot_matrix(
  mat = cor_mat,
  order = TRUE
)

```
