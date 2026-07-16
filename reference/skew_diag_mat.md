# Simulate a skewed diagonal variance matrix

Creates a diagonal `n x n` variance matrix with user-defined skewness
based on a gamma or inverse gamma distribution.

## Usage

``` r
skew_diag_mat(n = 5, shape = 1.5, scale = 1, inverse = FALSE, mean.var = NULL)
```

## Arguments

- n:

  A scalar defining the dimensions of the variance matrix.

- shape:

  A scalar defining the shape of the distribution.

- scale:

  A scalar defining the scale of the distribution.

- inverse:

  When `TRUE` (default is `FALSE`), the variances are sampled from the
  inverse gamma distribution instead of the gamma distribution.

- mean.var:

  An optional scalar defining the mean variance. . When supplied, the
  variances are scaled to achieve the defined mean.

## Value

A diagonal `n x n` variance matrix.

## Examples

``` r
# Simulate a random diagonal matrix with 10 columns and rows, and negatively skewed variances
# scaled to a mean of 0.1.
diag_mat <- skew_diag_mat(
  n = 10,
  shape = 1.5,
  scale = 1,
  mean.var = 0.1
)
```
