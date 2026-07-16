# Measure variances from a covariance matrix

Creates a table of variance components derived from a user-defined
covariance matrix.

## Usage

``` r
measure_variances(mat, correction = FALSE)
```

## Arguments

- mat:

  A symmetric `n x n` positive (semi)-definite variance matrix.

- correction:

  When `TRUE` (default is `FALSE`), the variance components are treated
  as estimates and sample corrections are applied.

## Value

A table which partitions the total variance into main effect and
interaction variance, heterogeneity of variance and lack of correlation,
and non-crossover and crossover interaction. When `correction = TRUE`, a
sample correction is applied to the variance components.

## Examples

``` r
# Generate a structured covariance matrix and then measure the partitioning of variance.

diag_mat <- rand_diag_mat(
  n = 10
)
cor_mat <- struc_cor_mat(
  n = 10,
)

cov_mat <- sqrt(diag_mat) %*% cor_mat %*% sqrt(diag_mat)

measure_variances(
  mat = cov_mat
)
#>                Component   Variance Proportion
#> 1            Main effect 0.24639625 0.50650027
#> 2            Interaction 0.24007189 0.49349973
#> 3 Heterogeniety of scale 0.03018576 0.06205084
#> 4    Lack of correlation 0.20988614 0.43144889
#> 5          Non-crossover 0.26990818 0.55483218
#> 6              Crossover 0.21655996 0.44516782
#> 7                  Total 0.48646814 1.00000000
```
