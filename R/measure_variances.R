#' Measure variances from a covariance matrix
#'
#' Creates a table of variance components derived from a user-defined covariance matrix.
#'
#' @param mat A symmetric \code{n x n} positive (semi)-definite variance matrix.
#' @param correction When \code{TRUE} (default is \code{FALSE}), the variance components are treated
#'   as estimates and sample corrections are applied.
#'
#' @return A table which partitions the total variance into main effect and interaction variance,
#'   heterogeneity of variance and lack of correlation, and non-crossover and crossover interaction.
#'   When \code{correction = TRUE}, a sample correction is applied to the variance components.
#'
#' @examples
#' # Generate a structured covariance matrix and then measure the partitioning of variance.
#'
#' diag_mat <- rand_diag_mat(
#'   n = 10
#' )
#' cor_mat <- struc_cor_mat(
#'   n = 10,
#' )
#'
#' cov_mat <- sqrt(diag_mat) %*% cor_mat %*% sqrt(diag_mat)
#'
#' measure_variances(
#'   mat = cov_mat
#' )
#'
#' @export
measure_variances <- function(mat,
                              correction = FALSE) {
  if (!is.matrix(mat)) stop("'mat' must be a matrix")
  if (!isSymmetric(mat)) stop("'mat' must be a symmetric matrix")
  if (any(diag(mat) < 0)) stop("All diagonal elements of 'mat' must >= 0")
  if (any(is.na(mat))) stop("'mat' must not contain missing values")
  if (!is.logical(correction)) stop("'correction' must be logical")

  n <- ncol(mat)
  total_var <- mean(diag(mat))

  if (correction) {
    main_eff <- mean(mat[upper.tri(mat)])
    if (main_eff < 0) stop("The mean upper triangular element of 'mat' must be >= 0")
    if (main_eff < 1e-8) warning("The main effect variance is zero")
    int_eff <- mean(diag(mat - main_eff))
    if (round(main_eff + int_eff, 8) != round(total_var, 8)) stop("Sum of main effect and interaction variances does not match total variance")

    sqrt_dg <- sqrt(diag(mat))
    het_scale <- sum((sqrt_dg - mean(sqrt_dg))^2) / (n - 1)
    lack_cor <- sum(outer(sqrt_dg, sqrt_dg, "*") * (1 - stats::cov2cor(mat))) / n / (n - 1)
    if (round(het_scale + lack_cor, 8) != round(int_eff, 8)) stop("Sum of heterogeneity of scale and lack of correlation does not match interaction variance")
    df <- data.frame(
      Component = c(
        "Main effect", "Interaction",
        "Heterogeniety of scale", "Lack of correlation",
        "Total"
      ),
      Variance = c(main_eff, int_eff, het_scale, lack_cor, total_var)
    )
  } else if (!correction) {
    main_eff <- mean(mat)
    if (main_eff < 0) stop("The mean element of 'mat' must be >= 0")
    if (main_eff < 1e-8) warning("The main effect variance is zero")
    int_eff <- mean(diag(mat - main_eff))
    if (round(main_eff + int_eff, 8) != round(total_var, 8)) stop("Sum of main effect and interaction variances does not match total variance")

    sqrt_dg <- sqrt(diag(mat))
    het_scale <- sum((sqrt_dg - mean(sqrt_dg))^2) / n
    lack_cor <- sum(outer(sqrt_dg, sqrt_dg, "*") * (1 - stats::cov2cor(mat))) / n^2
    if (round(het_scale + lack_cor, 8) != round(int_eff, 8)) stop("Sum of heterogeneity of scale and lack of correlation does not match interaction variance")

    cov_mat <- rowMeans(mat)
    adjust <- min(0, cov_mat)
    if (adjust < 0) {
      warning("Non-crossover and crossover interaction not completely independent, covariances will be attributed to crossover variance")
    }
    main_eff_adjust <- main_eff - adjust
    cov_mat_adjust <- cov_mat - adjust
    non_cross <- main_eff * t(cov_mat_adjust) %*% cov_mat_adjust / n / (main_eff_adjust^2)
    cross <- total_var - main_eff * t(cov_mat_adjust) %*% cov_mat_adjust / n / (main_eff_adjust^2)
    if (round(non_cross + cross, 8) != round(total_var, 8)) stop("Sum of non-crossover and crossover variance does not match total variance")

    df <- data.frame(
      Component = c(
        "Main effect", "Interaction",
        "Heterogeniety of scale", "Lack of correlation",
        "Non-crossover", "Crossover",
        "Total"
      ),
      Variance = c(main_eff, int_eff, het_scale, lack_cor, non_cross, cross, total_var)
    )
  }
  df$Proportion <- df$Variance / total_var
  return(df)
}
