#' Random correlation matrix
#'
#' Creates a general \code{p x p} correlation matrix with user-defined maximum and minimum
#' correlations. If the randomly generated correlation matrix is not positive definite,
#' the function \link[mbend]{bend} of the package 'mbend' is used with default arguments to turn
#' the correlation matrix into a positive definite correlation matrix.
#'
#' @param p A scalar defining the dimensions of the correlation matrix.
#' @param min_cor A scalar defining the minimum correlation. By default, min_cor = -1.
#' @param max_cor A scalar defining the maximum correlation. By default, max_cor = 1.
#'
#' @return A p x p correlation matrix.
#'
#' @examples
#' cor_A <- rand_cor_mat(10, min_cor = -0.2, max_cor = 0.8)
#' @export
rand_cor_mat <- function(p,
                         min_cor = -1,
                         max_cor = 1) {
  if (p < 1 | p %% 1 != 0) stop("'p' must be an integer > 0")

  if (min_cor < -1 | min_cor >= 1) stop("'min_cor' must be value >= -1 and < 1")
  if (max_cor <= -1 | min_cor > 1) stop("'max_cor' must be value > -1 and <= 1")
  if (max_cor < min_cor) stop("'max_cor' must not be smaller than 'min_cor'")

  n_cor <- sum(seq(1, (p - 1)))

  off_dg <- stats::runif(n_cor, min = min_cor, max = max_cor)
  cor_mat <- diag(p)
  cor_mat[lower.tri(cor_mat, diag = FALSE)] <- off_dg
  cor_mat <- t(cor_mat)
  cor_mat[lower.tri(cor_mat, diag = FALSE)] <- off_dg

  if (matrixcalc::is.positive.definite(cor_mat)) {
    return(cor_mat)
  } else {
    cor_mat <- mbend::bend(cor_mat)
    return(cor_mat$bent)
  }
}
