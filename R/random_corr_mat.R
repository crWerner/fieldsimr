#' Random correlation matrix
#'
#' Creates a symmetric \code{p x p} correlation matrix with user-defined minimum and maximum
#' correlation values.
#'
#' @param p A scalar defining the dimensions of the correlation matrix.
#' @param min_cor A scalar defining the minimum potential value.
#' @param max_cor A scalar defining the maximum potential value.
#' @param pos_def When \code{TRUE}, the function 'bend' of the package
#'   \href{https://cran.r-project.org/package=mbend}{`mbend'} is used with default arguments to
#'   bend a non-positive-definite correlation matrix to a positive-definite matrix (when
#'   appropriate). By default, \code{pos_def = FALSE}.
#'
#' @return A symmetric \code{p x p} correlation matrix.
#'
#' @examples
#' cor_A <- rand_cor_mat(10, min_cor = -0.2, max_cor = 0.8, pos_def = TRUE)
#' @export
rand_cor_mat <- function(p,
                         min_cor = -1,
                         max_cor = 1,
                         pos_def = FALSE) {
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

  if (pos_def && !matrixcalc::is.positive.definite(cor_mat)) {
    cor_mat <- mbend::bend(cor_mat)
    cor_mat <- round(cor_mat$bent, 12)
    cor_mat <- (cor_mat + t(cor_mat)) / 2
  }
  return(cor_mat)
}
