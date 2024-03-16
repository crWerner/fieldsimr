#' Simulate a random correlation matrix
#'
#' Creates a symmetric \code{n x n} correlation matrix with user-defined minimum and maximum
#' correlations.
#'
#' @param n A scalar defining the dimensions of the correlation matrix.
#' @param min.cor A scalar defining the minimum correlation.
#' @param max.cor A scalar defining the maximum correlation. \cr
#'  \strong{Note:} \code{-1 \le min.cor \le max.cor \le 1}.
#' @param pos.def When \code{TRUE} (default is \code{FALSE}), the function \code{bend} of the R package
#'   \href{https://cran.r-project.org/package=mbend}{`mbend`} is used with default arguments to
#'   bend a non-positive-definite matrix to a positive (semi)-definite matrix.
#' @param small.positive Argument passed to \code{bend} when \code{pos.def = TRUE} (default is 0.0001).
#'   Eigenvalues smaller than \code{small.positive} are replaced with this value. \cr
#'  \strong{Note:} \code{0 < small.positive \le 0.1}.
#'
#' @return A symmetric \code{n x n} correlation matrix. When \code{pos.def = TRUE},
#' the correlation matrix is guaranteed to be positive (semi)-definite.
#'
#' @examples
#' corA <- rand_cor_mat(
#'   n = 10,
#'   min.cor = -0.2,
#'   max.cor = 0.8,
#'   pos.def = TRUE,
#'   small.positive = 1e-4
#' )
#'
#' @export
rand_cor_mat <- function(n = 2,
                         min.cor = -1,
                         max.cor = 1,
                         pos.def = FALSE,
                         small.positive = NULL) {
  if (n < 1 || n %% 1 != 0) stop("'n' must be a positive integer")

  if (min.cor < -1 || min.cor > 1) stop("'min.cor' must be a value >= -1 and <= 1")
  if (max.cor < -1 || max.cor > 1) stop("'max.cor' must be a value >= -1 and <= 1")
  if (max.cor < min.cor) stop("'max.cor' must be >= 'min.cor'")

  ncor <- sum(seq(1, (n - 1)))

  off_dg <- stats::runif(ncor, min = min.cor, max = max.cor)
  cor_mat <- diag(n)
  cor_mat[lower.tri(cor_mat, diag = FALSE)] <- off_dg
  cor_mat <- t(cor_mat)
  cor_mat[lower.tri(cor_mat, diag = FALSE)] <- off_dg
  colnames(cor_mat) <- rownames(cor_mat) <- 1:n

  is_pos_def <- sum(eigen(cor_mat)$values > 1e-8) == n
  if (pos.def && !is_pos_def) {
    if (is.null(small.positive)) {
      small.positive <- 1e-4
    }
    if (small.positive <= 0 | small.positive > 0.1) stop("'small.positive' must be a positive value <= 0.1")
    cor_mat <- mbend::bend(cor_mat, small.positive = small.positive)
    cor_mat <- round(cor_mat$bent, 12)
    cor_mat <- (cor_mat + t(cor_mat)) / 2
  }
  return(cor_mat)
}
