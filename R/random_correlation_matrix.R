#' Simulate a random correlation matrix
#'
#' Creates a symmetric \code{n x n} correlation matrix with user-defined minimum and maximum
#' correlations based on a continuous uniform distribution.
#'
#' @param n A scalar defining the dimensions of the correlation matrix.
#' @param min.cor A scalar defining the minimum correlation.
#' @param max.cor A scalar defining the maximum correlation. \cr
#'  \strong{Note:} \code{-1 < min.cor < max.cor < 1}.
#' @param pos.def When \code{TRUE} (default is \code{FALSE}), the function \code{bend} of the R package
#'   \href{https://cran.r-project.org/package=mbend}{`mbend`} is used to
#'   bend a non-positive (semi)-definite matrix to be positive (semi)-definite.
#' @param small.positive Argument passed to \code{bend} when \code{pos.def = TRUE} (default is 1e-8).
#'   Eigenvalues smaller than \code{small.positive} are replaced by this. \cr
#'  \strong{Note:} \code{0 < small.positive < 0.1}.
#'
#' @return A symmetric \code{n x n} correlation matrix. When \code{pos.def = TRUE},
#' the correlation matrix is guaranteed to be positive (semi)-definite.
#'
#' @examples
#' # Simulate and visualise a random correlation matrix with 10 columns and rows.
#' cor_mat <- rand_cor_mat(
#'   n = 10,
#'   min.cor = -0.2,
#'   max.cor = 0.8,
#'   pos.def = TRUE
#' )
#'
#' plot_matrix(
#'   mat = cor_mat,
#'   order = TRUE
#' )
#'
#' @export
rand_cor_mat <- function(n = 5,
                         min.cor = -1,
                         max.cor = 1,
                         pos.def = FALSE,
                         small.positive = NULL) {
  if (!(is.atomic(n) && length(n) == 1L)) stop("'n' must be a scalar")
  if (n < 2 || n %% 1 != 0) stop("'n' must be an integer > 1")
  if (!(is.atomic(min.cor) && length(min.cor) == 1L)) stop("'min.cor' must be a scalar")
  if (min.cor < -1 || min.cor > 1) stop("'min.cor' must be a value >= -1 and <= 1")
  if (!(is.atomic(max.cor) && length(max.cor) == 1L)) stop("'max.cor' must be a scalar")
  if (max.cor < -1 || max.cor > 1) stop("'max.cor' must be a value >= -1 and <= 1")
  if (max.cor < min.cor) stop("'max.cor' must be a value >= 'min.cor'")
  if (!is.logical(pos.def)) stop("'pos.def' must be logical")

  ncor <- sum(seq(1, (n - 1)))

  off_dg <- stats::runif(ncor, min = min.cor, max = max.cor)
  cor_mat <- diag(n)
  cor_mat[lower.tri(cor_mat, diag = FALSE)] <- off_dg
  cor_mat <- t(cor_mat)
  cor_mat[lower.tri(cor_mat, diag = FALSE)] <- off_dg
  colnames(cor_mat) <- rownames(cor_mat) <- 1:n

  if (pos.def) {
    is_pos_semi_def <- all(eigen(cor_mat, only.values = TRUE)$values >= 0)
    if (is_pos_semi_def) {
      message("'cor_mat' is already positive (semi)-definite, matrix was not altered")
    } else if (!is_pos_semi_def) {
      if (is.null(small.positive)) {
        small.positive <- 1e-8
      }
      if (small.positive <= 0 | small.positive > 0.1) stop("'small.positive' must be a positive value <= 0.1")
      cor_mat <- mbend::bend(cor_mat, small.positive = small.positive)$bent
    }
  }
  if (any(abs(cor_mat) > 1)) stop("absolute values greater than 1 produced")
  return(cor_mat)
}

#' Simulate a structured correlation matrix with reduced rank
#'
#' Creates a symmetric \code{n x n} correlation matrix with user-defined structure and rank.
#'
#' @param n A scalar defining the dimensions of the correlation matrix.
#' @param base.cor A scalar defining the baseline correlation. \cr
#'   \strong{Note:} \code{-1 < base.cor < 1}.
#' @param range A scalar defining the range of correlations around the baseline. By default,
#'   \code{range = 1 - base.cor} which ensures the matrix is positive semi-definite with defined rank.
#'   \strong{Note:} \code{base.cor + range <= 1}.
#' @param rank A scalar defining the rank of the correlation matrix.
#' @param skew A scalar defining the skewness imposed on the correlations.
#'   \strong{Note:} \code{-1 < skew < 1}.
#' @param base.mat An optional \code{n x n} base correlation matrix.
#'   When supplied, \code{base.cor} and \code{skew} are ignored and noise is simulated based on \code{rank}.
#' @param pos.def When \code{TRUE} (default is \code{FALSE}), the function \code{bend} of the R package
#'   \href{https://cran.r-project.org/package=mbend}{`mbend`} is used to
#'   bend a non-positive (semi)-definite matrix to be positive (semi)-definite.
#' @param small.positive Argument passed to \code{bend} when \code{pos.def = TRUE} (default is 1e-8).
#'   Eigenvalues smaller than \code{small.positive} are replaced by this. \cr
#'   \strong{Note:} \code{0 < small.positive < 0.1}.
#'
#' @return A symmetric \code{n x n} correlation matrix with defined rank.
#'   When \code{pos.def = TRUE}, the correlation matrix is guaranteed to be positive (semi)-definite.
#'
#' @examples
#' # Simulate and visualise a correlation matrix with 10 columns and rows, rank equal to 4 and
#' # negatively skewed correlations.
#' cor_mat <- struc_cor_mat(
#'   n = 10,
#'   base.cor = 0.3,
#'   range = 0.7,
#'   rank = 4,
#'   skew = -0.5
#' )
#'
#' plot_matrix(
#'   mat = cor_mat,
#'   order = TRUE
#' )
#'
#' @export
struc_cor_mat <- function(n = 5,
                          base.cor = 0.5,
                          range = NULL,
                          rank = 3,
                          skew = 0,
                          base.mat = NULL,
                          pos.def = FALSE,
                          small.positive = NULL) {
  if (!(is.atomic(n) && length(n) == 1L)) stop("'n' must be a scalar")

  if (is.null(base.mat)) {
    if (n < 2 || n %% 1 != 0) stop("'n' must be an integer > 1")
    if (!is.logical(pos.def)) stop("'pos.def' must be logical")
    if (!(is.atomic(base.cor) && length(base.cor) == 1L)) stop("'base.cor' must be a scalar")
    if (base.cor < -1 || base.cor > 1) stop("'base.cor' must be a value >= -1 and <= 1")

    abs_base_cor <- abs(base.cor)
    if (is.null(range)) {
      range <- 1 - abs_base_cor
    }
    if (!(is.atomic(range) && length(range) == 1L)) stop("'range' must be a scalar")
    if (range < 0) stop("'range' must be a positive number")
    if (abs_base_cor + range > 1) stop("The absolute sum of 'base.cor' and 'range' must be <= 1")
    if (!(is.atomic(rank) && length(rank) == 1L)) stop("'rank' must be a scalar")
    if (rank < 1 || rank %% 1 != 0) stop("'rank' must be a positive integer")
    if (!(is.atomic(skew) && length(skew) == 1L)) stop("'skew' must be a scalar")
    if (skew < -1 || skew > 1) stop("'skew' must be a value >= -1 and <= 1")

    if (abs_base_cor + range < 1) warning("The absolute sum of 'base.cor' and 'range' is not 1, defined rank will not be achieved")
    if (!pos.def) {
      insertion <- paste0("matrix will be indefinite and ")
    } else {
      insertion <- NULL
    }
    if (base.cor < 0) warning("'base.cor' is less than 0, ", insertion, "defined rank will not be achieved")
    if (skew > 0) warning("'skew' is greater than 0, ", insertion, "defined rank will not be achieved")
    if (rank > n) warning("'rank' greater than 'n', defined rank is not possible")
    if (rank == 1 && base.cor != 1) stop("'rank' must be greater than 1 when 'base.cor' is not 1")
    if (rank != 1 && base.cor == 1) {
      warning("'rank' set to 1 since 'base.cor' is 1")
      rank <- 1
    }
    if (!(n == 2 && rank == 2) && base.cor != 0 || rank == 1) {
      rank <- rank - 1
    }
    base_mat <- matrix(rep(base.cor, n * n), ncol = n)
  } else if (!is.null(base.mat)) {
    if (is.null(range)) {
      range <- 0
    }
    if (!(is.atomic(range) && length(range) == 1L)) stop("'range' must be a scalar")
    if (range < 0) stop("'range' must be a positive number")
    if (!is.matrix(base.mat)) stop("'base.mat' must be a matrix")
    if (!isSymmetric(base.mat)) stop("'base.mat' must be a symmetric matrix")
    if (any(is.na(base.mat))) stop("'base.mat' must not contain missing values")
    max_abs_base_cor <- max(abs(base.mat[upper.tri(base.mat)]))
    if (max_abs_base_cor + range > 1) stop("The absolute sum of 'range' and the maximum correlation in 'base.mat' must be <= 1")
    warning("'base.mat' supplied, 'n' and 'base.cor' will be ignored")
    base_mat <- base.mat
    n <- ncol(base_mat)
  }

  skew_abs <- abs(skew)
  diag(base_mat) <- 1 - range
  if (rank > 0) {
    eigen_ls <- lapply(seq_len(n), function(x) cbind(stats::runif(rank, -1, 1 - skew_abs)))
    eigen_ls <- lapply(eigen_ls, function(x) x / sqrt(sum(x^2)))
    eigen_vect <- do.call(cbind, eigen_ls)
    cor_mat <- base_mat + range * t(eigen_vect) %*% eigen_vect

    if (skew > 0) {
      cor_mat <- -cor_mat
    }
  } else {
    cor_mat <- base_mat
  }
  diag(cor_mat) <- 1
  colnames(cor_mat) <- rownames(cor_mat) <- 1:n

  if (pos.def) {
    is_pos_semi_def <- all(eigen(cor_mat, symmetric = TRUE)$values >= 0)
    if (is_pos_semi_def) {
      message("'cor_mat' is already positive (semi)-definite, matrix was not altered")
    } else if (!is_pos_semi_def) {
      if (is.null(small.positive)) {
        small.positive <- 1e-8
      }
      if (!(is.atomic(small.positive) && length(small.positive) == 1L)) stop("'small.positive' must be a scalar")
      if (small.positive <= 0 | small.positive > 0.1) stop("'small.positive' must be a positive value <= 0.1")
      cor_mat <- mbend::bend(cor_mat, small.positive = small.positive)$bent
    }
  }
  if (any(abs(cor_mat) > 1)) stop("absolute values greater than 1 produced")
  return(cor_mat)
}

#' Simulate a reduced rank correlation matrix with multiple groups
#'
#' Creates a symmetric correlation matrix with user-defined structure, rank and groupings.
#'
#' @param n A vector defining the size of each group.
#' @param within.cor A vector defining the baseline correlation within each group.
#'   If only one value is supplied, all groups will be assigned the same correlation. \cr
#'   \strong{Note:} \code{-1 < within.cor < 1}.
#' @param between.cor A scalar defining the baseline correlation between groups. \cr
#'   \strong{Note:} \code{between.cor <= within.cor}.
#' @param range A scalar defining the range of correlations around the baseline. By default,
#'   \code{range = 1 - max(within.cor)} which ensures the matrix is positive semi-definite.
#'   \strong{Note:} \code{max(within.cor) + range <= 1}.
#' @param rank A scalar defining the rank of the correlation matrix.
#' @param skew A scalar defining the skewness imposed on the correlations.
#'   \strong{Note:} \code{-1 < skew < 1}.
#' @param pos.def When \code{TRUE} (default is \code{FALSE}), the function \code{bend} of the R package
#'   \href{https://cran.r-project.org/package=mbend}{`mbend`} is used to
#'   bend a non-positive (semi)-definite matrix to be positive (semi)-definite.
#' @param small.positive Argument passed to \code{bend} when \code{pos.def = TRUE} (default is 1e-8).
#'   Eigenvalues smaller than \code{small.positive} are replaced by this. \cr
#'   \strong{Note:} \code{0 < small.positive < 0.1}.
#' @param return.groups When \code{TRUE} (default is \code{FALSE}), a list is returned with additional
#'   entries containing the members of each group.
#'
#' @return A symmetric correlation matrix with defined rank and groupings.
#'   When \code{pos.def = TRUE}, the correlation matrix is guaranteed to be positive (semi)-definite.
#'   When \code{return.groups = TRUE}, a list is returned with additional entries containing the group members.
#'
#' @examples
#' # Simulate and visualise a correlation matrix with 2 groups containing 5 and 10 members,
#' # correlations of 0.4 within groups and 0 between groups and rank equal to 4
#' cor_ls <- group_cor_mat(
#'   n = c(5, 10),
#'   within.cor = 0.4,
#'   between.cor = 0,
#'   rank = 4,
#'   return.groups = TRUE
#' )
#'
#' plot_matrix(
#'   mat = cor_ls$cor.mat,
#'   group.df = cor_ls$group.df,
#'   order = TRUE
#' )
#'
#' @export
group_cor_mat <- function(n = c(5, 5),
                          within.cor = 0.5,
                          between.cor = 0.2,
                          range = NULL,
                          rank = 4,
                          skew = 0,
                          pos.def = FALSE,
                          small.positive = NULL,
                          return.groups = FALSE) {
  if (any(n %% 1 != 0)) stop("'n' must contain integers")
  ntotal <- sum(n)
  ngroups <- length(n)
  if (!is.logical(pos.def)) stop("'pos.def' must be logical")
  if (length(within.cor) == 1) {
    within.cor <- rep(within.cor, ngroups)
  }
  if (length(within.cor) != ngroups) stop("Length of 'within.cor' must match length of 'n'")
  if (any(within.cor < -1) || any(within.cor > 1)) stop("'within.cor' must be a value >= -1 and <= 1")
  if (!(is.atomic(between.cor) && length(between.cor) == 1L)) stop("'between.cor' must be a scalar")
  if (any(between.cor < -1) || any(between.cor > 1)) stop("'between.cor' must be a value >= -1 and <= 1")

  abs_within_cor <- abs(within.cor)
  abs_between_cor <- abs(between.cor)
  if (is.null(range)) {
    range <- 1 - max(abs_within_cor)
  }

  if (!(is.atomic(range) && length(range) == 1L)) stop("'range' must be a scalar")
  if (range < 0) stop("'range' must be a positive number")
  if (any(between.cor > within.cor)) {
    if (!pos.def) warning("between.cor' is greater than some values in 'within.cor', matrix will be indefinite")
    if (any(abs_between_cor + range > 1)) stop("The absolute sum of 'between.cor' and 'range' must be <= 1")
  } else if (any(abs_within_cor + range > 1)) stop("The absolute sum of each value in 'within.cor' and 'range' must be <= 1")
  if (!(is.atomic(rank) && length(rank) == 1L)) stop("'rank' must be a scalar")
  if (rank < 1 || rank %% 1 != 0) stop("'rank' must be a positive integer")
  if (!(is.atomic(skew) && length(skew) == 1L)) stop("'skew' must be a scalar")
  if (skew < -1 || skew > 1) stop("'skew' must be a value >= -1 and <= 1")

  if (length(unique(abs_within_cor)) > 1 && rank < ntotal) {
    warning("More than one unique value supplied in 'within.cor', defined rank will not be achieved")
  } else if (any(abs_within_cor + range < 1)) warning("The absolute sum of 'within.cor' and 'range' is not 1, defined rank will not be achieved")
  if (!pos.def) {
    insertion <- paste0("matrix will be indefinite and ")
  } else {
    insertion <- NULL
  }
  if (any(within.cor < 0)) warning("some values in 'within.cor' are less than 0, ", insertion, "defined rank will not be achieved")
  if (between.cor < 0) warning("'between.cor' is less than 0, ", insertion, "defined rank will not be achieved")
  if (skew > 0) warning("'skew' is greater than 0, ", insertion, "defined rank will not be achieved")
  if (rank > ntotal) warning("'rank' greater than sum of 'n', defined rank is not possible")

  if (all(within.cor == 1)) {
    if (rank != ngroups && between.cor != 1) {
      warning("'rank' set to number of groups since all values in 'within.cor' are 1")
    } else if (rank != 1 && between.cor == 1) {
      warning("'rank' set to 1 since all values in 'within.cor' and 'between.cor' are 1")
    }
    rank <- ngroups
  }
  if (rank <= ngroups && any(within.cor != 1) && any(within.cor != between.cor)) stop("'rank' must be greater than number of groups when all values in 'within.cor' do not equal 'between.cor'")
  if (rank == 1 && any(c(within.cor, between.cor) != 1)) stop("'rank' must be greater than 1 when all values in 'within.cor' and 'between.cor' are not 1")

  rank <- rank - ngroups
  if (all(between.cor == within.cor) && any(!within.cor %in% c(0, 1))) {
    rank <- rank + (ngroups - 1)
  }

  base_mat <- mapply(function(x, y) matrix(rep(x, y * y), ncol = y) - between.cor, x = within.cor, y = n, SIMPLIFY = FALSE)
  base_mat <- as.matrix(Reduce(Matrix::bdiag, base_mat)) + between.cor

  skew_abs <- abs(skew)
  diag(base_mat) <- 1 - range
  if (rank > 0) {
    eigen_ls <- lapply(seq_len(ntotal), function(x) cbind(stats::runif(rank, -1, 1 - skew_abs)))
    eigen_ls <- lapply(eigen_ls, function(x) x / sqrt(sum(x^2)))
    eigen_vect <- do.call(cbind, eigen_ls)
    cor_mat <- base_mat + range * t(eigen_vect) %*% eigen_vect

    if (skew > 0) {
      cor_mat <- -cor_mat
    }
  } else {
    cor_mat <- base_mat
  }
  diag(cor_mat) <- 1
  colnames(cor_mat) <- rownames(cor_mat) <- 1:ntotal

  if (pos.def) {
    is_pos_semi_def <- all(eigen(cor_mat, symmetric = TRUE)$values >= 0)
    if (is_pos_semi_def) {
      message("'cor_mat' is already positive (semi)-definite, matrix was not altered")
    } else if (!is_pos_semi_def) {
      if (is.null(small.positive)) {
        small.positive <- 1e-8
      }
      if (!(is.atomic(small.positive) && length(small.positive) == 1L)) stop("'small.positive' must be a scalar")
      if (small.positive <= 0 | small.positive > 0.1) stop("'small.positive' must be a positive value <= 0.1")
      cor_mat <- mbend::bend(cor_mat, small.positive = small.positive)$bent
    }
  }
  if (any(abs(cor_mat) > 1)) stop("absolute values greater than 1 produced")

  if (return.groups) {
    group_df <- data.frame(variable = 1:ntotal, group = rep(1:ngroups, n))
    cor_mat <- list(cor.mat = cor_mat, group.df = group_df)
  }
  return(cor_mat)
}
