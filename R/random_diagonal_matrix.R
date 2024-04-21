#' Simulate a random diagonal variance matrix
#'
#' Creates a diagonal \code{n x n} variance matrix with user-defined minimum and maximum
#' variances based on a continuous uniform distribution.
#'
#' @param n A scalar defining the dimensions of the variance matrix.
#' @param min.var A scalar defining the minimum variance.
#' @param max.var A scalar defining the maximum variance. \cr
#'  \strong{Note:} \code{0 < min.var < max.var}.
#'
#' @return A diagonal \code{n x n} variance matrix. 
#'
#' @examples
#' # Simulate a random diagonal matrix with 10 columns and rows.
#' diag_mat <- rand_diag_mat(
#'   n = 10,
#'   min.var = 0,
#'   max.var = 0.2
#' )
#'
#' @export
rand_diag_mat <- function(n = 5,
                          min.var = 0,
                          max.var = 1) {
  if (!(is.atomic(n) && length(n) == 1L)) stop("'n' must be a scalar")
  if (n < 2 || n %% 1 != 0) stop("'n' must be an integer > 1")
  if (!(is.atomic(min.var) && length(min.var) == 1L)) stop("'min.var' must be a scalar")
  if (min.var < 0) stop("'min.var' must be a value >= 0")
  if (!(is.atomic(max.var) && length(max.var) == 1L)) stop("'max.var' must be a scalar")
  if (max.var < min.var) stop("'max.var' must be a value >= 'min.var'")
  
  dg <- stats::runif(n, min = min.var, max = max.var)
  diag_mat <- diag(dg)
  colnames(diag_mat) <- rownames(diag_mat) <- 1:n
  return(diag_mat)
}

#' Simulate a skewed diagonal variance matrix
#'
#' Creates a diagonal \code{n x n} variance matrix with user-defined skewness
#' based on a gamma or inverse gamma distribution.
#'
#' @param n A scalar defining the dimensions of the variance matrix.
#' @param shape A scalar defining the shape of the distribution.
#' @param scale A scalar defining the scale of the distribution.
#' @param inverse When \code{TRUE} (default is \code{FALSE}), the variances are 
#'   sampled from the inverse gamma distribution instead of the gamma distribution.
#' @param mean.var An optional scalar defining the mean variance. 
#'.  When supplied, the variances are scaled to achieve the defined mean.
#'
#' @return A diagonal \code{n x n} variance matrix. 
#'
#' @examples
#' # Simulate a random diagonal matrix with 10 columns and rows, and negatively skewed variances
#' # scaled to a mean of 0.1.
#' diag_mat <- skew_diag_mat(
#'   n = 10,
#'   shape = 1.5,
#'   scale = 1,
#'   mean.var = 0.1
#' )
#'
#' @export
skew_diag_mat <- function(n = 5,
                          shape = 1.5,
                          scale = 1,
                          inverse = FALSE,
                          mean.var = NULL) {
  if (!(is.atomic(n) && length(n) == 1L)) stop("'n' must be a scalar")
  if (n < 2 || n %% 1 != 0) stop("'n' must be an integer > 1")
  if (!(is.atomic(shape) && length(shape) == 1L)) stop("'shape' must be a scalar")
  if (shape <= 0) stop("'shape' must be a positive value")
  if (!(is.atomic(scale) && length(scale) == 1L)) stop("'scale' must be a scalar")
  if (scale <= 0) stop("'scale' must be a positive value")
  if (!is.logical(inverse)) stop("'inverse' must be logical")
  
  if (!inverse) {
    dg <- stats::rgamma(n, shape = shape, scale = scale)
  } else if (inverse) {
    dg <- 1/stats::rgamma(n, shape = shape, rate = scale)
  }
  
  if (!is.null(mean.var)) {
    if (!(is.atomic(mean.var) && length(mean.var) == 1L)) stop("'mean.var' must be a scalar")
    if (mean.var <= 0) stop("'mean.var' must be a positive value")
    if (mean(dg) == 0) {
      warning("Mean diagonal element is zero, 'mean.var' cannot be applied")
      } else dg <- dg * mean.var/mean(dg)
  }
  
  diag_mat <- diag(dg)
  colnames(diag_mat) <- rownames(diag_mat) <- 1:n
  return(diag_mat)
}