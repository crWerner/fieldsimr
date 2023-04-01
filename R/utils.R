#' Spline interpolation
#'
#' Interpolate missing values in a vector using natural splines. Linear extrapolation of
#' missing values is also executed.
#' @noRd
spline_interp <- function(y) {
  finite_vals <- is.finite(y)
  if (sum(finite_vals) == 0) {
    return(rep(NA, length(y)))
  }
  if (sum(finite_vals) == 1) {
    return(rep(y[finite_vals], length(y)))
  }
  stats::spline(x = seq_along(y), y = y, xout = seq_along(y), method = "natural")$y
}


#' Matrix interpolation
#'
#' Interpolate and extrapolate missing values in a matrix. First, missing values in each
#' column and row of the matrix are interpolated and extrapolated individually. Then, missing
#' values in each matrix element are replaced using the mean of the intersecting column-wise and row-wise interpolation and/or
#' extrapolation values.
#' @noRd
fill_matrix <- function(mat) {
  rowwise <- t(apply(mat, 1, spline_interp))
  colwise <- apply(mat, 2, spline_interp)
  colwise[is.na(colwise)] <- rowwise[is.na(colwise)]
  rowwise[is.na(rowwise)] <- colwise[is.na(rowwise)]
  mat[is.na(mat)] <- ((colwise + rowwise) / 2)[is.na(mat)]
  mat
}
