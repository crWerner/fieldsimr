# ’ Linear vector interpolation
# ’
# ’ Linear interpolation of a vector using the "stats::approx" function.
# ’ @noRd
interp_vector <- function(y) {
  finite_vals <- is.finite(y)
  if (sum(finite_vals) == 0) {
    return(rep(NA, length(y)))
  }
  if (sum(finite_vals) == 1) {
    return(rep(y[finite_vals], length(y)))
  }
  stats::approx(seq_along(y), y, seq_along(y))$y
}


# ’ Linear vector extrapolation: head of a vector
# ’
# ’ Extrapolation of the head of a vector using the "stats::predict" function.
# ’ @noRd
extrap_head <- function(y) {
  x <- seq_along(y)
  finite_head <- which(!is.na(y))[1:2]
  head_NA <- which(is.na(y))
  head_NA <- head_NA[head_NA < finite_head[1]]

  y[head_NA] <- stats::predict(stats::lm(y ~ x, data = data.frame(x, y)[finite_head, ]), list(x = x[head_NA]))
  round(y, 6) # round is necessary. Check length!
}


# ’ Linear vector extrapolation: tail of a vector
# ’
# ’ Extrapolation of the tail of a vector using the "stats::predict" function.
# ’ @noRd
extrap_tail <- function(y) {
  x <- seq_along(y)
  finite_tail <- utils::tail(which(!is.na(y)), 2)
  tail_NA <- which(is.na(y))
  tail_NA <- tail_NA[tail_NA > finite_tail[2]]

  y[tail_NA] <- stats::predict(stats::lm(y ~ x, data = data.frame(x, y)[finite_tail, ]), list(x = x[tail_NA]))
  round(y, 6) # round is necessary. Check length!
}


# ’ Linear interpolation and extrapolation of a vector
# ’
# ’ Linear interpolation and extrapolation of a vector using the functions "interp_vector",
# ’ "extrap_head", and "extrap_tail".
# ’ @noRd
fill_vector <- function(y) {
  y <- interp_vector(y)
  x <- seq_along(y)
  if (all(is.na(y))) {
    return(y)
  }
  if (is.na(y[1])) y <- extrap_head(y)
  if (is.na(utils::tail(y, 1))) y <- extrap_tail(y)
  y
}


# ’ Bilinear matrix interpolation and extrapolation
# ’
# ’ Bilinear interpolation and extrapolation of a matrix. First, missing values of each row and
# ’ column of the matrix are interpolated or extrapolated individually. Then, missing values are
# ’ replaced using the mean of the row-wise and column-wise interpolation or extrapolation.
# ’ @noRd
fill_matrix <- function(mat) {
  rowwise <- t(apply(mat, 1, fill_vector))
  colwise <- apply(mat, 2, fill_vector)
  colwise[is.na(colwise)] <- rowwise[is.na(colwise)]
  rowwise[is.na(rowwise)] <- colwise[is.na(rowwise)]
  mat[is.na(mat)] <- ((colwise + rowwise) / 2)[is.na(mat)]
  mat
}
