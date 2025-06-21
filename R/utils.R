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


#' Check object class
#'
#' Checks the class of objects, with particular focus on numeric vs character factors.
#' @noRd
check_class <- function(object) {
  if (is.numeric(object)) {
    return("numeric")
  } else if (is.character(object)) {
    suppressWarnings(character_is_numeric <- all(!is.na(as.numeric(unique(object)))))
    if (character_is_numeric) {
      return("numeric character")
    } else {
      return("character")
    }
  } else if (is.factor(object)) {
    suppressWarnings(levels_are_numeric <- all(!is.na(as.numeric(levels(object)))))
    if (levels_are_numeric) {
      return("numeric factor")
    } else {
      return("character factor")
    }
  } else {
    return(class(object))
  }
}


#' Make factor
#'
#' Converts an object to a factor, preserving the order for character objects.
#' @noRd
make_factor <- function(object) {
  if (check_class(object) %in% c("numeric", "numeric character")) {
    return(factor(as.numeric(object)))
  } else if (check_class(object) %in% "character") {
    return(factor(object, levels = unique(object)))
  } else {
    return(object)
  }
}


#' Pad out a data frame
#'
#' Pads out a data frame to include all combinations of the two columns specified,
#' creating dummy rows for missing values.
#' @noRd
padout_df <- function(df, var1 = "env", var2 = "id") {
  df[[var1]] <- make_factor(df[[var1]])
  df[[var2]] <- make_factor(df[[var2]])
  all_combos <- expand.grid(
    levels(df[[var1]]),
    levels(df[[var2]]),
    stringsAsFactors = FALSE
  )
  names(all_combos) <- c(var1, var2)
  padded_df <- merge(all_combos, df, by = c(var1, var2), all.x = TRUE, sort = FALSE)
  all_cols <- union(names(df), names(padded_df))
  padded_df <- padded_df[, all_cols]
  padded_df[[var1]] <- make_factor(padded_df[[var1]])
  padded_df[[var2]] <- make_factor(padded_df[[var2]])
  return(padded_df)
}
