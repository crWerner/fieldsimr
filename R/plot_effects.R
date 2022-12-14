#' Graphics for plot effects
#'
#' Graphically displays plot effects (e.g., phenotypic values, genetic values, errors) onto
#' a field array, where the colour gradient ranges from red (low value) to green (high value). \cr
#' This function requires a data frame generated with \link[FieldSimR]{field_trial_error} as an
#' input, or any data frame with columns named "env", "col", "row", and the effect to be displayed.
#' If the data frame contains a column named "block", then block borders will distinguish the
#' blocks if \code{blocks = TRUE}.
#'
#' @param df A data frame containing the columns "env", "row", "col", and the effect to be plotted.
#'   If \code{df} contains a column named "block", then block borders will distinguish the blocks
#'   if \code{blocks = TRUE}. If \code{df} is a list, only the first entry will be used unless
#'   otherwise specified.
#' @param env The name of the environment to be plotted.
#' @param effect The name of the effect to be plotted.
#' @param blocks When TRUE (default), blocks are distinguished with block borders.
#'
#' @return Graphic of the field array, where the colour gradient ranges from red (low value) to
#'   green (high value) of the effect
#'
#' @examples
#' # Plot the simulated error for trait 2 in environment 2 of the field trial error example data
#' # frame 'df_error_bivar'.
#'
#' error_df <- df_error_bivar
#'
#' plot_effects(
#'   error_df,
#'   env = 2,
#'   effect = "e.Trait.2"
#' )
#' @export
plot_effects <- function(df,
                         env,
                         effect,
                         blocks = TRUE) {
  if (inherits(df, "list")) df <- data.frame(df[[1]])

  colnames(df) <- toupper(colnames(df))
  effect <- toupper(effect)
  if (any(!c("ENV", "COL", "ROW") %in% colnames(df))) {
    stop("'df' must contain columns 'env', 'row', 'col', and the effect to be plotted.")
  }

  eff <- which(colnames(df) == effect)
  df <- df[df[["ENV"]] == env, ]
  n_rows <- length(unique(df$ROW))
  n_cols <- length(unique(df$COL))

  plot_mat <- matrix(numeric(), nrow = n_rows, ncol = n_cols)

  for (i in 1:nrow(df)) {
    r <- df$ROW[i]
    c <- df$COL[i]
    plot_mat[r, c] <- df[i, eff]
  }

  if (length(unique(df$BLOCK)) > 1) {
    df1 <- df[df[["BLOCK"]] == 1, ]
    df2 <- df[df[["BLOCK"]] == 2, ]

    if (any(unique(df1$ROW) == unique(df2$ROW)) == FALSE) {
      nx <- 0
      ny <- length(unique(df$BLOCK))
    } else if (any(unique(df1$COL) == unique(df2$COL)) == FALSE) {
      nx <- length(unique(df$BLOCK))
      ny <- 0
    } else {
      stop("Check row and column assignment within blocks")
    }
  }

  x_labs <- seq(2, ncol(plot_mat), 2)
  x_ticks <- (seq(2, ncol(plot_mat), 2) - 0.5)

  y_labs <- seq(2, nrow(plot_mat), 2)
  y_ticks <- (seq(2, nrow(plot_mat), 2) - 0.5)

  fields::image.plot(
    x = 0:n_cols[1], y = 0:n_rows[1],
    z = t(plot_mat), zlim = range(plot_mat),
    ylim = rev(range(0:n_rows[1])),
    col = grDevices::hcl.colors(n = 10000, "RdYlGn"),
    xlab = "Column", ylab = "Row", axes = FALSE
  )

  graphics::box()
  graphics::axis(1, at = x_ticks, labels = x_labs)
  graphics::axis(2, at = y_ticks, labels = y_labs)

  if (blocks == TRUE & length(unique(df$BLOCK)) > 1) {
    graphics::grid(nx = nx, ny = ny, lty = 1, col = "#000000", lwd = 5)
    graphics::grid(nx = nx, ny = ny, lty = 1, col = "#FFFFFF", lwd = 3)
  }
}
