#' Graphics for plot effects
#'
#' Graphically displays plot effects (e.g., phenotypic values, genetic values, error terms) onto a
#' field array, in which the colour gradient ranges from red (low value) to green (high value).
#' The function requires a data frame generated with field_trial_error as an input, or any data
#' frame with columns “env”, “col”, “row”, and the effect to be displayed. If the data frame
#' contains a “block” column, the field array is split into blocks if \code{blocks = TRUE}.
#'
#' @param df A data frame with columns "env", "col", "row", and the effect to be plotted.
#'   If \code{df} contains a “block” column, the field array is split into blocks if
#'   \code{blocks = TRUE}. \cr
#'   \strong{Note:} If \code{df} is a list, only the first entry will be used unless specified
#'   otherwise.
#' @param env The ID of the environment to be plotted.
#' @param effect The effect to be plotted.
#' @param blocks When TRUE (default), the field array is split into blocks.
#'
#' @return A graphical field array, in which the colour gradient ranges from red (low value)
#' to green (high value).
#'
#' @examples
#' # Plot the simulated total error term for trait 2 in environment 2 provided in the example data
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
    stop("'df' must contain columns 'env', 'col', 'row', and the effect to be plotted.")
  }

  df <- df[df[["ENV"]] == env, ]
  colnames(df)[colnames(df) %in% effect] <- "EFF"

  n_rows <- length(unique(df$ROW))
  n_cols <- length(unique(df$COL))
  n_blocks <- length(unique(df$BLOCK))

  if (n_blocks > 1) {
    df1 <- df[df[["BLOCK"]] == 1, ]
    df2 <- df[df[["BLOCK"]] == 2, ]

    if (any(unique(df[df[["BLOCK"]] == 1, ]$ROW) == unique(df[df[["BLOCK"]] == 2, ]$ROW)) == FALSE) {
      dist <- (n_rows / n_blocks)
      x_min <- rep(0.5, n_blocks)
      y_min <- (seq(0, n_rows, dist) + 0.5)[1:n_blocks]
      x_max <- rep((n_cols + 0.5), n_blocks)
      y_max <- (seq(0, n_rows, dist) + 0.5)[2:(n_blocks + 1)]

      block_x_min <- rep(0.5, n_blocks)
      block_y_min <- y_max
      block_x_max <- rep((n_cols + 0.5), n_blocks)
      block_y_max <- y_max

      block_x_min_2 <- rep(0, n_blocks)
      block_y_min_2 <- y_max
      block_x_max_2 <- rep((n_cols + 1), n_blocks)
      block_y_max_2 <- y_max
    } else if (any(unique(df[df[["BLOCK"]] == 1, ]$COL) == unique(df[df[["BLOCK"]] == 2, ]$COL)) == FALSE) {
      dist <- (n_cols / n_blocks)
      x_min <- (seq(0, n_cols, dist) + 0.5)[1:n_blocks]
      y_min <- rep(0.5, n_blocks)
      x_max <- (seq(0, n_cols, dist) + 0.5)[2:(n_blocks + 1)]
      y_max <- rep((n_rows + 0.5), n_blocks)

      block_x_min <- x_max
      block_y_min <- rep(0.5, n_blocks)
      block_x_max <- x_max
      block_y_max <- rep((n_rows + 0.5), n_blocks)

      block_x_min_2 <- x_max
      block_y_min_2 <- rep(0, n_blocks)
      block_x_max_2 <- x_max
      block_y_max_2 <- rep((n_rows + 1), n_blocks)
    } else {
      stop("Check row and column assignment within blocks")
    }
  }

  COL <- ROW <- EFF <- NULL

  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = COL, y = ROW)) +
    ggplot2::geom_tile(ggplot2::aes(fill = EFF)) +
    ggplot2::scale_fill_gradient2(low = "#A51122", mid = "#FEFDBE", high = "#006228") +
    ggplot2::xlab("Columns") +
    ggplot2::ylab("Rows") +
    ggplot2::theme_grey(base_size = 10) +
    ggplot2::ggtitle(effect) +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 12, colour = "gray40")
    ) +
    ggplot2::annotate(
      geom = "rect", xmin = 0.5, ymin = 0.5,
      xmax = n_cols + 0.5, ymax = n_rows + 0.5,
      fill = "transparent", col = "black", lwd = 0.5
    )

  if (blocks == TRUE & length(unique(df$BLOCK)) > 1) {
    for (i in 1:(n_blocks - 1)) {
      p <- p + ggplot2::geom_segment(
        x = block_x_min[i],
        y = block_y_min[i],
        xend = block_x_max[i],
        yend = block_y_max[i],
        size = 1.5
      ) +
        ggplot2::geom_segment(
          x = block_x_min_2[i],
          y = block_y_min_2[i],
          xend = block_x_max_2[i],
          yend = block_y_max_2[i],
          size = 1,
          col = "white"
        )
    }
  }
  return(p)
}
