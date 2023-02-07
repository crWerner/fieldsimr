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

#' Construct a qqplot
#'
#' @param df A data frame containing the columns "col", "row", and the effect to be plotted.
#' @param effect The name of the effect to be plotted.
#' @param labels When TRUE (default), the column and row labels are inserted onto the qqplot.
#'   Otherwise, data points without labels are plotted. In this case, the data frame does not require
#'   the columns "col" and "row", just the effect to be plotted.
#' @param plot When TRUE (default), the qqplot is displayed graphically.
#'   Otherwise, a data frame is returned.
#'
#' @return Graphic of the qqplot, where the x- and y- axes display the theoretical and
#'   sample quantiles. When \code{plot = FALSE}, a data frame is returned with the column and row displacements
#'   as well as the estimated semi-variances.
#'
#' @examples
#' # Plot the qqplot for the spatial error component in the data frame error_df.
#'
#' qq_plot(
#'   error_df,
#'   effect = "e_spat",
#'   labels = TRUE,
#'   plot = TRUE,
#' )
#' @export
#'
#'

qq_plot <- function(df,
                    effect,
                    labels = TRUE,
                    plot = TRUE) {
  if (!labels) {
    qq.df <- data.frame(effect = df[[effect]])
    p1 <- ggplot2::ggplot(qq.df, ggplot2::aes(sample = effect)) +
      ggplot2::stat_qq()
    qq.df <- data.frame(
      sample = ggplot2::ggplot_build(p1)$data[[1]]["sample"],
      theo = ggplot2::ggplot_build(p1)$data[[1]]["theoretical"]
    )
    if (!plot) {
      colnames(qq.df)[2] <- "theoretical"
      print(qq.df)
    }
    if (plot) {
      p1 <- ggplot2::ggplot(data = qq.df, ggplot2::aes(x = theo, y = sample)) +
        ggplot2::stat_qq_line(data = qq.df, ggplot2::aes(sample = sample), colour = "steelblue", linewidth = 0.75, inherit.aes = F) +
        ggplot2::geom_point(size = 2) +
        ggplot2::labs(y = "Sample quantiles", x = "Theoretical quantiles") +
        ggplot2::theme(
          title = ggplot2::element_text(size = 10),
          axis.title.x = ggplot2::element_text(size = 12),
          axis.text.x = ggplot2::element_text(size = 9),
          axis.title.y = ggplot2::element_text(size = 12),
          axis.text.y = ggplot2::element_text(size = 9)
        )
      print(p1)
    }
  }

  if (labels) {
    colnames(df) <- tolower(colnames(df))
    if (any(!c("col", "row") %in% colnames(df))) {
      stop("'df' must contain columns 'col' and 'row' in order to plot labels")
    }
    qq.df <- data.frame(
      col = df[["col"]],
      row = df[["row"]],
      effect = df[[effect]]
    )
    qq.df$col <- factor(as.numeric(trimws(qq.df$col)))
    qq.df$row <- factor(as.numeric(trimws(qq.df$row)))
    p1 <- ggplot2::ggplot(qq.df, ggplot2::aes(sample = effect)) +
      ggplot2::stat_qq()
    qq.df <- data.frame(
      col = qq.df$col[order(qq.df$effect)],
      row = qq.df$row[order(qq.df$effect)],
      sample = ggplot2::ggplot_build(p1)$data[[1]]["sample"],
      theo = ggplot2::ggplot_build(p1)$data[[1]]["theoretical"]
    )
    qq.df <- qq.df[order(qq.df$col, qq.df$row), ]
    rownames(qq.df) <- NULL
    if (!plot) {
      colnames(qq.df)[2] <- "theoretical"
      print(qq.df)
    }
    if (plot) {
      qq.df$ColRowLabel <- paste0(qq.df$col, ":", qq.df$row)
      p1 <- ggplot2::ggplot(data = qq.df, ggplot2::aes(x = theo, y = sample, label = ColRowLabel)) +
        ggplot2::stat_qq_line(data = qq.df, ggplot2::aes(sample = sample), colour = "steelblue", linewidth = 0.75, inherit.aes = F) +
        ggplot2::geom_text(size = 4) +
        ggplot2::labs(y = "Sample quantiles", x = "Theoretical quantiles") +
        ggplot2::ggtitle(label = "Residuals indexed as Col:Row") +
        ggplot2::theme(
          title = ggplot2::element_text(size = 10),
          axis.title.x = ggplot2::element_text(size = 12),
          axis.text.x = ggplot2::element_text(size = 9),
          axis.title.y = ggplot2::element_text(size = 12),
          axis.text.y = ggplot2::element_text(size = 9)
        )
      print(p1)
    }
  }
}

#' Construct a sample variogram
#'
#' @param df A data frame containing the columns "col", "row", and the effect to be plotted.
#' @param effect The name of the effect to be plotted.
#' @param plot When TRUE (default), the sample variogram is displayed graphically.
#'   Otherwise, a data frame is returned.
#' @param min.np The minimum number of pairs that semi-variances will be displayed for.
#'   By default, \code{min.np = 30}.
#'
#' @return Graphic of the sample variogram, where the x- and y- axes display the row and
#'   column displacements and the z-axis displays the semi-variance (variogram ordinates).
#'   When \code{plot = FALSE}, a data frame is returned with the column and row displacements
#'   as well as the estimated semi-variances.
#'
#' @examples
#' # Plot the sample variogram for the spatial error component in the data frame error_df.
#'
#' sample_variogram(
#'   error_df,
#'   effect = "e_spat",
#'   plot = TRUE,
#' )
#' @export
#'
sample_variogram <- function(df,
                             effect,
                             plot = TRUE,
                             min.np = 30) {
  colnames(df) <- tolower(colnames(df))
  if (any(!c("col", "row") %in% colnames(df))) {
    stop("'df' must contain columns 'col' and 'row', and the effect to be plotted.")
  }

  sample.df <- data.frame(
    col = df[["col"]],
    row = df[["row"]],
    effect = df[[effect]]
  )
  sample.df <- sample.df[order(sample.df$col, sample.df$row), ]

  col.dis <- abs(outer(as.numeric(trimws(sample.df$col)), as.numeric(trimws(sample.df$col)), FUN = "-"))
  row.dis <- abs(outer(as.numeric(trimws(sample.df$row)), as.numeric(trimws(sample.df$row)), FUN = "-"))
  var.mat <- outer(sample.df$effect, sample.df$effect, FUN = "-")^2 / 2
  sample.df <- data.frame(
    col.dis = col.dis[upper.tri(col.dis, diag = T)],
    row.dis = row.dis[upper.tri(row.dis, diag = T)],
    semi.var = var.mat[upper.tri(var.mat, diag = T)]
  )
  sample.df <- sample.df[order(sample.df$col.dis, sample.df$row.dis), ]

  sample.df <- data.frame(
    col.dis = rep(unique(sample.df$col.dis), each = length(unique(sample.df$row.dis))),
    row.dis = unique(sample.df$row.dis),
    semi.var = c(with(sample.df, tapply(semi.var, list(row.dis, col.dis), function(x) mean(x, na.rm = T)))),
    np = c(with(sample.df, tapply(semi.var, list(row.dis, col.dis), function(x) length(x[!is.na(x)]))))
  )

  if (!plot) {
    sample.df$col.dis <- factor(sample.df$col.dis)
    sample.df$row.dis <- factor(sample.df$row.dis)
    print(sample.df)
  }

  if (plot) {
    lattice::lattice.options(
      layout.heights = list(bottom.padding = list(x = -1), top.padding = list(x = -1.5)),
      layout.widths = list(left.padding = list(x = -1.25), right.padding = list(x = -3))
    )
    graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
    p1 <- lattice::wireframe(semi.var ~ row.dis * col.dis,
      data = sample.df[sample.df$np >= min.np, ], drape = T, colorkey = F, zoom = 0.97, cuts = 30,
      screen = list(z = 30, x = -60, y = 0), aspect = c(1, 0.66),
      scales = list(distance = c(1.2, 1.2, 0.5), arrows = F, cex = 0.7, col = "black"),
      zlab = list(label = paste("Semi-variance"), cex = 0.9, rot = 90, just = c(0.5, -2.25)),
      xlab = list(label = paste("Row displacement"), cex = 0.9, rot = 19, just = c(0.5, -0.75)),
      ylab = list(label = paste("Column displacement"), cex = 0.9, rot = -49, just = c(0.5, -0.75)),
      par.settings = list(axis.line = list(col = "transparent"), clip = list(panel = "off"))
    )
    print(p1)
  }
}

#' Construct a theoretical variogram
#'
#' @param n_cols The total number of columns.
#' @param n_rows The total number of rows.
#' @param var_R The total error variance. By default, \code{var_R = 1}.
#' @param prop_spatial The proportion of spatial error variance to total error
#'   variance (spatial + random). By default, \code{prop_spatial = 0.5}.
#' @param col_cor The column autocorrelation value.
#' @param row_cor The row autocorrelation value.
#' @param plot When TRUE (default), the theoretical variogram is displayed graphically.
#'   Otherwise, a data frame is returned.
#'
#' @return Graphic of the theoretical variogram, where the x- and y- axes display the row and
#'   column displacements and the z-axis displays the semi-variance (variogram ordinates).
#'   When \code{plot = FALSE}, a data frame is returned with the column and row displacements
#'   as well as the theoretical semi-variances.
#'
#' @examples
#' # Plot a theoretical variogram for a field with 10 columns and 20 rows,
#' # using column and row autocorrelations of 0.4 and 0.8.
#'
#' theoretical_variogram(
#'   n_cols = 10,
#'   n_rows = 20,
#'   var_R = 1,
#'   prop_spatial = 0.5,
#'   col_cor = 0.4,
#'   row_cor = 0.8,
#'   plot = TRUE
#' )
#' @export
#'
theoretical_variogram <- function(n_cols,
                                  n_rows,
                                  var_R = 1,
                                  prop_spatial = 0.5,
                                  col_cor,
                                  row_cor,
                                  plot = TRUE) {
  prop_rand <- 1 - prop_spatial
  col.displacement <- rep(0:(n_cols - 1), each = n_rows)
  row.displacement <- rep(0:(n_rows - 1), times = n_cols)
  theoretical.df <- data.frame(
    col.dis = col.displacement,
    row.dis = row.displacement,
    semi.var = var_R * (prop_rand + prop_spatial * (1 - col_cor^(col.displacement) * row_cor^(row.displacement)))
  )
  theoretical.df$semi.var[1] <- 0
  theoretical.df$col.dis <- as.numeric(trimws(theoretical.df$col.dis))
  theoretical.df$row.dis <- as.numeric(trimws(theoretical.df$row.dis))

  if (!plot) {
    theoretical.df$col.dis <- factor(theoretical.df$col.dis)
    theoretical.df$row.dis <- factor(theoretical.df$row.dis)
    print(theoretical.df)
  }

  if (plot) {
    require(lattice)
    lattice::lattice.options(
      layout.heights = list(bottom.padding = list(x = -1), top.padding = list(x = -1.5)),
      layout.widths = list(left.padding = list(x = -1.25), right.padding = list(x = -3))
    )
    graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
    p1 <- lattice::wireframe(semi.var ~ row.dis * col.dis,
      data = theoretical.df, drape = T, colorkey = F, zoom = 0.97, cuts = 30,
      screen = list(z = 30, x = -60, y = 0), aspect = c(1, 0.66),
      scales = list(distance = c(1.2, 1.2, 0.5), arrows = F, cex = 0.7, col = "black"),
      zlab = list(label = paste("Semi-variance"), cex = 0.9, rot = 90, just = c(0.5, -2.25)),
      xlab = list(label = paste("Row displacement"), cex = 0.9, rot = 19, just = c(0.5, -0.75)),
      ylab = list(label = paste("Column displacement"), cex = 0.9, rot = -49, just = c(0.5, -0.75)),
      par.settings = list(axis.line = list(col = "transparent"), clip = list(panel = "off"))
    )
    print(p1)
  }
}
