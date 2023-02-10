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
#' error_df <- df_error_bivar[df_error_bivar$env == 2,]
#'
#' plot_effects(
#'   df = error_df,
#'   effect = "e.Trait.2"
#' )
#' @export
plot_effects <- function(df,
                         effect,
                         blocks = TRUE) {
  if (inherits(df, "list")) df <- data.frame(df[[1]])

  colnames(df) <- tolower(colnames(df))
  if (any(!c("col", "row") %in% colnames(df))) {
    stop("'df' must contain columns 'col', 'row', and the effect to be plotted.")
  }
  
  effect <- tolower(effect)
  colnames(df)[colnames(df) %in% effect] <- "eff"

  n_cols <- length(unique(df$col))
  n_rows <- length(unique(df$row))
  n_blocks <- length(unique(df$block))

  if (n_blocks > 1) {
    df1 <- df[df[["block"]] == 1, ]
    df2 <- df[df[["block"]] == 2, ]

    if (any(unique(df[df[["block"]] == 1, ]$row) == unique(df[df[["block"]] == 2, ]$row)) == FALSE) {
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
    } else if (any(unique(df[df[["block"]] == 1, ]$col) == unique(df[df[["block"]] == 2, ]$col)) == FALSE) {
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

  col <- row <- eff <- NULL

  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = col, y = row)) +
    ggplot2::geom_tile(ggplot2::aes(fill = eff)) +
    ggplot2::scale_fill_gradient2(low = "#A51122", mid = "#FEFDBE", high = "#006228") +
    ggplot2::xlab("Column") +
    ggplot2::ylab("Row") +
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
#' @param df A data frame containing the effect to be plotted.
#' @param effect The name of the effect to be plotted.
#' @param labels When FALSE (default) data points without labels are plotted. When TRUE, 
#'   column and row labels are inserted onto the qqplot. This requires the additional columns 
#'   "col" and "row" in the data frame.
#' @param plot When TRUE (default), the qqplot is displayed graphically. When FALSE, a 
#'   data frame is returned.
#'
#' @return Graphic of the qqplot, where the x- and y- axes display the theoretical and
#'   sample quantiles. When \code{plot = FALSE}, a data frame is returned with the "theoretical" 
#'   and "sample" quantiles, as well as the columns "col and "row" when \code{labels=TRUE}.
#'
#' @examples
#' # Plot the simulated total error for trait 2 in environment 2 provided in the example data
#' # frame 'df_error_bivar'.
#'
#' error_df <- df_error_bivar[df_error_bivar$env == 2,]
#' 
#' qq_plot(
#'   df = error_df,
#'   effect = "e.Trait.2"
#'   labels = TRUE,
#'   plot = TRUE,
#' )
#' @export
#'
#'

qq_plot <- function(df,
                    effect,
                    labels = FALSE,
                    plot = TRUE) {
   colnames(df) <- tolower(colnames(df))
   effect <- tolower(effect)
   colnames(df)[colnames(df) %in% effect] <- "eff"
  
  if (!labels) {
    qq_df <- data.frame(effect = df[["eff"]])
    p <- ggplot2::ggplot(qq_df, ggplot2::aes(sample = effect)) +
      ggplot2::stat_qq()
    qq_df <- data.frame(
      ggplot2::ggplot_build(p)$data[[1]]["sample"],
      ggplot2::ggplot_build(p)$data[[1]]["theoretical"]
    )
    if (!plot) {
      return(qq_df)
    }
    if (plot) {
      p <- ggplot2::ggplot(data = qq_df, ggplot2::aes(x = theoretical, y = sample)) +
        ggplot2::stat_qq_line(data = qq_df, ggplot2::aes(sample = sample), colour = "steelblue", linewidth = 0.75, inherit.aes = F) +
        ggplot2::geom_point(size = 2) +
        ggplot2::labs(y = "Sample quantiles", x = "Theoretical quantiles") +
        ggplot2::theme(
          title = ggplot2::element_text(size = 10),
          axis.title.x = ggplot2::element_text(size = 12),
          axis.text.x = ggplot2::element_text(size = 9),
          axis.title.y = ggplot2::element_text(size = 12),
          axis.text.y = ggplot2::element_text(size = 9)
        )
      return(p)
    }
  }

  if (labels) {
    if (any(!c("col", "row") %in% colnames(df))) {
      stop("'df' must contain columns 'col' and 'row' in order to plot labels")
    }
    qq_df <- data.frame(
      col = df[["col"]],
      row = df[["row"]],
      effect = df[["eff"]]
    )
    qq_df$col <- factor(as.numeric(trimws(qq_df$col)))
    qq_df$row <- factor(as.numeric(trimws(qq_df$row)))
    p <- ggplot2::ggplot(qq_df, ggplot2::aes(sample = effect)) +
      ggplot2::stat_qq()
    qq_df <- data.frame(
      col = qq_df$col[order(qq_df$effect)],
      row = qq_df$row[order(qq_df$effect)],
      ggplot2::ggplot_build(p)$data[[1]]["sample"],
      ggplot2::ggplot_build(p)$data[[1]]["theoretical"]
    )
    qq_df <- qq_df[order(qq_df$col, qq_df$row), ]
    rownames(qq_df) <- NULL
    if (!plot) {
      return(qq_df)
    }
    if (plot) {
      qq_df$ColRowLabel <- paste0(qq_df$col, ":", qq_df$row)
      theoretical <- ColRowLabel <- NULL
      p <- ggplot2::ggplot(data = qq_df, ggplot2::aes(x = theoretical, y = sample, label = ColRowLabel)) +
        ggplot2::stat_qq_line(data = qq_df, ggplot2::aes(sample = sample), colour = "steelblue", linewidth = 0.75, inherit.aes = F) +
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
      return(p)
    }
  }
}

#' Construct a sample variogram
#'
#' @param df A data frame containing the columns "col", "row", and the effect to be plotted.
#' @param effect The name of the effect to be plotted.
#' @param plot When TRUE (default), the sample variogram is displayed graphically.
#'   Otherwise, a data frame is returned.
#' @param min_np Only semi variances based on at least \code{min_np} pairs of plots will be displayed.
#'   By default, \code{min_np = 30}.
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
                             min_np = 30) {
  colnames(df) <- tolower(colnames(df)) # align with plot_effects
  if (any(!c("col", "row") %in% colnames(df))) {
    stop("'df' must contain columns 'col' and 'row', and the effect to be plotted.")
  }

  sample_df <- data.frame(
    col = df[["col"]],
    row = df[["row"]],
    effect = df[[effect]]
  )
  sample_df <- sample_df[order(sample_df$col, sample_df$row), ]

  col_dis <- abs(outer(as.numeric(trimws(sample_df$col)), as.numeric(trimws(sample_df$col)), FUN = "-"))
  row_dis <- abs(outer(as.numeric(trimws(sample_df$row)), as.numeric(trimws(sample_df$row)), FUN = "-"))
  var_mat <- outer(sample_df$effect, sample_df$effect, FUN = "-")^2 / 2
  sample_df <- data.frame(
    col_dis = col_dis[upper.tri(col_dis, diag = T)],
    row_dis = row_dis[upper.tri(row_dis, diag = T)],
    semi_var = var_mat[upper.tri(var_mat, diag = T)]
  )
  sample_df <- sample_df[order(sample_df$col_dis, sample_df$row_dis), ]

  sample_df <- data.frame(
    col_dis = rep(unique(sample_df$col_dis), each = length(unique(sample_df$row_dis))),
    row_dis = unique(sample_df$row_dis),
    semi_var = c(with(sample_df, tapply(semi_var, list(row_dis, col_dis), function(x) mean(x, na.rm = T)))),
    np = c(with(sample_df, tapply(semi_var, list(row_dis, col_dis), function(x) length(x[!is.na(x)]))))
  )

  if (!plot) {
    sample_df$col_dis <- factor(sample_df$col_dis)
    sample_df$row_dis <- factor(sample_df$row_dis)
    return(sample_df)
  }

  if (plot) {
    lattice::lattice.options(
      layout.heights = list(bottom.padding = list(x = -1), top.padding = list(x = -1.5)),
      layout.widths = list(left.padding = list(x = -1.25), right.padding = list(x = -3))
    )
    graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
    p <- lattice::wireframe(semi_var ~ row_dis * col_dis,
      data = sample_df[sample_df$np >= min_np, ], drape = T, colorkey = F, zoom = 0.97, cuts = 30,
      screen = list(z = 30, x = -60, y = 0), aspect = c(1, 0.66),
      scales = list(distance = c(1.2, 1.2, 0.5), arrows = F, cex = 0.7, col = "black"),
      zlab = list(label = paste("Semi-variance"), cex = 0.9, rot = 90, just = c(0.5, -2.25)),
      xlab = list(label = paste("Row displacement"), cex = 0.9, rot = 19, just = c(0.5, -0.75)),
      ylab = list(label = paste("Column displacement"), cex = 0.9, rot = -49, just = c(0.5, -0.75)),
      par.settings = list(axis.line = list(col = "transparent"), clip = list(panel = "off"))
    )
    return(p)
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
  col_displacement <- rep(0:(n_cols - 1), each = n_rows)
  row_displacement <- rep(0:(n_rows - 1), times = n_cols)
  theoretical_df <- data.frame(
    col_dis = col_displacement,
    row_dis = row_displacement,
    semi_var = var_R * (prop_rand + prop_spatial * (1 - col_cor^(col_displacement) * row_cor^(row_displacement)))
  )
  theoretical_df$semi_var[1] <- 0
  theoretical_df$col_dis <- as.numeric(trimws(theoretical_df$col_dis))
  theoretical_df$row_dis <- as.numeric(trimws(theoretical_df$row_dis))

  if (!plot) {
    theoretical_df$col_dis <- factor(theoretical_df$col_dis)
    theoretical_df$row_dis <- factor(theoretical_df$row_dis)
    return(theoretical_df)
  }

  if (plot) {
    lattice::lattice.options(
      layout.heights = list(bottom.padding = list(x = -1), top.padding = list(x = -1.5)),
      layout.widths = list(left.padding = list(x = -1.25), right.padding = list(x = -3))
    )
    graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
    p <- lattice::wireframe(semi_var ~ row_dis * col_dis,
      data = theoretical_df, drape = T, colorkey = F, zoom = 0.97, cuts = 30,
      screen = list(z = 30, x = -60, y = 0), aspect = c(1, 0.66),
      scales = list(distance = c(1.2, 1.2, 0.5), arrows = F, cex = 0.7, col = "black"),
      zlab = list(label = paste("Semi-variance"), cex = 0.9, rot = 90, just = c(0.5, -2.25)),
      xlab = list(label = paste("Row displacement"), cex = 0.9, rot = 19, just = c(0.5, -0.75)),
      ylab = list(label = paste("Column displacement"), cex = 0.9, rot = -49, just = c(0.5, -0.75)),
      par.settings = list(axis.line = list(col = "transparent"), clip = list(panel = "off"))
    )
    return(p)
  }
}
