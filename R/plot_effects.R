#' Graphics for plot effects
#'
#' Graphically displays plot effects (e.g., phenotypes, genetic values, plot errors) onto a
#' field array, in which the colour gradient ranges from red (low value) to green (high value).
#' The function requires a data frame generated with field_trial_error as an input, or any data
#' frame with columns 'col', 'row', and the effect to be displayed. When the data frame contains a
#' 'block' column, the field array is split into blocks if \code{blocks = TRUE}.
#'
#' @param df A data frame with columns 'col', 'row', and the effect to be plotted. When \code{df}
#'   also contains a 'block' column, the field array is split into blocks if \code{blocks = TRUE}.
#' @param effect The effect to be plotted.
#' @param blocks When \code{TRUE} (default), the field array is split into blocks.
#' @param labels When \code{TRUE}, row and column numbers are plotted. By default,
#'   \code{labels = FALSE}.
#'
#'
#' @return A graphical field array, in which the colour gradient ranges from red (low value)
#' to green (high value).
#'
#' @examples
#' # Plot the simulated plot errors for trait 2 in environment 2 provided in the example data
#' # frame 'df_error_bivar'.
#'
#' error_df <- df_error_bivar[df_error_bivar$env == 2, ]
#'
#' plot_effects(
#'   df = error_df,
#'   effect = "e.Trait.2"
#' )
#' @export
plot_effects <- function(df,
                         effect,
                         blocks = TRUE,
                         labels = FALSE) {
  colnames(df) <- tolower(colnames(df))
  effect <- tolower(effect)

  if (any(!c("col", "row", effect) %in% colnames(df))) {
    stop("'df' must contain columns 'col', 'row', and the effect to be plotted.")
  }

  if (blocks) {
    if (!"block" %in% colnames(df)) {
      stop("'df' must contain column 'block' if blocks = TRUE.")
    }
  }

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
  mid_pt <- mean(df$eff)
  max_pt <- max(abs(df$eff))

  if (labels) {
    p <- ggplot2::ggplot(data = df, ggplot2::aes(x = col, y = row)) +
      ggplot2::geom_tile(ggplot2::aes(fill = eff)) +
      ggplot2::scale_fill_gradient2(
        low = "#A51122", mid = "#FEFDBE", high = "#006228",
        midpoint = mid_pt, limits=c(-max_pt, max_pt)
      ) +
      ggplot2::scale_y_discrete(limits = rev) +
      ggplot2::xlab("Column") +
      ggplot2::ylab("Row") +
      ggplot2::theme_grey(base_size = 10) +
      ggplot2::ggtitle(effect) +
      ggplot2::labs(fill = "Effect") +
      ggplot2::theme(
        axis.ticks = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(size = 11),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8, r = 0, b = 0, l = 0)),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 8, b = 0, l = 0)),
        axis.title = ggplot2::element_text(size = 12),
        panel.background = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 12, colour = "gray40")
      ) +
      ggplot2::annotate(
        geom = "rect", xmin = 0.5, ymin = 0.5,
        xmax = n_cols + 0.5, ymax = n_rows + 0.5,
        fill = "transparent", col = "black", lwd = 0.5
      )
  } else {
    p <- ggplot2::ggplot(data = df, ggplot2::aes(x = col, y = row)) +
      ggplot2::geom_tile(ggplot2::aes(fill = eff)) +
      ggplot2::scale_fill_gradient2(
        low = "#A51122", mid = "#FEFDBE", high = "#006228",
        midpoint = mid_pt, limits=c(-max_pt, max_pt)
      ) +
      ggplot2::scale_y_discrete(limits = rev) +
      ggplot2::xlab("Column") +
      ggplot2::ylab("Row") +
      ggplot2::theme_grey(base_size = 10) +
      ggplot2::ggtitle(effect) +
      ggplot2::labs(fill = "Effect") +
      ggplot2::theme(
        axis.ticks = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(size = 12),
        panel.background = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 12, colour = "gray40")
      ) +
      ggplot2::annotate(
        geom = "rect", xmin = 0.5, ymin = 0.5,
        xmax = n_cols + 0.5, ymax = n_rows + 0.5,
        fill = "transparent", col = "black", lwd = 0.5
      )
  }

  if (blocks == TRUE & length(unique(df$block)) > 1) {
    for (i in 1:(n_blocks - 1)) {
      p <- p + ggplot2::geom_segment(
        x = block_x_min[i],
        y = block_y_min[i],
        xend = block_x_max[i],
        yend = block_y_max[i],
        linewidth = 1.5
      ) +
        ggplot2::geom_segment(
          x = block_x_min_2[i],
          y = block_y_min_2[i],
          xend = block_x_max_2[i],
          yend = block_y_max_2[i],
          linewidth = 1,
          col = "white"
        )
    }
  }
  return(p)
}

#' Q-Q plot
#'
#' Creates a quantile-quantile (Q-Q) plot which compares the theoretical quantiles of a normal
#' distribution with the sample quantiles of the distribution of user effects.
#'
#' @param df A data frame containing the effect to be plotted.
#' @param effect The name of the effect to be plotted.
#' @param labels When \code{FALSE} (default), data points without labels are plotted. When
#'   \code{TRUE}, column and row labels are shown in the Q-Q plot. This requires additional
#'   columns 'col' and 'row' in the data frame.
#'
#' @return A Q-Q plot with the x- and y-axes displaying the theoretical and sample quantiles of
#'   the effect to be plotted, respectively.
#'
#' @examples
#' # Q-Q plot of the simulated plot errors for trait 2 in environment 2 provided in the example
#' # data frame 'df_error_bivar'.
#'
#' error_df <- df_error_bivar[df_error_bivar$env == 2, ]
#'
#' qq <- qq_plot(
#'   df = error_df,
#'   effect = "e.Trait.2",
#'   labels = TRUE
#' )
#'
#' # Q-Q plot
#' qq
#'
#' # Extraction of a data frame containing the theoretical and sample quantiles of
#' # the effect to be plotted.
#' qq_df <- qq$data
#'
#' @export
#'
#'

qq_plot <- function(df,
                    effect,
                    labels = FALSE) {
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

#' Sample variogram
#'
#' Creates a sample variogram. The x- and y-axes display the row and column displacements,
#' respectively. The z-axis displays the semi-variance (variogram ordinates).
#'
#' @param df A data frame containing the columns 'col', 'row', and the effect to be plotted.
#' @param effect The name of the effect to be plotted.
#' @param min_np Only semi variances based on at least \code{min_np} pairs of plots will be
#'   displayed. By default, \code{min_np = 30}.
#'
#' @return Graphic of the sample variogram, where the x- and y-axes display the row and
#'   column displacements and the z-axis displays the semi-variance (variogram ordinates).
#'
#' @examples
#' # Sample variogram of the simulated plot errors for trait 2 in environment 2 provided in the
#' # example data frame 'df_error_bivar'.
#'
#' error_df <- df_error_bivar[df_error_bivar$env == 2, ]
#'
#' vario <- sample_variogram(
#'   df = error_df,
#'   effect = "e.Trait.2",
#' )
#'
#' # Sample variogram
#' vario
#'
#' # Extraction of a data frame containing the column and row displacements as well as the
#' # semi-variances (sample variogram ordinates).
#'
#' sample_df <- vario$data
#'
#' @export
#'
sample_variogram <- function(df,
                             effect,
                             min_np = 30) {
  colnames(df) <- tolower(colnames(df))
  effect <- tolower(effect)
  colnames(df)[colnames(df) %in% effect] <- "eff"

  if (any(!c("col", "row") %in% colnames(df))) {
    stop("'df' must contain columns 'col' and 'row', and the effect to be plotted.")
  }

  sample_df <- data.frame(
    col = df[["col"]],
    row = df[["row"]],
    effect = df[["eff"]]
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
  sample_df$col_dis <- factor(sample_df$col_dis)
  sample_df$row_dis <- factor(sample_df$row_dis)
  p$data <- sample_df
  return(p)
}

#' Theoretical variogram
#'
#' Creates a theoretical variogram. The x- and y-axes display the row and column displacements,
#' respectively. The z-axis displays the semi-variance (variogram ordinates).
#'
#' @param n_cols A scalar defining the number of columns.
#' @param n_rows A scalar defining the number of rows.
#' @param var_R A scalar defining the total error variance.
#' @param prop_spatial A scalar defining the proportion of spatial error variance to total error
#'   variance (spatial + random).
#' @param col_cor A scalar defining the column autocorrelation value.
#' @param row_cor A scalar defining the row autocorrelation value.
#'
#' @return Graphic of the theoretical variogram, where the x- and y- axes display the row and
#'   column displacements and the z-axis displays the semi-variance (variogram ordinates).
#'
#' @examples
#' # Theoretical variogram for a field with 10 columns and 20 rows, using column and row
#' # autocorrelations of 0.4 and 0.8.
#'
#' vario <- theoretical_variogram(
#'   n_cols = 10,
#'   n_rows = 20,
#'   var_R = 1,
#'   prop_spatial = 0.5,
#'   col_cor = 0.4,
#'   row_cor = 0.8
#' )
#'
#' # Theoretical variogram
#' vario
#'
#' # Extraction of a data frame containing the column and row displacements as well as the
#' # semi-variances (theoretical variogram ordinates).
#'
#' theoretical_df <- vario$data
#'
#' @export
#'
theoretical_variogram <- function(n_cols,
                                  n_rows,
                                  var_R = 1,
                                  prop_spatial = 0.5,
                                  col_cor,
                                  row_cor) {
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
  theoretical_df$col_dis <- factor(theoretical_df$col_dis)
  theoretical_df$row_dis <- factor(theoretical_df$row_dis)
  p$data <- theoretical_df
  return(p)
}
