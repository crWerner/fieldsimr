#' Graphics for plot effects
#'
#' Creates a graphical field array for a set of plot effects (e.g., phenotypes, genetic values, or plot errors).
#' Requires a data frame generated with the functions \link[FieldSimR]{field_trial_error} and \link[FieldSimR]{make_phenotypes},
#' or any data frame matching the description below.
#'
#' @param df A data frame with the columns 'col', 'row', and the effects to be plotted.
#' @param effect The name of the effects to be plotted.
#' @param blocks When \code{TRUE} (default), the field array is split into blocks.
#'   This requires an additional column 'block' in the data frame.
#' @param labels When \code{TRUE} (default is \code{FALSE}), column and row labels are displayed.
#'
#' @return A graphical field array with x- and y-axes displaying the column and row numbers,
#'  and colour gradient ranging from red (low value) to green (high value).
#'
#' @examples
#' # Display the simulated plot errors in the example data frame 'error_df_bivar' for Trait 1 in Environment 1.
#'
#' error_df <- error_df_bivar[error_df_bivar$env == 1, ]
#'
#' plot_effects(df = error_df,
#'              effect = "e.Trait1", labels = T)
#'
#' @export
plot_effects <- function(df,
                         effect,
                         blocks = TRUE,
                         labels = FALSE) {
  colnames(df)[grep("block|col|row", tolower(colnames(df)))] <- tolower(colnames(df))[grep("block|col|row", tolower(colnames(df)))]

  colnames(df)[grep("col", colnames(df))] <- "col"
  if (any(!c("col", "row", effect) %in% colnames(df))) {
    stop("'df' must contain the columns 'col', 'row', and the effect to be plotted")
  }
  df$col <- factor(as.numeric(as.character(df$col)))
  df$row <- factor(as.numeric(as.character(df$row)))

  if (blocks) {
    if (!"block" %in% colnames(df)) {
      stop("'df' must contain the column 'block' if blocks are to be plotted")
    }
    df$block <- factor(as.numeric(as.character(df$block)))
  }

  colnames(df)[colnames(df) %in% effect] <- "eff"

  ncols <- length(unique(df$col))
  nrows <- length(unique(df$row))
  nblocks <- length(unique(df$block))

  if (nblocks > 1) {
    df1 <- df[df[["block"]] == 1, ]
    df2 <- df[df[["block"]] == 2, ]

    if (any(unique(df[df[["block"]] == 1, ]$row) == unique(df[df[["block"]] == 2, ]$row)) == FALSE) {
      dist <- (nrows / nblocks)
      x_min <- rep(0.5, nblocks)
      y_min <- (seq(0, nrows, dist) + 0.5)[1:nblocks]
      x_max <- rep((ncols + 0.5), nblocks)
      y_max <- (seq(0, nrows, dist) + 0.5)[2:(nblocks + 1)]

      block_x_min <- rep(0.5, nblocks)
      block_y_min <- y_max
      block_x_max <- rep((ncols + 0.5), nblocks)
      block_y_max <- y_max

      block_x_min_2 <- rep(0, nblocks)
      block_y_min_2 <- y_max
      block_x_max_2 <- rep((ncols + 1), nblocks)
      block_y_max_2 <- y_max
    } else if (any(unique(df[df[["block"]] == 1, ]$col) == unique(df[df[["block"]] == 2, ]$col)) == FALSE) {
      dist <- (ncols / nblocks)
      x_min <- (seq(0, ncols, dist) + 0.5)[1:nblocks]
      y_min <- rep(0.5, nblocks)
      x_max <- (seq(0, ncols, dist) + 0.5)[2:(nblocks + 1)]
      y_max <- rep((nrows + 0.5), nblocks)

      block_x_min <- x_max
      block_y_min <- rep(0.5, nblocks)
      block_x_max <- x_max
      block_y_max <- rep((nrows + 0.5), nblocks)

      block_x_min_2 <- x_max
      block_y_min_2 <- rep(0, nblocks)
      block_x_max_2 <- x_max
      block_y_max_2 <- rep((nrows + 1), nblocks)
    } else {
      stop("Check column and row assignment within blocks")
    }
  }

  col <- row <- eff <- NULL
  mid_pt <- mean(df$eff, na.rm = TRUE)
  max_pt <- max(abs(c(mid_pt - min(df$eff, na.rm = TRUE), max(df$eff, na.rm = TRUE) - mid_pt)), na.rm = TRUE) + 1e-8

  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = col, y = row)) +
    ggplot2::geom_tile(ggplot2::aes(fill = eff)) +
    ggplot2::scale_fill_gradient2(
      low = "#A51122", mid = "#FEFDBE", high = "#006228",
      midpoint = mid_pt, limits = c(mid_pt - max_pt, mid_pt + max_pt)
    ) +
    ggplot2::scale_x_discrete(expand = c(0.0001,0.0001)) +
    ggplot2::scale_y_discrete(limits = rev, expand = c(0.0001,0.0001)) +
    ggplot2::xlab("Column") +
    ggplot2::ylab("Row") +
    ggplot2::theme_grey(base_size = 10) +
    ggplot2::ggtitle(effect) +
    ggplot2::labs(fill = "Effect") +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 9),
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 11),
      panel.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(t = 4, r = 0, b = 6, l = 0), size = 12, colour = "gray40")
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = rep(seq(0.5, ncols - 0.5, 1), each = nrows), xmax = rep(seq(1.5, ncols + 0.5, 1), each = nrows),
                   ymin = rep(seq(0.5, nrows - 0.5, 1), ncols), ymax = rep(seq(1.5, nrows + 0.5, 1), ncols)),
      fill = "transparent", colour = "black", linewidth = 0.05, inherit.aes = FALSE
    ) +
    ggplot2::annotate(
      geom = "rect", xmin = 0.5, ymin = 0.5,
      xmax = ncols + 0.5, ymax = nrows + 0.5,
      fill = "transparent", col = "black", lwd = 1.25
    )

  if (labels) {
    p <- p + ggplot2::theme(axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8, r = 0, b = 0, l = 0)),
                            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 7, b = 0, l = 0)))
  } else {
    p <- p + ggplot2::theme(axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6, r = 0, b = 0, l = 0)),
                            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 4, b = 0, l = 0)),
                            axis.ticks = ggplot2::element_blank(),
                            axis.text = ggplot2::element_blank())
  }

  if (blocks == TRUE & length(unique(df$block)) > 1) {
    for (i in 1:(nblocks - 1)) {
      p <- p + ggplot2::geom_segment(
        x = block_x_min[i],
        y = block_y_min[i],
        xend = block_x_max[i],
        yend = block_y_max[i],
        linewidth = 2
      ) +
        ggplot2::geom_segment(
          x = block_x_min_2[i],
          y = block_y_min_2[i],
          xend = block_x_max_2[i],
          yend = block_y_max_2[i],
          linewidth = 0.7,
          col = "white"
        )
    }
  }
  return(p)
}

#' Graphics for matrices
#'
#' Creates a heatmap for a symmetric matrix (e.g., correlation or covariance matrix).
#'
#' @param mat A symmetric matrix.
#' @param order When \code{TRUE} (default is \code{FALSE}), the function \code{agnes} of the package
#'   \href{https://cran.r-project.org/package=cluster}{cluster} is used with default arguments to
#'   order the matrix based on a dendrogram.
#' @param labels When \code{TRUE} (default is \code{FALSE}), variable labels are displayed.
#'
#' @return A graphical field array with x- and y-axes displaying the variable numbers,
#'  and colour gradient ranging from blue (low value) to red (high value).
#'
#' @examples
#' # Display a simulated correlation matrix.
#'
#' corA <- rand_cor_mat(n = 10,
#'                      min.cor = -1,
#'                      max.cor = 1)
#'
#' plot_matrix(mat = corA,
#'             order = TRUE,
#'             labels = TRUE)
#'
#' @export
plot_matrix <- function(mat,
                        order = FALSE,
                        labels = FALSE) {
  mat <- round(mat, 12)
  if (!isSymmetric(mat)) stop("'mat' must be a symmetric matrix")

  n <- ncol(mat)
  if (is.null(colnames(mat)) & !is.null(rownames(mat))) {
    colnames(mat) <- rownames(mat)
  } else if(!is.null(colnames(mat)) & is.null(rownames(mat))) {
    rownames(mat) <- colnames(mat)
  } else {
    colnames(mat) <- rownames(mat) <- 1:n
  }

  is_cor_mat <- TRUE
  effect <- "Correlation matrix"
  effect_short <- "Cor."
  if (any(diag(mat) != 1)) {
    is_cor_mat <- FALSE
    effect <- "Covariance matrix"
    effect_short <- "Cov."
  }

  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("var1", "var2", "eff")

  if(is_cor_mat) {
    df$eff[df$var1 == df$var2] <- NA
  }
  df$var1 <- factor(as.numeric(trimws(df$var1)))
  df$var2 <- factor(as.numeric(trimws(df$var2)))

  if(order) {
    if (!is_cor_mat) {
      mat <- cov2cor(mat)
    }
    dis_mat <- 1 - mat
    order2 <- cluster::agnes(x = dis_mat, diss = TRUE, method = "average")$order
    df$var1 <- factor(df$var1, levels = order2)
    df$var2 <- factor(df$var2, levels = order2)
  }

  var1 <- var2 <- eff <- NULL
  if (is_cor_mat) {
    mid_pt <- 0
    max_pt <- 1.1
  } else {
    mid_pt <- mean(df$eff, na.rm = TRUE)
    max_pt <- max(abs(c(mid_pt - min(df$eff, na.rm = TRUE), max(df$eff, na.rm = TRUE) - mid_pt)), na.rm = TRUE) + 1e-8
  }

  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = var1, y = var2)) +
    ggplot2::geom_tile(ggplot2::aes(fill = eff)) +
    ggplot2::scale_fill_gradient2(
      low = "midnightblue", mid = "#FEFDBE", high = "#A51122", na.value = "transparent",
      midpoint = mid_pt, limits = c(mid_pt - max_pt, mid_pt + max_pt)
    ) +
    ggplot2::scale_x_discrete(expand = c(0.0001,0.0001)) +
    ggplot2::scale_y_discrete(limits = rev, expand = c(0.0001,0.0001)) +
    ggplot2::xlab("Variable") +
    ggplot2::ylab("Variable") +
    ggplot2::theme_grey(base_size = 10) +
    ggplot2::ggtitle(effect) +
    ggplot2::labs(fill = effect_short) +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 9),
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 11),
      panel.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(t = 4, r = 0, b = 6, l = 0), size = 12, colour = "gray40")
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = rep(seq(0.5, n - 0.5, 1), each = n), xmax = rep(seq(1.5, n + 0.5, 1), each = n),
                   ymin = rep(seq(0.5, n - 0.5, 1), n), ymax = rep(seq(1.5, n + 0.5, 1), n)),
      fill = "transparent", colour = "black", linewidth = 0.05, inherit.aes = FALSE
    ) +
    ggplot2::annotate(
      geom = "rect", xmin = 0.5, ymin = 0.5,
      xmax = n + 0.5, ymax = n + 0.5,
      fill = "transparent", col = "black", lwd = 1.25
    )

  if (labels) {
    p <- p + ggplot2::theme(axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8, r = 0, b = 0, l = 0)),
                            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 6, b = 0, l = 0)))
  } else {
    p <- p + ggplot2::theme(axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6, r = 0, b = 0, l = 0)),
                            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 4, b = 0, l = 0)),
                            axis.ticks = ggplot2::element_blank(),
                            axis.text = ggplot2::element_blank())
  }
  return(p)
}

#' Q-Q plot
#'
#' Creates a normal quantile-quantile (Q-Q) plot for a set of effects (e.g., phenotypes, genetic values, or plot errors).
#'
#' @param df A data frame with the effects to be plotted.
#' @param effect The name of the effects to be plotted.
#' @param labels When \code{TRUE} (default is \code{FALSE}), column and row labels are displayed.
#'   This requires additional columns 'col' and 'row' in the data frame.
#'
#' @return A Q-Q plot with x- and y-axes displaying the theoretical and sample quantiles of
#'   the effects, respectively.
#'
#' @examples
#' # Q-Q plot of the simulated plot errors in the example data frame 'error_df_bivar' for Trait 1 in Environment 1.
#'
#' error_df <- error_df_bivar[error_df_bivar$env == 1, ]
#'
#' qq <- qq_plot(df = error_df,
#'               effect = "e.Trait1",
#'               labels = TRUE)
#'
#' # Q-Q plot
#' qq
#'
#' # Extract the data frame with the theoretical and sample quantiles of
#' # the user-defined effects.
#' qq_df <- qq$data
#'
#' @export
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
    colnames(df)[grep("col", colnames(df))] <- "col"
    if (any(!c("col", "row") %in% colnames(df))) {
      stop("'df' must contain the columns 'col' and 'row' if labels are to be plotted")
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

    qq_df$cr.label <- factor(paste0(qq_df$col, ":", qq_df$row))
    theoretical <- cr.label <- NULL
    p <- ggplot2::ggplot(data = qq_df, ggplot2::aes(x = theoretical, y = sample, label = cr.label)) +
      ggplot2::stat_qq_line(data = qq_df, ggplot2::aes(sample = sample), colour = "steelblue", linewidth = 0.75, inherit.aes = F) +
      ggplot2::geom_text(size = 4) +
      ggplot2::labs(y = "Sample quantiles", x = "Theoretical quantiles") +
      ggplot2::ggtitle(label = "Effects indexed as col:row") +
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
#' Creates a sample variogram for a set of effects (e.g., plot errors).
#'
#' @param df A data frame with the columns 'col', 'row', and the effects to be plotted.
#' @param effect The name of the effects to be plotted.
#' @param min.np Minimum number of pairs for which semi-variances are displayed (default is 30).
#'
#' @return A sample variogram with x- and y-axes displaying the row and
#'   column displacements, and the z-axis displaying the average semi-variances (variogram ordinates)
#'   for the effects.
#'
#' @examples
#' # Sample variogram of the simulated plot errors in the example data frame 'error_df_bivar' for Trait 1 in Environment 1.
#'
#' error_df <- error_df_bivar[error_df_bivar$env == 1, ]
#'
#' variogram <- sample_variogram(df = error_df,
#'                               effect = "e.Trait1")
#'
#' # Sample variogram
#' variogram
#'
#' # Extract the data frame with the column and row displacements, and the
#' # sample semi-variances.
#'
#' variogram_df <- variogram$data
#'
#' @export
sample_variogram <- function(df,
                             effect,
                             min.np = 30) {
  colnames(df) <- tolower(colnames(df))
  effect <- tolower(effect)
  colnames(df)[colnames(df) %in% effect] <- "eff"

  colnames(df)[grep("col", colnames(df))] <- "col"
  if (any(!c("col", "row") %in% colnames(df))) {
    stop("'df' must contain the columns 'col' and 'row', and the effect to be plotted")
  }

  variogram_df <- data.frame(
    col = df[["col"]],
    row = df[["row"]],
    effect = df[["eff"]]
  )
  variogram_df <- variogram_df[order(variogram_df$col, variogram_df$row), ]

  col_dis <- abs(outer(as.numeric(trimws(variogram_df$col)), as.numeric(trimws(variogram_df$col)), FUN = "-"))
  row_dis <- abs(outer(as.numeric(trimws(variogram_df$row)), as.numeric(trimws(variogram_df$row)), FUN = "-"))
  var_mat <- outer(variogram_df$effect, variogram_df$effect, FUN = "-")^2 / 2
  variogram_df <- data.frame(
    col_dis = col_dis[upper.tri(col_dis, diag = T)],
    row_dis = row_dis[upper.tri(row_dis, diag = T)],
    semi_var = var_mat[upper.tri(var_mat, diag = T)]
  )
  variogram_df <- variogram_df[order(variogram_df$col_dis, variogram_df$row_dis), ]

  variogram_df <- data.frame(
    col.dis = rep(unique(variogram_df$col_dis), each = length(unique(variogram_df$row_dis))),
    row.dis = unique(variogram_df$row_dis),
    semi.var = c(with(variogram_df, tapply(semi_var, list(row_dis, col_dis), function(x) mean(x, na.rm = T)))),
    np = c(with(variogram_df, tapply(semi_var, list(row_dis, col_dis), function(x) length(x[!is.na(x)]))))
  )

  lattice::lattice.options(
    layout.heights = list(bottom.padding = list(x = -1), top.padding = list(x = -1.5)),
    layout.widths = list(left.padding = list(x = -1.25), right.padding = list(x = -3))
  )
  graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
  p <- lattice::wireframe(semi.var ~ row.dis * col.dis,
    data = variogram_df[variogram_df$np >= min.np, ], drape = T, colorkey = F, zoom = 0.97, cuts = 30,
    screen = list(z = 30, x = -60, y = 0), aspect = c(1, 0.66),
    scales = list(distance = c(1.2, 1.2, 0.5), arrows = F, cex = 0.7, col = "black"),
    zlab = list(label = paste("Semi-variance"), cex = 0.9, rot = 90, just = c(0.5, -2.25)),
    xlab = list(label = paste("Row displacement"), cex = 0.9, rot = 19, just = c(0.5, -0.75)),
    ylab = list(label = paste("Column displacement"), cex = 0.9, rot = -49, just = c(0.5, -0.75)),
    par.settings = list(axis.line = list(col = "transparent"), clip = list(panel = "off"))
  )
  variogram_df$col.dis <- factor(variogram_df$col.dis)
  variogram_df$row.dis <- factor(variogram_df$row.dis)
  p$data <- variogram_df
  return(p)
}

#' Theoretical variogram
#'
#' Creates a theoretical variogram for a separable first order autoregressive process.
#'
#' @param ncols A scalar defining the number of columns.
#' @param nrows A scalar defining the number of rows.
#' @param varR A scalar defining the total error variance.
#' @param prop.spatial A scalar defining the proportion of spatial error variance to total (spatial + random) error
#'   variance.
#' @param col.cor A scalar defining the column autocorrelation,
#' @param row.cor A scalar defining the row autocorrelation.
#'
#' @return A theoretical variogram with x- and y-axes displaying the row and column displacements,
#'  and the z-axis displaying the semi-variances (variogram ordinates) for a separable autoregressive process.
#'
#' @examples
#' # Theoretical variogram for a field trial with 10 columns and 20 rows, based on column and row
#' # autocorrelations of 0.5 and 0.7, and a proportion of spatial error variance of 0.5.
#'
#' variogram <- theoretical_variogram(ncols = 10,
#'                                    nrows = 20,
#'                                    varR = 1,
#'                                    prop.spatial = 0.5,
#'                                    col.cor = 0.5,
#'                                    row.cor = 0.7)
#'
#' # Theoretical variogram
#' variogram
#'
#' # Extract the data frame with the column and row displacements, and the
#' # theoretical semi-variances.
#'
#' variogram_df <- variogram$data
#'
#' @export
theoretical_variogram <- function(ncols = 10,
                                  nrows = 10,
                                  varR = 1,
                                  prop.spatial = 1,
                                  col.cor = 0.5,
                                  row.cor = 0.5) {
  prop_rand <- 1 - prop.spatial
  col_dis <- rep(0:(ncols - 1), each = nrows)
  row_dis <- rep(0:(nrows - 1), times = ncols)
  variogram_df <- data.frame(
    col.dis = col_dis,
    row.dis = row_dis,
    semi.var = varR * (prop_rand + prop.spatial * (1 - col.cor^(col_dis) * row.cor^(row_dis)))
  )
  variogram_df$semi.var[1] <- 0
  variogram_df$col.dis <- as.numeric(trimws(variogram_df$col.dis))
  variogram_df$row.dis <- as.numeric(trimws(variogram_df$row.dis))

  lattice::lattice.options(
    layout.heights = list(bottom.padding = list(x = -1), top.padding = list(x = -1.5)),
    layout.widths = list(left.padding = list(x = -1.25), right.padding = list(x = -3))
  )
  graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
  p <- lattice::wireframe(semi.var ~ row.dis * col.dis,
    data = variogram_df, drape = T, colorkey = F, zoom = 0.97, cuts = 30,
    screen = list(z = 30, x = -60, y = 0), aspect = c(1, 0.66),
    scales = list(distance = c(1.2, 1.2, 0.5), arrows = F, cex = 0.7, col = "black"),
    zlab = list(label = paste("Semi-variance"), cex = 0.9, rot = 90, just = c(0.5, -2.25)),
    xlab = list(label = paste("Row displacement"), cex = 0.9, rot = 19, just = c(0.5, -0.75)),
    ylab = list(label = paste("Column displacement"), cex = 0.9, rot = -49, just = c(0.5, -0.75)),
    par.settings = list(axis.line = list(col = "transparent"), clip = list(panel = "off"))
  )
  variogram_df$col.dis <- factor(variogram_df$col.dis)
  variogram_df$row.dis <- factor(variogram_df$row.dis)
  p$data <- variogram_df
  return(p)
}
