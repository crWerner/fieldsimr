#' Graphics for plot effects
#'
#' Creates a graphical field array for a set of plot effects (e.g., phenotypes, genetic values, or plot errors).
#' Requires a data frame generated with the functions \link[FieldSimR]{field_trial_error} or \link[FieldSimR]{make_phenotypes},
#' or any data frame matching the description below.
#'
#' @param df A data frame with the columns 'col', 'row', and the effects to be plotted.
#' @param effect The name of the effects to be plotted.
#' @param blocks When \code{TRUE} (default), the field array is split into blocks.
#'   This requires an additional column 'block' in the data frame.
#' @param labels When \code{TRUE} (default), column and row labels are displayed.
#'
#' @return A graphical field array with x- and y-axes displaying the column and row numbers,
#'  and colour gradient ranging from red (low value) to green (high value).
#'
#' @examples
#' # Display the simulated plot errors in the example data frame 'error_df_bivar'
#' # for Trait 1 in Environment 1.
#'
#' error_df <- error_df_bivar[error_df_bivar$env == 1, ]
#'
#' plot_effects(
#'   df = error_df,
#'   effect = "e.Trait1",
#'   labels = TRUE
#' )
#'
#' @export
plot_effects <- function(df,
                         effect,
                         blocks = TRUE,
                         labels = TRUE) {
  if (!is.data.frame(df)) {
    stop("'df' must be a data frame")
  }
  colnames(df)[grep("block|col|row", tolower(colnames(df)))] <- tolower(colnames(df))[grep("block|col|row", tolower(colnames(df)))]

  if (any(!c("col", "row", effect) %in% colnames(df))) {
    stop("'df' must contain the columns 'col', 'row', and the effect to be plotted")
  }
  df$col <- factor(as.numeric(as.character(df$col)))
  df$row <- factor(as.numeric(as.character(df$row)))

  nblocks <- 1
  if (blocks) {
    if (!"block" %in% colnames(df)) {
      stop("'df' must contain the column 'block' if blocks are to be plotted")
    }
    df$block <- factor(as.numeric(as.character(df$block)))
  }

  ncols <- nlevels(df$col)
  nrows <- nlevels(df$row)
  nblocks <- nlevels(df$block)

  plot_x_min <- rep(seq(0.5, ncols - 0.5, 1), each = nrows)
  plot_y_min <- rep(seq(0.5, nrows - 0.5, 1), ncols)
  plot_x_max <- rep(seq(1.5, ncols + 0.5, 1), each = nrows)
  plot_y_max <- rep(seq(1.5, nrows + 0.5, 1), ncols)

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

  col <- row <- NULL
  if (effect != "block") {
    mid_pt <- mean(df[[effect]], na.rm = TRUE)
    max_pt <- max(abs(c(mid_pt - min(df[[effect]], na.rm = TRUE), max(df[[effect]], na.rm = TRUE) - mid_pt)), na.rm = TRUE) + 1e-8
    p <- ggplot2::ggplot(data = df, ggplot2::aes(x = col, y = row)) +
      ggplot2::geom_tile(ggplot2::aes(fill = get(effect))) +
      ggplot2::scale_fill_gradient2(
        low = "#A51122", mid = "#FEFDBE", high = "#006228",
        midpoint = mid_pt, limits = c(mid_pt - max_pt, mid_pt + max_pt)
      ) +
      ggplot2::ggtitle(label = effect) +
      ggplot2::labs(fill = "Effect")
  } else if (effect == "block") {
    p <- ggplot2::ggplot(data = df, ggplot2::aes(x = col, y = row)) +
      ggplot2::geom_tile(ggplot2::aes(fill = get(effect)), alpha = 0.6) +
      ggplot2::scale_fill_manual(values = c("#888888", "#6699CC", "#882255", "#117733", "#332288")) +
      ggplot2::labs(fill = "Block")
  }
  p <- p + ggplot2::scale_x_discrete(expand = c(0.0001, 0.0001)) +
    ggplot2::scale_y_discrete(limits = rev, expand = c(0.0001, 0.0001)) +
    ggplot2::xlab("Column") +
    ggplot2::ylab("Row") +
    ggplot2::theme_grey(base_size = 10) +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 9),
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 11),
      panel.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(t = 4, r = 0, b = 6, l = 0), size = 12, colour = "gray40")
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = plot_x_min, xmax = plot_x_max,
        ymin = plot_y_min, ymax = plot_y_max
      ),
      fill = "transparent", colour = "black", linewidth = 0.05, inherit.aes = FALSE
    ) +
    ggplot2::annotate(
      geom = "rect", xmin = 0.5, ymin = 0.5,
      xmax = ncols + 0.5, ymax = nrows + 0.5,
      fill = "transparent", col = "black", lwd = 1.5
    )

  if (labels) {
    p <- p + ggplot2::theme(
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8, r = 0, b = 0, l = 0)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 7, b = 0, l = 0))
    )
  } else {
    p <- p + ggplot2::theme(
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6, r = 0, b = 0, l = 0)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 4, b = 0, l = 0)),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank()
    )
  }

  if (nblocks > 1) {
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
#' Creates a heatmap of a symmetric matrix (e.g., correlation or covariance matrix).
#'
#' @param mat A symmetric matrix.
#' @param order When \code{TRUE} (default is \code{FALSE}), the function \code{agnes} of the R package
#'   \href{https://cran.r-project.org/package=cluster}{`cluster`} is used with default arguments to
#'   order the matrix based on a dendrogram.
#' @param group.df An optional data frame with columns containing the variable names followed by the group numbers.
#'   When supplied, the heatmap is split into groups and then ordered (when \code{order = TRUE}).
#' @param labels When \code{TRUE} (default), variable labels are displayed.
#'
#' @return A heatmap with x- and y-axes displaying the variable numbers,
#'  and colour gradient ranging from blue (low value) to red (high value).
#'
#' @examples
#' # Display a random correlation matrix.
#'
#' corA <- rand_cor_mat(
#'   n = 10,
#'   min.cor = -1,
#'   max.cor = 1
#' )
#'
#' # Define groups.
#'
#' group_df <- data.frame(variable = 1:10, group = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4))
#'
#' plot_matrix(
#'   mat = corA,
#'   order = TRUE,
#'   group.df = group_df,
#'   labels = TRUE
#' )
#'
#' @export
plot_matrix <- function(mat,
                        order = FALSE,
                        group.df = NULL,
                        labels = TRUE) {
  if (!is.matrix(mat)) stop("'mat' must be a matrix")
  if (!isSymmetric(mat)) stop("'mat' must be a symmetric matrix")
  if (any(colnames(mat) != rownames(mat))) stop("colnames and rownames of 'mat' must match")

  n <- ncol(mat)
  if (is.null(colnames(mat)) && !is.null(rownames(mat))) {
    colnames(mat) <- rownames(mat)
  } else if (!is.null(colnames(mat)) && is.null(rownames(mat))) {
    rownames(mat) <- colnames(mat)
  } else {
    colnames(mat) <- rownames(mat) <- 1:n
  }

  groups <- FALSE
  ngroups <- 1
  if (!is.null(group.df)) {
    if (!is.data.frame(group.df)) stop("'group.df' must be a data frame")
    if (ncol(group.df) < 2) stop("'group.df' must be a data frame with at least two columns")
    colnames(group.df)[1:2] <- c("variable", "group")
    if (any(!colnames(mat) %in% group.df$variable)) stop("'group.df' must contain all variables in 'mat'")

    group.df$variable <- factor(as.numeric(as.character(group.df$variable)))
    group.df$group <- factor(as.numeric(as.character(group.df$group)))
    group.df <- group.df[order(group.df$group), ]
    rownames(group.df) <- NULL
    ord <- as.character(group.df$variable)
    mat <- mat[ord, ord]
    groups <- TRUE
    ngroups <- nlevels(group.df$group)

    dist <- cumsum(table(group.df$group))
    group_x_min <- c(0.52, dist + 0.5)[1:ngroups]
    group_y_min <- (n + 1) - group_x_min
    group_x_max <- (dist + 0.5)[1:ngroups]
    group_x_max[ngroups] <- n + 0.48
    group_y_max <- (n + 1) - group_x_max
  }

  is_cor_mat <- TRUE
  effect <- "Correlation matrix"
  effect_short <- "cor"
  effect_short2 <- "Cor."
  if (any(diag(mat) != 1)) {
    is_cor_mat <- FALSE
    effect <- "Covariance matrix"
    effect_short <- "cov"
    effect_short2 <- "Cov."
  }

  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("var1", "var2", effect_short)

  if (is_cor_mat) {
    df[[effect_short]][df$var1 == df$var2] <- NA
  }
  levs <- unique(df$var1)
  df$var1 <- factor(as.numeric(as.character(df$var1)), levels = levs)
  df$var2 <- factor(as.numeric(as.character(df$var2)), levels = levs)

  if (order) {
    if (!is_cor_mat) {
      mat <- stats::cov2cor(mat)
    }
    dis_mat <- 1 - mat

    if (ngroups == 1) {
      order2 <- cluster::agnes(x = dis_mat, diss = TRUE, method = "average")$order.lab
    } else if (ngroups > 1) {
      order2 <- c()
      for (i in 1:ngroups) {
        ord <- as.character(group.df$variable[group.df$group == i])
        dis_mat_tmp <- as.matrix(dis_mat[ord, ord])
        if (ncol(dis_mat_tmp) == 1) {
          order2 <- c(order2, ord)
        } else if (ncol(dis_mat_tmp) > 1) {
          order_tmp <- cluster::agnes(x = dis_mat_tmp, diss = TRUE, method = "average")$order.lab
          order2 <- c(order2, order_tmp)
        }
      }
    }
    df$var1 <- factor(df$var1, levels = order2)
    df$var2 <- factor(df$var2, levels = order2)
  }

  var1 <- var2 <- NULL
  if (is_cor_mat) {
    mid_pt <- 0
    max_pt <- 1.1
  } else {
    mid_pt <- mean(df[[effect_short]], na.rm = TRUE)
    max_pt <- max(abs(c(mid_pt - min(df[[effect_short]], na.rm = TRUE), max(df[[effect_short]], na.rm = TRUE) - mid_pt)), na.rm = TRUE) + 1e-8
  }

  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = var1, y = var2)) +
    ggplot2::geom_tile(ggplot2::aes(fill = get(effect_short))) +
    ggplot2::scale_fill_gradient2(
      low = "#195696", mid = "#fcfce1", high = "#A51122", na.value = "transparent",
      midpoint = mid_pt, limits = c(mid_pt - max_pt, mid_pt + max_pt)
    ) +
    ggplot2::scale_x_discrete(expand = c(0.0001, 0.0001)) +
    ggplot2::scale_y_discrete(limits = rev, expand = c(0.0001, 0.0001)) +
    ggplot2::xlab("Variable") +
    ggplot2::ylab("Variable") +
    ggplot2::theme_grey(base_size = 10) +
    ggplot2::ggtitle(label = effect) +
    ggplot2::labs(fill = effect_short2) +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 9),
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 11),
      panel.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(t = 4, r = 0, b = 6, l = 0), size = 12, colour = "gray40")
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = rep(seq(0.5, n - 0.5, 1), each = n), xmax = rep(seq(1.5, n + 0.5, 1), each = n),
        ymin = rep(seq(0.5, n - 0.5, 1), n), ymax = rep(seq(1.5, n + 0.5, 1), n)
      ),
      fill = "transparent", colour = "black", linewidth = 0.05, inherit.aes = FALSE
    ) +
    ggplot2::annotate(
      geom = "rect", xmin = 0.5, ymin = 0.5,
      xmax = n + 0.5, ymax = n + 0.5,
      fill = "transparent", col = "black", lwd = 1.6
    )

  if (labels) {
    p <- p + ggplot2::theme(
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8, r = 0, b = 0, l = 0)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 6, b = 0, l = 0))
    )
  } else {
    p <- p + ggplot2::theme(
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6, r = 0, b = 0, l = 0)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 4, b = 0, l = 0)),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank()
    )
  }

  if (ngroups > 1) {
    for (i in 1:ngroups) {
      p <- p + ggplot2::geom_rect(
        xmin = group_x_min[i],
        ymin = group_y_min[i],
        xmax = group_x_max[i],
        ymax = group_y_max[i],
        color = "black",
        linewidth = 0.5,
        alpha = 0
      )
    }
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
#' # Q-Q plot of the simulated plot errors in the example data frame 'error_df_bivar'
#' # for Trait 1 in Environment 1.
#'
#' error_df <- error_df_bivar[error_df_bivar$env == 1, ]
#'
#' qq <- qq_plot(
#'   df = error_df,
#'   effect = "e.Trait1",
#'   labels = TRUE
#' )
#'
#' # Q-Q plot
#' qq
#'
#' # Extract the data frame with the theoretical and sample quantiles of the
#' # user-defined effects.
#'
#' qq_df <- qq$data
#'
#' @export
qq_plot <- function(df,
                    effect,
                    labels = FALSE) {
  if (is.vector(df)) {
    df <- data.frame(effect = c(df))
    effect <- "effect"
  }
  if (!is.data.frame(df)) {
    stop("'df' must be a data frame")
  }
  colnames(df)[grep("block|col|row", tolower(colnames(df)))] <- tolower(colnames(df))[grep("block|col|row", tolower(colnames(df)))]

  if (!(effect %in% colnames(df))) {
    stop("'df' must contain the effect to be plotted")
  }

  if (!labels) {
    qq_df <- data.frame(effect = df[[effect]])
    p <- ggplot2::ggplot(qq_df, ggplot2::aes(sample = effect)) +
      ggplot2::stat_qq()
    qq_df <- data.frame(
      ggplot2::ggplot_build(p)$data[[1]]["sample"],
      ggplot2::ggplot_build(p)$data[[1]]["theoretical"]
    )
    mid_pt_x <- mean(qq_df$theoretical, na.rm = TRUE)
    max_pt_x <- max(abs(c(mid_pt_x - min(qq_df$theoretical, na.rm = TRUE), max(qq_df$theoretical, na.rm = TRUE) - mid_pt_x)), na.rm = TRUE) + 1e-8

    p <- ggplot2::ggplot(data = qq_df, ggplot2::aes(x = theoretical, y = sample)) +
      ggplot2::stat_qq_line(data = qq_df, ggplot2::aes(sample = sample), colour = "steelblue", linewidth = 0.75, inherit.aes = F) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(y = "Sample quantiles", x = "Theoretical quantiles") +
      ggplot2::ggtitle(label = effect) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(margin = ggplot2::margin(t = 4, r = 0, b = 6, l = 0), size = 12, colour = "gray40"),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6, r = 0, b = 0, l = 0), size = 11),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 4, b = 0, l = 0), size = 11),
        axis.text = ggplot2::element_text(size = 10)
      ) +
      ggplot2::lims(x = c(mid_pt_x - max_pt_x, mid_pt_x + max_pt_x))
    return(p)
  }

  if (labels) {
    if (any(!c("col", "row") %in% colnames(df))) {
      stop("'df' must contain the columns 'col' and 'row' if labels are to be plotted")
    }
    qq_df <- data.frame(
      col = df[["col"]],
      row = df[["row"]],
      effect = df[[effect]]
    )
    qq_df$col <- factor(as.numeric(as.character(qq_df$col)))
    qq_df$row <- factor(as.numeric(as.character(qq_df$row)))
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

    mid_pt_x <- mean(qq_df$theoretical, na.rm = TRUE)
    max_pt_x <- max(abs(c(mid_pt_x - min(qq_df$theoretical, na.rm = TRUE), max(qq_df$theoretical, na.rm = TRUE) - mid_pt_x)), na.rm = TRUE) + 1e-8

    qq_df$cr.label <- factor(paste0(qq_df$col, ":", qq_df$row))
    theoretical <- cr.label <- NULL
    p <- ggplot2::ggplot(data = qq_df, ggplot2::aes(x = theoretical, y = sample, label = cr.label)) +
      ggplot2::stat_qq_line(data = qq_df, ggplot2::aes(sample = sample), colour = "steelblue", linewidth = 0.75, inherit.aes = F) +
      ggplot2::geom_text(size = 4) +
      ggplot2::labs(
        y = "Sample quantiles", x = "Theoretical quantiles",
        subtitle = "Effects indexed as col:row"
      ) +
      ggplot2::ggtitle(label = effect) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(margin = ggplot2::margin(t = 4, r = 0, b = 6, l = 0), size = 12, colour = "gray40"),
        plot.subtitle = ggplot2::element_text(size = 10, colour = "gray40"),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6, r = 0, b = 0, l = 0), size = 11),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 4, b = 0, l = 0), size = 11),
        axis.text = ggplot2::element_text(size = 10)
      ) +
      ggplot2::lims(x = c(mid_pt_x - max_pt_x, mid_pt_x + max_pt_x))
    return(p)
  }
}

#' Sample variogram
#'
#' Creates a sample variogram for a set of effects (e.g., plot errors).
#'
#' @param df A data frame with the columns 'col', 'row', and the effects to be plotted.
#' @param effect The name of the effects to be plotted.
#' @param min.np Minimum number of pairs for which semivariances are displayed (default is 30).
#'
#' @return A sample variogram with x- and y-axes displaying the row and
#'   column displacements, and the z-axis displaying the average semivariances (variogram ordinates)
#'   for the effects.
#'
#' @examples
#' # Sample variogram of plot errors simulated using a separable first order
#' # autoregressive (AR1) process.
#'
#' error_df <- field_trial_error(
#'   ntraits = 1,
#'   nenvs = 1,
#'   spatial.model = "AR1"
#' )
#'
#' variogram <- sample_variogram(
#'   df = error_df,
#'   effect = "e.Trait1"
#' )
#'
#' # Sample variogram
#'
#' variogram
#'
#' # Extract the data frame with the column and row displacements, and the
#' # average semivariances.
#'
#' variogram_df <- variogram$data
#'
#' @export
sample_variogram <- function(df,
                             effect,
                             min.np = 30) {
  if (!is.data.frame(df)) {
    stop("'df' must be a data frame")
  }
  colnames(df)[grep("block|col|row", tolower(colnames(df)))] <- tolower(colnames(df))[grep("block|col|row", tolower(colnames(df)))]

  if (any(!c("col", "row", effect) %in% colnames(df))) {
    stop("'df' must contain the columns 'col' and 'row', and the effect to be plotted")
  }

  variogram_df <- data.frame(
    col = df[["col"]],
    row = df[["row"]],
    effect = df[[effect]]
  )
  variogram_df <- variogram_df[order(variogram_df$col, variogram_df$row), ]

  col_dis <- abs(outer(as.numeric(as.character(variogram_df$col)), as.numeric(as.character(variogram_df$col)), FUN = "-"))
  row_dis <- abs(outer(as.numeric(as.character(variogram_df$row)), as.numeric(as.character(variogram_df$row)), FUN = "-"))
  var_mat <- outer(variogram_df$effect, variogram_df$effect, FUN = "-")^2 / 2
  variogram_df <- data.frame(
    col.dis = col_dis[upper.tri(col_dis, diag = TRUE)],
    row.dis = row_dis[upper.tri(row_dis, diag = TRUE)],
    semivar = var_mat[upper.tri(var_mat, diag = TRUE)]
  )
  variogram_df <- variogram_df[order(variogram_df$col.dis, variogram_df$row.dis), ]

  variogram_df <- data.frame(
    col.dis = rep(unique(variogram_df$col.dis), each = length(unique(variogram_df$row.dis))),
    row.dis = unique(variogram_df$row.dis),
    np = c(with(variogram_df, tapply(semivar, list(row.dis, col.dis), function(x) length(x[!is.na(x)])))),
    semivar = c(with(variogram_df, tapply(semivar, list(row.dis, col.dis), function(x) mean(x, na.rm = T))))
  )

  lattice::lattice.options(
    layout.heights = list(bottom.padding = list(x = -1), top.padding = list(x = -1.5)),
    layout.widths = list(left.padding = list(x = -1.25), right.padding = list(x = -3))
  )
  graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
  p <- lattice::wireframe(semivar ~ row.dis * col.dis,
    data = variogram_df[variogram_df$np >= min.np, ], drape = T, colorkey = F, zoom = 0.97, cuts = 30,
    screen = list(z = 30, x = -60, y = 0), aspect = c(1, 0.66),
    scales = list(distance = c(1.2, 1.2, 0.5), arrows = F, cex = 0.7, col = "black"),
    zlab = list(label = paste("Semivariance"), cex = 0.9, rot = 90, just = c(0.5, -2.25)),
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
#' Creates a theoretical variogram for a separable first order autoregressive (AR1) process.
#'
#' @param ncols A scalar defining the number of columns.
#' @param nrows A scalar defining the number of rows.
#' @param varR A scalar defining the error variance.
#' @param col.cor A scalar defining the column autocorrelation,
#' @param row.cor A scalar defining the row autocorrelation.
#' @param prop.spatial A scalar defining the proportion of spatial trend.
#'
#' @return A theoretical variogram with x- and y-axes displaying the row and column displacements,
#'  and the z-axis displaying the semivariances (variogram ordinates) for a separable autoregressive process.
#'
#' @examples
#' # Theoretical variogram for a field trial with 10 columns and 20 rows, based
#' # on column and row autocorrelations of 0.5 and 0.7, and a proportion of
#' # spatial trend of 0.5. The remaining proportion represents random error.
#'
#' variogram <- theoretical_variogram(
#'   ncols = 10,
#'   nrows = 20,
#'   varR = 1,
#'   col.cor = 0.5,
#'   row.cor = 0.7,
#'   prop.spatial = 0.5
#' )
#'
#' # Theoretical variogram
#'
#' variogram
#'
#' # Extract the data frame with the column and row displacements, and the
#' # theoretical semivariances.
#'
#' variogram_df <- variogram$data
#'
#' @export
theoretical_variogram <- function(ncols = 10,
                                  nrows = 20,
                                  varR = 1,
                                  col.cor = 0.5,
                                  row.cor = 0.7,
                                  prop.spatial = 1) {
  prop_rand <- 1 - prop.spatial
  col_dis <- rep(0:(ncols - 1), each = nrows)
  row_dis <- rep(0:(nrows - 1), times = ncols)
  variogram_df <- data.frame(
    col.dis = col_dis,
    row.dis = row_dis,
    semivar = varR * (prop_rand + prop.spatial * (1 - col.cor^(col_dis) * row.cor^(row_dis)))
  )
  variogram_df$semivar[1] <- 0
  variogram_df$col.dis <- as.numeric(as.character(variogram_df$col.dis))
  variogram_df$row.dis <- as.numeric(as.character(variogram_df$row.dis))

  lattice::lattice.options(
    layout.heights = list(bottom.padding = list(x = -1), top.padding = list(x = -1.5)),
    layout.widths = list(left.padding = list(x = -1.25), right.padding = list(x = -3))
  )
  graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
  p <- lattice::wireframe(semivar ~ row.dis * col.dis,
    data = variogram_df, drape = T, colorkey = F, zoom = 0.97, cuts = 30,
    screen = list(z = 30, x = -60, y = 0), aspect = c(1, 0.66),
    scales = list(distance = c(1.2, 1.2, 0.5), arrows = F, cex = 0.7, col = "black"),
    zlab = list(label = paste("Semivariance"), cex = 0.9, rot = 90, just = c(0.5, -2.25)),
    xlab = list(label = paste("Row displacement"), cex = 0.9, rot = 19, just = c(0.5, -0.75)),
    ylab = list(label = paste("Column displacement"), cex = 0.9, rot = -49, just = c(0.5, -0.75)),
    par.settings = list(axis.line = list(col = "transparent"), clip = list(panel = "off"))
  )
  variogram_df$col.dis <- factor(variogram_df$col.dis)
  variogram_df$row.dis <- factor(variogram_df$row.dis)
  p$data <- variogram_df
  return(p)
}
