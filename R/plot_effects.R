#' Graphics for plot-level effects
#'
#' Graphically displays plot-level effects (e.g., phenotypic values, genetic values, errors) onto
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
#' # Simulation of plot-level errors for two traits in three environments using a bivariate
#' # interpolation model for spatial variation.
#'
#' n_envs <- 3 # Number of simulated environments.
#' n_traits <- 2 # Number of simulated traits.
#'
#' # Field layout
#' n_cols <- 10 # Total number of columns in each environment.
#' n_rows <- c(20, 30, 20) # Total number of rows in each environment.
#' plot_length <- 5 # Plot length set to 5 meters in each environment.
#' plot_width <- 2 # Plot width set to 2 meters in each environment.
#' n_reps <- c(2, 3, 2) # Number of complete replicates (blocks) per environment.
#'
#' # Error variances for traits 1 and 2.
#' var_R <- c(0.4, 15)
#'
#' # Spatial error correlations between traits 1 and 2.
#' S_cor_R <- matrix(c(
#'   1.0, 0.2,
#'   0.2, 1.0
#' ),
#' ncol = 2
#' )
#'
#' error_df <- field_trial_error(
#'   n_envs = n_envs,
#'   n_traits = n_traits,
#'   n_cols = n_cols,
#'   n_rows = n_rows,
#'   plot_length = plot_length,
#'   plot_width = plot_width,
#'   n_reps = n_reps,
#'   rep_dir = "row",
#'   var_R = var_R,
#'   S_cor_R = S_cor_R,
#'   spatial_model = "bivariate",
#'   prop_spatial = 0.6,
#'   complexity = 14,
#'   return_effects = TRUE
#' )
#'
#' # Display the simulated error for trait 2 in environment 2.
#' plot_effects(error_df,
#'   env = 2,
#'   effect = "e.Trait.2"
#' )
#' @export
plot_effects <- function(df,
                         env,
                         effect,
                         blocks = TRUE) {
  if (inherits(df, "list")) df <- data.frame(df[[1]])

  eff <- which(colnames(df) == effect)
  df <- df[df[["env"]] == env, ]
  n_rows <- length(unique(df$row))
  n_cols <- length(unique(df$col))

  plot_mat <- matrix(numeric(), nrow = n_rows, ncol = n_cols)

  for (i in 1:nrow(df)) {
    r <- df$row[i]
    c <- df$col[i]
    plot_mat[r, c] <- df[i, eff]
  }

  if (length(unique(df$block)) > 1) {
    df1 <- df[df[["block"]] == 1, ]
    df2 <- df[df[["block"]] == 2, ]

    if (any(unique(df1$row) == unique(df2$row)) == FALSE) {
      nx <- 0
      ny <- max(df$block)
    } else if (any(unique(df1$col) == unique(df2$col)) == FALSE) {
      nx <- max(df$block)
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

  if (blocks == TRUE & length(unique(df$block)) > 1) {
    graphics::grid(nx = nx, ny = ny, lty = 1, col = "#000000", lwd = 5)
    graphics::grid(nx = nx, ny = ny, lty = 1, col = "#FFFFFF", lwd = 3)
  }
}
