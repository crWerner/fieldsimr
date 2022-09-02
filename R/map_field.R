#' Simulate a correlation matrix between traits or environments
#'
#' Creates a \code{p x p} correlation matrix, which can be used for traits or
#' environments between genotypes or residuals. An upper and lower bound of the
#' sampled correlations can be set.
#'
#' @param p A scalar defining the dimension of the correlation matrix
#' @param min_cor A scalar defining the minimum correlation. By default,
#'   min_cor = -1.
#' @param max_cor A scalar defining the maximum correlation. By default,
#'   max_cor = 1.
#' @param n_digits Number of decimal digits. By default, n_digits = 2.
#'
#' @return A p x p correlation matrix.
#'
#' @examples
#' # cor_A <- rand_cor_mat(10, min_cor = -0.2, max_cor = 0.8)
#'
#' @export
rand_cor_mat <- function(p,
                         min_cor = -1,
                         max_cor = 1,
                         n_digits = 2) {

  if (p < 1 | p %% 1 != 0) stop("'p' must be an integer > 0")

  if (min_cor < -1 | min_cor >= 1) stop("'min_cor' must be value >= -1 and < 1")
  if (max_cor <= -1 | min_cor > 1) stop("'max_cor' must be value > -1 and <= 1")
  if (max_cor < min_cor) stop("'max_cor' must not be smaller than 'min_cor'")

  n_cor <- sum(seq(1, (p - 1)))

  off_dg <- round(stats::runif(n_cor, min = min_cor, max = max_cor), n_digits)
  cor_mat <- diag(p)
  cor_mat[lower.tri(cor_mat, diag = FALSE)] <- off_dg
  cor_mat <- t(cor_mat)
  cor_mat[lower.tri(cor_mat, diag = FALSE)] <- off_dg

  return(cor_mat)
}




#' Map field trial
#'
#' Plots the values of some input trait in the field represented by a colour
#' gradient going from red (low value) to green (high value). The input  trait
#' must have a unique value for each row x column combination within an
#' environment. The function \code{map_field} was created to take data frames
#' generated
#' with \link[FieldSimR]{field_error} as an input, but can work with every data
#' frame that contains at least the columns "env", "col", "row" and the trait to
#' be plotted (e.g., the simulated plot-level residual "e.Trait.1").  If the
#' input data frame contains a column "blocks", the blocks of full replicated
#' will be indicates in the field plot.
#'
#' @param df Data frame containing at least columns "env", "row", "col" and
#'   the trait to be plotted. If df contains a column "block", the blocks of
#'   full replicates are also indicated, if \code{borders = TRUE}. If df is a
#'   list, the first element of the list will be used by default. To use a
#'   different list element, this has to be specified.
#' @param env Environment id of the field to be plotted.
#' @param trait Column identifier of the trait to be plotted.
#' @param borders When true, blocks of full replicates are indicated.
#'
#' @return A field plot of the input trait represented by a colour gradient from
#'   red (low value) to green (high value).
#'
#' @examples
#' # Simulation of residuals for two traits tested in three environments using
#' # bivariate interpolation to model spatial variation.
#'
#' n_envs <- 3 # number of simulated environments.
#' n_traits <- 2 # number of simulated traits.
#'
#' # Field layout
#' n_cols <- 10 # total number of columns in each environment.
#' n_rows <- c(20, 30, 20) # total number of rows per environment.
#' plot_length <- 5 # plot length of 5 meters.
#' plot_width <- 2 # plot width of 2 meters.
#' n_reps <- c(2, 3, 2) # number of complete replicates per environment.
#'
#' # Residual variances for traits 1 and 2
#' var_R <- c(0.4, 15)
#'
#' # Residual cor_Relations between traits 1 and 2, with regards to spatial model
#' cor_R <- matrix(c(1.0, 0.2, 0.2, 1.0), ncol = 2)
#'
#' plot_df <- field_error(
#'   n_envs = n_envs, n_traits = n_traits,
#'   n_cols = n_cols, n_rows = n_rows, plot_length = plot_length,
#'   plot_width = plot_width, n_reps = n_reps, rep_dir = "row",
#'   var_R = var_R, cor_R = cor_R, spatial_model = "bivariate",
#'   prop_spatial = 0.6, complexity = 14, effects = FALSE
#' )
#'
#' # Plot the simulated error for trait 2 in environment 2.
#' map_field(plot_df, env = 2, trait = "e.Trait.2")
#'
#' @export
map_field <- function(df,
                      env,
                      trait,
                      borders = TRUE) {

  if (inherits(df, "list")) df <- data.frame(df[[1]])

  trt <- which(colnames(df) == trait)
  df <- subset(df, env == env)
  n_rows <- length(unique(df$row))
  n_cols <- length(unique(df$col))

  plot_mat <- matrix(numeric(), nrow = n_rows, ncol = n_cols)

  for (i in 1:nrow(df)) {
    r <- df$row[i]
    c <- df$col[i]
    plot_mat[r, c] <- df[i, trt]
  }

  if (length(unique(df$block)) > 1) {
    df1 <- subset(df, df$block == 1)
    df2 <- subset(df, df$block == 2)

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
    col = grDevices::hcl.colors(n = 100000, "RdYlGn"),
    xlab = "Column", ylab = "Row", axes = FALSE
  )

  graphics::box()
  graphics::axis(1, at = x_ticks, labels = x_labs)
  graphics::axis(2, at = y_ticks, labels = y_labs)

  if (borders == TRUE & length(unique(df$block)) > 1) {
    graphics::grid(nx = nx, ny = ny, lty = 1, col = "#000000", lwd = 5)
    graphics::grid(nx = nx, ny = ny, lty = 1, col = "#FFFFFF", lwd = 3)
  }
}
