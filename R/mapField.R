#' Simulate a correlation matrix between traits or environments
#'
#' Creates a \code{p x p} correlation matrix, which can be used for traits or
#' environments between genotypes or residuals. An upper and lower bound of the
#' sampled correlations can be set.
#'
#' @param p A scalar defining the dimension of the correlation matrix
#' @param minCor A scalar defining the minimum correlation. By default,
#'   minCor = -1.
#' @param maxCor A scalar defining the maximum correlation. By default,
#'   maxCor = 1.
#' @param digits Decimal digits. By default, digits = 2.
#'
#' @return A p x p correlation matrix.
#'
#' @examples
#' # corA <- rand_cor_mat(10, minCor = -0.2, maxCor = 0.8)
#'
#' @export
rand_cor_mat <- function(p, minCor = -1, maxCor = 1, digits = 2) {
  if (p < 1 | p %% 1 != 0) stop("'p' must be an integer > 0")

  if (minCor < -1 | minCor >= 1) stop("'minCor' must be value >= -1 and < 1")
  if (maxCor <= -1 | minCor > 1) stop("'maxCor' must be value > -1 and <= 1")
  if (maxCor < minCor) stop("'maxCor' must not be smaller than 'minCor'")

  nCor <- sum(seq(1, (p - 1)))

  offDg <- round(stats::runif(nCor, min = minCor, max = maxCor), digits = digits)
  corMat <- diag(p)
  corMat[lower.tri(corMat, diag = FALSE)] <- offDg
  corMat <- t(corMat)
  corMat[lower.tri(corMat, diag = FALSE)] <- offDg

  return(corMat)
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
#' nEnvs <- 3 # number of simulated environments.
#' nTraits <- 2 # number of simulated traits.
#'
#' # Field layout
#' nCols <- 10 # total number of columns in each environment.
#' nRows <- c(20, 30, 20) # total number of rows per environment.
#' plotLength <- 5 # plot length of 5 meters.
#' plotWidth <- 2 # plot width of 2 meters.
#' nReps <- c(2, 3, 2) # number of complete replicates per environment.
#'
#' # Residual variances for traits 1 and 2
#' varR <- c(0.4, 15)
#'
#' # Residual correlations between traits 1 and 2, with regards to spatial model
#' corR <- matrix(c(1.0, 0.2, 0.2, 1.0), ncol = 2)
#'
#' plot_df <- field_error(
#'   nEnvs = nEnvs, nTraits = nTraits,
#'   nCols = nCols, nRows = nRows, plotLength = plotLength,
#'   plotWidth = plotWidth, nReps = nReps, repDir = "row",
#'   varR = varR, corR = corR, spatialModel = "bivariate",
#'   propSpatial = 0.6, complexity = 14, effects = FALSE
#' )
#'
#' # Plot the simulated error for trait 2 in environment 2.
#' map_field(plot_df, env = 2, trait = "e.Trait.2")
#'
#' @export
map_field <- function(df, env, trait, borders = TRUE) {
  if (inherits(df, "list")) df <- data.frame(df[[1]])

  trt <- which(colnames(df) == trait)
  df <- subset(df, env == env)
  nRows <- length(unique(df$row))
  nCols <- length(unique(df$col))

  plotMat <- matrix(numeric(), nrow = nRows, ncol = nCols)

  for (i in 1:nrow(df)) {
    r <- df$row[i]
    c <- df$col[i]
    plotMat[r, c] <- df[i, trt]
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

  xLabs <- seq(2, ncol(plotMat), 2)
  xTicks <- (seq(2, ncol(plotMat), 2) - 0.5)

  yLabs <- seq(2, nrow(plotMat), 2)
  yTicks <- (seq(2, nrow(plotMat), 2) - 0.5)

  fields::image.plot(
    x = 0:nCols[1], y = 0:nRows[1],
    z = t(plotMat), zlim = range(plotMat),
    ylim = rev(range(0:nRows[1])),
    col = grDevices::hcl.colors(n = 100000, "RdYlGn"),
    xlab = "Column", ylab = "Row", axes = FALSE
  )

  graphics::box()
  graphics::axis(1, at = xTicks, labels = xLabs)
  graphics::axis(2, at = yTicks, labels = yLabs)

  if (borders == TRUE & length(unique(df$block)) > 1) {
    graphics::grid(nx = nx, ny = ny, lty = 1, col = "#000000", lwd = 5)
    graphics::grid(nx = nx, ny = ny, lty = 1, col = "#FFFFFF", lwd = 3)
  }
}
