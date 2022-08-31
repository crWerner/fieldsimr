#' Simulate plot-level error in a multi-environment field trial
#'
#' Creates a data frame of residuals for each plot in a simulated multi-environment
#' field trial testing one or multiple traits. Residuals are based on a spatial
#' component and a random component. The spatial component can be simulated
#' using bivariate interpolation based on akima's \link[akima]{interp} function
#' or based on a separable autoregressive process (AR1:AR1) withinmenvironments.
#' If multiple traits are simulated, a correlated residual betweenmtraits can be
#' simulated based assuming a correlation of the spatial residual between traits,
#' a correlation of the random residual between traits, or ancombination of both.
#' Between traits and environments, a separable covariance structure is assumed.
#'
#' @param nEnvs Number of environments to be simulated (same as used in
#'   \code{compsym_asr_input} or \code{unstr_asr_output}, where applicable).
#' @param nTraits Number of traits to be simulated.
#' @param nCols  A vector containing the number of columns in each environment.
#'   If only one value is provided, all environments will be assigned the same
#'   value.
#' @param nRows A vector containing the number of rows in each environment. If
#'   only one value is provided, all environments will be assigned the same value.
#' @param plotLength A vector of plot lengths for each environment. If only
#'   one value is provided, all environments will be assigned the same value.
#' @param plotWidth A vector of plot widths for each environment. If only one
#'   value is provided, all environments will be assigned the same value.
#' @param nReps A vector containing the number of complete replicates in each
#'   environment. If only one value is provided,all environments will be assigned
#'   the same value.
#' @param repDir Character string specifying the direction of replicate blocks.
#'   One of either "column" (side-by-side, the default) or "row" (above-and-below).
#'   Ignored when \code{nReps = 1}.
#' @param varR A vector of desired residual variances for each trait-by-environment
#'   combination (ordered as environments within traits). If the length of
#'   \code{varR} corresponds to \code{nTraits}, the traits will be assigned the
#'   same residual variance for each environment.
#' @param corR A matrix of residual correlations between more than one traits
#'   with regards to the spatial model only. If not defined and
#'   \code{nTraits > 1}, a diagonal matrix is assigned.
#' @param RcorR A matrix of residual correlations between more than one traits
#'   with regards to the random error model only. If not defined and
#'   \code{nTraits > 1}, a diagonal matrix is assigned.
#' @param spatialModel Character string specifying the model used to simulate
#'   spatial variation. One of either "Bivariate" (bivariate interpolation, the
#'   default) or "AR1:AR1" (two-dimensional autoregressive process of order one).
#' @param propSpatial A vector containing the proportion of residual variance for
#'   the spatial model compared to the total (spatial + random) error variance.
#'   If only one value is provided, all environments will be assigned the same
#'   value. By default, \code{propSpatial = 0.5}.
#' @param complexity A single number indicating the complexity of the bivariate
#'   interpolation model. By default, \code{complexity = 12}. Note that much
#'   lower values might result in convergence problems. For more information on
#'   the complexity parameter check \link[akima]{interp}.
#' @param colCor A vector of column autocorrelations for each environment used in
#'   the AR1:AR1 model. If only one value is provided, all environments will be
#'   assigned the same value.
#' @param rowCor A vector of column autocorrelations for each environment used in
#'   the AR1:AR1 model. If only one value is provided, all environments will be
#'   assigned the same value.
#' @param effects When true, a list is returned with additional list entries for
#'   each trait containing the simulated residuals corresponding to the spatial
#'   model and random error within each environment.
#'
#' @return A data-frame containing environment number, block number, column number,
#'   row number and simulated residuals for each trait. When \code{effects = TRUE},
#'   a list is returned with additional list entries for each trait containing
#'   the simulated residuals corresponding to the spatial model and random error
#'   within each environment.
#'
#' @examples
#' # Simulation of residuals for two traits tested in three environments using
#' # bivariate interpolation to model spatial variation.
#'
#' nEnvs <- 3        # number of simulated environments.
#' nTraits <- 2      # number of simulated traits.
#'
#' # Field layout
#' nCols <- 10              # total number of columns in each environment.
#' nRows <- c(20, 30, 20)   # total number of rows per environment.
#' plotLength <- 5          # plot length of 5 meters.
#' plotWidth <- 2           # plot width of 2 meters.
#' nReps <- c(2, 3, 2)      # number of complete replicates per environment.
#'
#' # Residual variances for traits 1 and 2
#' varR <- c(0.4, 15)
#'
#' # Residual correlations between traits 1 and 2, with regards to spatial model
#' corR <- matrix(c(1.0,  0.2, 0.2, 1.0), ncol = 2)
#'
#' plot_df <- field_error(nEnvs = nEnvs, nTraits = nTraits,
#'                       nCols = nCols, nRows = nRows, plotLength = plotLength,
#'                       plotWidth = plotWidth, nReps = nReps, repDir = "row",
#'                       varR = varR, corR = corR, spatialModel = "bivariate",
#'                       propSpatial = 0.6, complexity = 14, effects = TRUE)
#'
#' @export
field_error <- function(nEnvs,
                        nTraits,
                        nReps,
                        nCols,
                        nRows,
                        plotLength,
                        plotWidth,
                        repDir = "column",
                        varR,
                        corR = NULL,
                        RcorR = NULL,
                        spatialModel = "bivariate",
                        propSpatial = 0.5,
                        complexity = 12,
                        colCor,
                        rowCor,
                        effects = FALSE) {

  if (nEnvs < 1 | nEnvs %% 1 != 0) stop("'nEnvs' must be an integer > 0")
  if (nTraits < 1 | nTraits %% 1 != 0) stop("'nTraits' must be an integer > 0")

  if (min(nCols) < 1 | any(nCols %% 1 != 0)) {
    stop("'nCols' must contain integers > 0")
  }
  if (length(nCols) == 1) nCols <- rep(nCols, nEnvs)
  if (length(nCols) != nEnvs) {
    stop("Length of vector 'nCols' does not match total number of environments")
  }
  if (min(nRows) < 1 | any(nRows %% 1 != 0)) {
    stop("'nRows' must contain integers > 0")
  }
  if (length(nRows) == 1) nRows <- rep(nRows, nEnvs)
  if (length(nRows) != nEnvs) {
    stop("Length of vector 'nRows' does not match total number of environments")
  }

  if (min(plotLength) <= 0) stop("'plotLength' must contain values > 0")
  if (length(plotLength) == 1) plotLength <- rep(plotLength, nEnvs)
  if (length(plotLength) != nEnvs) {
    stop("Length of vector 'plotLength' does not match total number of environments")
  }
  if (min(plotWidth) <= 0) stop("'plotWidth' must be > 0")
  if (length(plotWidth) == 1) plotWidth <- rep(plotWidth, nEnvs)
  if (length(plotWidth) != nEnvs) {
    stop("Length of vector 'plotWidth' does not match total number of environments")
  }

  if (min(nReps) < 1 | any(nReps %% 1 != 0)) {
    stop("'nReps' must contain integers > 0")
  }

  repDir <- tolower(repDir)
  if (repDir == "column"){
    if (any((nCols / nReps) %% 1 != 0)) {
      stop("Number of columns not divisible by number of reps in at least one
         environment. Review your trial design!")
    }
  } else if (repDir == "row"){
    if (any((nRows / nReps) %% 1 != 0)) {
      stop("Number of rows not divisible by number of reps in at least one
           environment. Review your trial design!")
    }
  } else {
    stop("'repDir' must be 'row' or 'column'")
  }

  if (length(nReps) == 1) nReps <- rep(nReps, nEnvs)
  if (length(nReps) != nEnvs) {
    stop("Length of vector 'nReps' does not match total number of environments")
  }

  if (length(varR) == nTraits) varR <- rep(varR, nEnvs)
  if (any(varR <= 0)) {
    stop("'varR' must contain values greater than 0")
  }
  if (length(varR) != (nEnvs * nTraits)) {
    stop("Number of values in argument 'varR' must either match number of traits
         or number of trait x environment combinations")
  }

  if (is.null(corR)) corR <- diag(nTraits)
  if (is.null(RcorR)) RcorR <- diag(nTraits)

  spatialModel <- tolower(spatialModel)
  if (spatialModel != "bivariate" & spatialModel != "ar1:ar1") {
    stop("'spatialModel' must be 'bivariate' or 'AR1:AR1'")
  }

  if (any(propSpatial < 0) | any(propSpatial > 1)) {
    stop("'propSpatial' must contain values between 0 and 1")
  }
  if (length(propSpatial) == 1) propSpatial <- rep(propSpatial, nEnvs)
  if (length(propSpatial) != nEnvs) {
    stop("Length of vector 'propSpatial' does not match total number of environments")

  }

  envs <- rep(1:nEnvs, times = nCols * nRows)
  reps <- c(unlist(mapply(function(x, y, z) rep(1:z, each = c(x * y / z)), x = nCols, y = nRows, z = nReps)))
  cols <- c(unlist(mapply(function(x, y) rep(1:x, each = y), x = nCols, y = nRows)))
  rows <- c(unlist(mapply(function(x, y) rep(1:x, times = y), y = nCols, x = nRows)))

  if (repDir == "row"){
    cols <- c(unlist(mapply(function(x, y) rep(1:x, times = y), x = nCols, y = nRows)))
    rows <- c(unlist(mapply(function(x, y) rep(1:x, each = y), y = nCols, x = nRows)))
  }

  plot.df <- data.frame(env = envs,
                        block = reps,
                        col = cols,
                        row = rows)

  plot.df <- plot.df[order(plot.df$env, plot.df$col, plot.df$row), ]
  rownames(plot.df) <- NULL

  if(spatialModel == "ar1:ar1"){

    if (any(colCor < 0) | any(colCor > 1)) {
      stop("'colCor' must contain values between 0 and 1'")
    }
    if (length(colCor) == 1) colCor <- rep(colCor, nEnvs)
    if (length(colCor) != nEnvs) {
      stop("Length of vector 'colCor' does not match total number of environments")
    }
    if (any(rowCor < 0) | any(rowCor > 1)) stop('rowCor must be between 0 and 1')
    if (length(rowCor) == 1) rowCor <- rep(rowCor, nEnvs)
    if (length(rowCor) != nEnvs) {
      stop("Length of vector 'rowCor' does not match total number of environments")
    }

    power.lst <- lapply(nCols, function(x) abs(outer(1:x, 1:x, "-")))
    colAR1 <- mapply(function(x, y) x^y, x = colCor, y = power.lst, SIMPLIFY = FALSE)

    power.lst <- lapply(nRows, function(x) abs(outer(1:x, 1:x, "-")))
    rowAR1 <- mapply(function(x, y) x^y, x = rowCor, y = power.lst, SIMPLIFY = FALSE)

    corMat.lst <- mapply(function(x, y) kronecker(x, y), x = colAR1, y = rowAR1, SIMPLIFY = FALSE)
    corMat.lst <- mapply(function(x) kronecker(corR, x), x = corMat.lst, SIMPLIFY = FALSE)

    L.lst1 <- lapply(corMat.lst, function(x) chol(x))
    plotError.lst1 <- mapply(function(x, y) matrix(c(stats::rnorm(x) %*% y), ncol = nTraits),
                             x = nCols * nRows * nTraits, y = L.lst1, SIMPLIFY = FALSE)

  }

  if(spatialModel == "bivariate"){

    if (complexity <= 0) stop("'complexity' must be an integer > 0")

    nPlots <- nCols * nRows
    cols.lst <- with(plot.df, tapply(col, env, function(x) x))
    colcentres.lst <- mapply(function(x, y, z) rep(y, z)* (x - 0.5), x = cols.lst,
                             y = plotLength, z = nPlots, SIMPLIFY = FALSE)

    colcentres <- unlist(colcentres.lst)
    colcentres.lst <- lapply(colcentres.lst, function(x) unique(x))

    rows.lst <- with(plot.df, tapply(row, env, function(x) x))
    rowcentres.lst <- mapply(function(x, y, z) rep(y, z)* (x - 0.5), x = rows.lst,
                             y = plotWidth, z = nPlots, SIMPLIFY = FALSE)
    rowcentres <- unlist(rowcentres.lst)
    rowcentres.lst <- lapply(rowcentres.lst, function(x) unique(x))

    colGap <- plotLength / 4
    rowGap <- plotWidth / 4

    xInterp.list <- mapply(function(x, y, z) c(0 - z, x * y + z, 0 - z, x * y + z, sample(stats::runif(n = complexity, min = 0, max = (x * y)))),
                           x = nCols, y = plotLength, z = colGap, SIMPLIFY = FALSE)
    yInterp.list <- mapply(function(x, y, z) c(0 - z, 0 - z, x * y + z, x * y + z, sample(stats::runif(n = complexity, min = 0, max = (x * y)))),
                           x = nRows, y = plotWidth, z = rowGap, SIMPLIFY = FALSE)
    zInterp.list <- lapply(nCols, function(x) scale(matrix(stats::rnorm((4 + complexity) * nTraits), ncol = nTraits)) %*% chol(corR))


    for(i in 1:nTraits) {
      tmp <- mapply(function(x, y, z, xo, yo)
        c(t(akima::interp(x = x, y = y, z = z[, i], xo = xo, yo = yo, linear = F, extrap = T, duplicate = "mean")$z)),
        x = xInterp.list, y = yInterp.list, z = zInterp.list, xo = colcentres.lst, yo = rowcentres.lst, SIMPLIFY = FALSE)
      if(i == 1) {plotError.lst1 <- tmp}
      if(i > 1) {plotError.lst1 <- Map("cbind", plotError.lst1, tmp)}
    }
  }

  plotError.lst2 <- mapply(function(x) matrix(c(stats::rnorm(x)), ncol = nTraits) %*% chol(RcorR),
                           x = nCols * nRows * nTraits, SIMPLIFY = FALSE)

  varR <- as.data.frame(t(matrix(varR, ncol = nTraits, byrow = TRUE)))

  if(nTraits == 1) {
    varR <- lapply(X = varR, FUN = as.matrix)
  } else {
    varR <- lapply(X = varR, FUN = c)
  }

  e.Spat <- mapply(function(w, x) (scale(w) * sqrt(x)), w = plotError.lst1, x = propSpatial, SIMPLIFY = F)
  e.Rand <- mapply(function(x, y) (scale(y) * sqrt(1 - x)), x = propSpatial, y = plotError.lst2, SIMPLIFY = F)
  e.scale <- mapply(function(x, y) sqrt(diag(1 / diag(as.matrix(stats::var(x + y))))), x = e.Spat, y = e.Rand, SIMPLIFY = F)
  e.Spat <- mapply(function(x, y, z) x %*% y %*% diag(sqrt(z)), x = e.Spat, y = e.scale, z = varR, SIMPLIFY = F)
  e.Rand <- mapply(function(x, y, z) x %*% y %*% diag(sqrt(z)), x = e.Rand, y = e.scale, z = varR, SIMPLIFY = F)
  plotError.lst <- mapply(function(x, y) x + y, x = e.Spat, y = e.Rand, SIMPLIFY = F)
  plotError <- do.call(what = "rbind", plotError.lst)
  colnames(plotError) <- paste0("e.Trait.", 1:nTraits)
  plot.df <- cbind(plot.df, plotError)

  if (effects) {
    e.Spat <- do.call("rbind", e.Spat)
    e.Rand <- do.call("rbind", e.Rand)
    e.all <- lapply(seq_len(ncol(e.Spat)), function(i) cbind(e.Spat[,i], e.Rand[,i]))
    resids <- lapply(e.all, function(x) data.frame(plot.df[,1:4],
                                                   e.Spatial = x[,1],
                                                   e.Random = x[,2]))

    listNames <- c("plot_df", paste0("Trait.", 1:nTraits))
    plot.df <- list(plot.df)
    plot.df <- c(plot.df, resids)
    names(plot.df) <- listNames
  }

  return(plot.df)
}

