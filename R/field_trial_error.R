#' Simulate plot errors in plant breeding field trials
#'
#' Creates a data frame of simulated plot errors in multi-environment field trials for one or
#' more traits. The plot errors capture spatial trend, random error (noise), and extraneous variation.
#' Spatial trend is simulated using bivariate interpolation or a separable first-order
#' autoregressive (AR1) process. Random error is simulated using an independent process. Extraneous
#' variation is simulated using random or zig-zag ordering between neighbouring columns and/or rows.
#' The three error components are combined at a user-defined ratio. \cr
#' Correlated plot errors can be simulated between traits by setting different correlation structures
#' for each error component. A separable structure is assumed between traits and plots within environments,
#' but different error variances can be specified for each environment-within-trait combination.
#'
#' @param ntraits Number of traits to be simulated.
#' @param nenvs Number of environments to be simulated.
#' @param nblocks A vector defining the number of blocks in each environment.
#'   If only one value is specified, all environments will be assigned the same number.
#' @param block.dir A vector defining the block direction in each environment. Use
#'   'col' for side-by-side (default), 'row' for above-and-below, or
#'   NA if only one block is simulated. If only one value is specified, all environments
#'   will be assigned the same block direction.
#' @param ncols A vector defining the number of columns in each environment. If only one
#'   value is specified, all environments will be assigned the same number.
#' @param nrows A vector defining the number of rows in each environment. If only one
#'   value is specified, all environments will be assigned the same number.
#' @param varR A vector of error variances for each environment-within-trait combination. If only
#'   one value is specified, all combinations will be assigned the same error variance.
#' @param ScorR A matrix of spatial error correlations between traits. If not specified and
#'   spatial trend is simulated, a diagonal matrix is constructed.
#' @param RcorR A matrix of random error correlations between traits. If not specified and
#'   random error is simulated, a diagonal matrix is constructed.
#' @param EcorR A matrix of extraneous error correlations between traits. If not specified and
#'   extraneous variation is simulated, a diagonal matrix is constructed. \cr
#'   \strong{Note:} the same correlation between traits is used for the column and row errors.
#'   Currently only implemented when \code{ext.ord = "random"}.
#' @param spatial.model A character string defining the model used to simulate spatial trend.
#'   Use 'Bivariate' for bivariate interpolation (default) or 'AR1' for a separable first-order
#'   autoregressive process. Bivariate interpolation is implemented with the \code{interp}
#'   function of the R package \href{https://CRAN.R-project.org/package=interp}{`interp`}.
#' @param complexity A vector defining the complexity of the simulated spatial trend in each
#'   environment when \code{spatial.model = "Bivariate"}. If only one value is specified,
#'   all environments will be assigned the same complexity.
#'   If not specified and \code{spatial.model = "Bivariate"}, the complexity is set to half the maximum
#'   number of columns and rows in each environment.
#' @param plot.length A vector of plot lengths for each environment (column direction). If only one value is specified,
#'   all environments will be assigned the same plot length. Only required when \code{spatial.model = "Bivariate"}.
#' @param plot.width A vector of plot widths for each environment (row direction). If only one value is specified,
#'   all environments will be assigned the same plot width. Only required when \code{spatial.model = "Bivariate"}.
#' @param col.cor A vector of column autocorrelations for each environment.
#'   If only one value is specified, all environments will be assigned the same column autocorrelation.
#'   Only required when \code{spatial.model = "AR1"}.
#' @param row.cor A vector of row autocorrelations for each environment.
#'   If only one value is specified, all environments will be assigned the same row
#'   autocorrelation. Only required when \code{spatial.model = "AR1"}.
#' @param prop.spatial A vector defining the proportion of spatial trend for each environment-within-trait
#'   combination. If only one value is specified, all combinations will be assigned the proportion.
#' @param ext.ord A character string defining the method used to simulate extraneous variation.
#'   Use 'random' (default) for random variation between neighbouring columns and/or rows or
#'   'zig-zag' for alternating positive and negative values.
#' @param ext.dir A vector defining the direction of extraneous variation for each environment.
#'   Use 'row' (default) for row variation, 'col' for column variation, 'both' for variation in both directions,
#'   or NA if no extraneous variation is simulated. When \code{ext.dir = "both"}, half the variance is assigned to the columns and
#'   half is assigned to the rows. If only one value is specified, all environments will be
#'   assigned the same direction.
#' @param prop.ext A vector defining the proportion of extraneous variation for each environment-within-trait
#'   combination. If only one value is specified, all combinations will be assigned the same proportion.
#' @param return.effects When \code{TRUE} (default is \code{FALSE}), a list is returned with additional entries
#'   containing the spatial, random, and extraneous error terms for each trait.
#'
#' @return A data frame with columns 'env', 'block', 'col', and 'row', followed by the
#'   simulated plot errors for each trait. When \code{return.effects = TRUE}, a list is returned with additional entries
#'   containing the spatial, random, and extraneous error terms for each trait.
#'
#' @examples
#' # Simulate plot errors for two traits in two environments using an AR1 model
#' # for spatial variation.
#'
#' # Error variances for the four environment-within-trait combinations.
#' varR <- c(0.2, 0.4, 10, 15) # Trait 1 x 2 environments, Trait 2 x 2 environments
#'
#' # Spatial error correlations between the two simulated traits.
#' ScorR <- matrix(c(
#'   1.0, 0.2,
#'   0.2, 1.0
#' ), ncol = 2)
#'
#' error_ls <- field_trial_error(
#'   ntraits = 2,
#'   nenvs = 2,
#'   nblocks = 2,
#'   block.dir = "row",
#'   ncols = 10,
#'   nrows = 20,
#'   varR = varR,
#'   ScorR = ScorR,
#'   spatial.model = "AR1",
#'   col.cor = 0.5,
#'   row.cor = 0.7,
#'   prop.spatial = 0.4,
#'   ext.ord = "zig-zag",
#'   ext.dir = "row",
#'   prop.ext = 0.2,
#'   return.effects = TRUE
#' )
#'
#' @export
field_trial_error <- function(ntraits = 1,
                              nenvs = 1,
                              nblocks = 2,
                              block.dir = "col",
                              ncols = 10,
                              nrows = 20,
                              varR = 1,
                              ScorR = NULL,
                              RcorR = NULL,
                              EcorR = NULL,
                              spatial.model = "Bivariate",
                              complexity = NULL,
                              plot.length = 8,
                              plot.width = 2,
                              col.cor = 0.5,
                              row.cor = 0.7,
                              prop.spatial = 0.5,
                              ext.ord = "random",
                              ext.dir = "row",
                              prop.ext = 0,
                              return.effects = FALSE) {
  if (ntraits < 1 | ntraits %% 1 != 0) stop("'ntraits' must be a positive integer")
  if (nenvs < 1 | nenvs %% 1 != 0) stop("'nenvs' must be a positive integer")

  if (min(ncols) < 1 | any(ncols %% 1 != 0)) {
    stop("'ncols' must contain postive integers")
  }
  if (length(ncols) == 1) ncols <- rep(ncols, nenvs)
  if (length(ncols) != nenvs) {
    stop("Length of 'ncols' must match number of environments")
  }
  if (min(nrows) < 1 | any(nrows %% 1 != 0)) {
    stop("'nrows' must contain postive integers")
  }
  if (length(nrows) == 1) nrows <- rep(nrows, nenvs)
  if (length(nrows) != nenvs) {
    stop("Length of 'nrows' must match number of environments")
  }

  if (min(plot.length) <= 0) stop("'plot.length' must contain positive values")
  if (length(plot.length) == 1) plot.length <- rep(plot.length, nenvs)
  if (length(plot.length) != nenvs) {
    stop("Length of 'plot.length' must match number of environments")
  }
  if (min(plot.width) <= 0) stop("'plot.width' must be > 0")
  if (length(plot.width) == 1) plot.width <- rep(plot.width, nenvs)
  if (length(plot.width) != nenvs) {
    stop("Length of 'plot.width' must match number of environments")
  }

  if (min(nblocks) < 1 | any(nblocks %% 1 != 0)) {
    stop("'nblocks' must contain positive integers")
  }

  block.dir[is.na(block.dir)] <- block.dir[tolower(block.dir) == "column"] <- "col"
  block.dir <- tolower(block.dir)
  if (length(block.dir) == 1) block.dir <- rep(block.dir, nenvs)
  if (any(block.dir == "col")) {
    if (any(((ncols / nblocks) %% 1 != 0)[block.dir == "col"])) {
      stop("Number of columns must be divisible by number of blocks. Review your trial design!")
    }
  }
  if (any(block.dir == "row")) {
    if (any(((nrows / nblocks) %% 1 != 0)[block.dir == "row"])) {
      stop("Number of rows must be divisible by number of blocks. Review your trial design!")
    }
  }

  if (length(nblocks) == 1) nblocks <- rep(nblocks, nenvs)
  if (length(nblocks) != nenvs) {
    stop("Length of 'nblocks' must match number of environments")
  }

  if (length(varR) == 1) varR <- rep(varR, ntraits * nenvs)
  if (any(varR < 0)) {
    stop("'varR' must contain values greater than 0")
  }
  if (length(varR) != (nenvs * ntraits)) {
    stop("Length of 'varR' must match number of environment-within-trait combinations")
  }
  if (any(varR == 0)) {
    warning("The error variance is zero for some environments")
  }

  if (is.null(ScorR)) ScorR <- diag(ntraits)
  if (is.null(RcorR)) RcorR <- diag(ntraits)
  if (is.null(EcorR)) EcorR <- diag(ntraits)

  if (any(eigen(ScorR)$values < 1e-8)) {
    stop("'ScorR' must be positive definite")
  }
  if (any(eigen(RcorR)$values < 1e-8)) {
    stop("'RcorR' must be positive definite")
  }
  if (any(eigen(EcorR)$values < 1e-8)) {
    stop("'EcorR' must be positive definite")
  }

  spatial.model <- tolower(spatial.model)
  spatial.model[spatial.model == "ar1:ar1"] <- "ar1"
  if (spatial.model != "bivariate" & spatial.model != "ar1") {
    stop("'spatial.model' must be either 'Bivariate' or 'AR1'")
  }

  if (any(prop.spatial < 0) | any(prop.spatial > 1)) {
    stop("'prop.spatial' must contain values between 0 and 1")
  }
  if (length(prop.spatial) == 1) prop.spatial <- rep(prop.spatial, ntraits * nenvs)
  if (length(prop.spatial) != ntraits * nenvs) {
    stop("Length of 'prop.spatial' must match number of environment-within-trait combinations")
  }

  if (any(prop.ext < 0) | any(prop.ext > 1)) {
    stop("'prop.ext' must contain values between 0 and 1")
  }
  if (length(prop.ext) == 1) prop.ext <- rep(prop.ext, ntraits * nenvs)
  if (length(prop.ext) != ntraits * nenvs) {
    stop("Length of 'prop.ext' must match number of environment-within-trait combinations")
  }

  block.dir[is.na(block.dir)] <- "neither"
  if (length(ext.dir) == 1) ext.dir <- rep(ext.dir, nenvs)
  if (length(ext.dir) != nenvs) {
    stop("Length of 'ext.dir' must match number of environments")
  }
  if (any(prop.ext == 0)) {
    ext.dir[prop.ext == 0] <- "neither"
  }

  if (any(prop.spatial + prop.ext > 1)) {
    stop("The sum of 'prop.spatial' and 'prop.ext' must be between 0 and 1")
  }
  prop_rand <- 1 - prop.spatial - prop.ext

  ext.dir <- tolower(ext.dir)
  if (any(prop.ext > 0) & !all(ext.dir %in% c("col", "row", "both", "neither"))) {
    stop("'ext.dir' must be one of 'col', 'row', or 'both'")
  }

  if (any(ext.dir %in% c("col", "both") & ncols <= 3 & prop.ext > 0)) {
    stop("'ncols' must be greater than 3 when simulating column extraneous variation")
  }

  if (any(ext.dir %in% c("row", "both") & ncols <= 3 & prop.ext > 0)) {
    stop("'nrows' must be greater than 3 when simulating row extraneous variation")
  }

  prop.ext_col <- prop.ext_row <- prop.ext
  prop.ext_col[ext.dir == "row"] <- prop.ext_row[ext.dir == "col"] <- 0
  prop.ext_col[ext.dir == "both"] <- prop.ext_row[ext.dir == "both"] <- prop.ext[ext.dir == "both"] / 2

  ext.ord <- tolower(ext.ord)
  if (!all(ext.ord %in% c("random", "zig-zag"))) {
    stop("'ext.ord' must be one of 'random' or 'zig-zag'")
  }

  envs <- rep(1:nenvs, times = ncols * nrows)
  blocks <- c(unlist(mapply(function(x, y, z) rep(1:z, each = c(x * y / z)), x = ncols, y = nrows, z = nblocks)))

  cols1 <- c(unlist(mapply(function(x, y) rep(1:x, each = y), x = ncols, y = nrows)))
  rows1 <- c(unlist(mapply(function(x, y) rep(1:x, times = y), y = ncols, x = nrows)))
  cols2 <- c(unlist(mapply(function(x, y) rep(1:x, times = y), x = ncols, y = nrows)))
  rows2 <- c(unlist(mapply(function(x, y) rep(1:x, each = y), y = ncols, x = nrows)))

  n_plots <- ncols * nrows
  n_plots1 <- rep(as.numeric(factor(block.dir, levels = c("row", "col"))) - 1, times = n_plots)
  n_plots2 <- 1 - n_plots1
  cols <- n_plots1 * cols1 + n_plots2 * cols2
  rows <- n_plots1 * rows1 + n_plots2 * rows2

  error_df <- data.frame(
    env = factor(envs),
    block = factor(blocks),
    col = factor(cols),
    row = factor(rows)
  )
  error_df <- error_df[order(error_df$env, error_df$col, error_df$row), ]
  rownames(error_df) <- NULL

  plot_error_lst1 <- lapply(n_plots, function(x) matrix(0, nrow = x, ncol = ntraits))
  if (any(prop.spatial > 0)) {
    if (spatial.model == "ar1") {
      if (length(col.cor) == 1) col.cor <- rep(col.cor, nenvs)
      if (length(col.cor) != nenvs) {
        stop("Length of 'col.cor' must match number of environments")
      }
      if (any(col.cor < -1) | any(col.cor > 1)) {
        stop("'col.cor' must contain values between -1 and 1'")
      }
      if (any(abs(col.cor) == 1)) {
        col.cor[abs(col.cor) == 1] <- sign(col.cor[abs(col.cor) == 1]) * (1 - 1e-7)
      }

      if (length(row.cor) == 1) row.cor <- rep(row.cor, nenvs)
      if (length(row.cor) != nenvs) {
        stop("Length of 'row.cor' must match number of environments")
      }
      if (any(row.cor < -1) | any(row.cor > 1)) {
        stop("'row.cor' must contain values between -1 and 1")
      }
      if (any(abs(row.cor) == 1)) {
        row.cor[abs(row.cor) == 1] <- sign(row.cor[abs(row.cor) == 1]) * (1 - 1e-7)
      }

      power_lst <- lapply(ncols, function(x) abs(outer(1:x, 1:x, "-")))
      col_ar1 <- mapply(function(x, y) x^y, x = col.cor, y = power_lst, SIMPLIFY = FALSE)

      power_lst <- lapply(nrows, function(x) abs(outer(1:x, 1:x, "-")))
      row_ar1 <- mapply(function(x, y) x^y, x = row.cor, y = power_lst, SIMPLIFY = FALSE)

      cor_mat_lst <- mapply(function(x, y) kronecker(t(chol(x)), t(chol(y))), x = col_ar1, y = row_ar1, SIMPLIFY = FALSE)
      x <- T
      y <- 0
      while (x | y == 100) {
        y <- y + 1
        plot_error_lst1 <- mapply(function(x, y) y %*% scale(matrix(c(stats::rnorm(x)), ncol = ntraits)),
          x = ncols * nrows * ntraits, y = cor_mat_lst, SIMPLIFY = FALSE
        )

        x <- any(unlist(lapply(plot_error_lst1, function(x) eigen(stats::var(x))$values < 1e-7)))
        if (y == 100) {
          stop("Appropriate AR1 spatial error not obtained in 100 iterations")
        }
      }
      if (ntraits > 1) {
        plot_error_lst1 <- mapply(function(x) scale(x %*% solve(chol(stats::var(x))) %*% chol(ScorR)),
          x = plot_error_lst1, SIMPLIFY = FALSE
        )
      }
      if (ntraits == 1) {
        plot_error_lst1 <- mapply(function(x) scale(x),
          x = plot_error_lst1, SIMPLIFY = FALSE
        )
      }
    }

    if (spatial.model == "bivariate") {
      if (is.null(complexity)) complexity <- apply(cbind(ncols, nrows), 1, function(x) ceiling(max(x) / 2))
      if (length(complexity) == 1) complexity <- rep(complexity, nenvs)
      if (length(complexity) != nenvs) {
        stop("Length of 'complexity' must match number of environments")
      }
      if (any(complexity < 0)) stop("'complexity' must cotain non-negative integers")

      cols_lst <- with(error_df, tapply(col, env, function(x) c(1:max(as.numeric(trimws(x))), max(as.numeric(trimws(x))) + 1)))
      rows_lst <- with(error_df, tapply(row, env, function(x) c(1:max(as.numeric(trimws(x))), max(as.numeric(trimws(x))) + 1)))
      max_lst <- mapply(function(x, y) 1:max(x, y), x = cols_lst, y = rows_lst, SIMPLIFY = FALSE)

      col_centres_lst <- mapply(function(x, y) round(y * (x - 0.5), 8),
        x = max_lst, y = plot.length, SIMPLIFY = FALSE
      )

      row_centres_lst <- mapply(function(x, y) round(y * (x - 0.5), 8),
        x = max_lst, y = plot.width, SIMPLIFY = FALSE
      )

      gap <- min(plot.length, plot.width) / 2

      plot_error_lst1 <- NA
      y <- 0
      while (sum(is.na(unlist(plot_error_lst1))) > 0 | y == 100) {
        y <- y + 1
        xInterp_list <- mapply(function(w, x, y, z) c(0 - z, x * y + z, 0 - z, x * y + z, sample(stats::runif(n = w + 2, min = 0, max = (x * y)), size = w)),
          w = complexity, x = ncols, y = plot.length, z = gap, SIMPLIFY = FALSE
        )
        yInterp_list <- mapply(function(w, x, y, z) c(0 - z, 0 - z, x * y + z, x * y + z, sample(stats::runif(n = w + 2, min = 0, max = (x * y)), size = w)),
          w = complexity, x = nrows, y = plot.width, z = gap, SIMPLIFY = FALSE
        )
        x <- T
        w <- 0
        while (x | w == 100) {
          w <- w + 1
          zInterp_list <- mapply(function(w) scale(matrix(stats::rnorm((4 + w) * ntraits), ncol = ntraits)),
            w = complexity, SIMPLIFY = FALSE
          )

          x <- any(unlist(lapply(zInterp_list, function(x) eigen(stats::var(x))$values < 1e-7)))
          if (y == 100) {
            stop("Appropriate bivariate spatial error not obtained in 100 iterations")
          }
        }
        if (ntraits > 1) {
          zInterp_list <- mapply(function(x) scale(x %*% solve(chol(stats::var(x))) %*% chol(ScorR)),
            x = zInterp_list, SIMPLIFY = FALSE
          )
        }

        if (ntraits == 1) {
          zInterp_list <- mapply(function(x) scale(x),
            x = zInterp_list, SIMPLIFY = FALSE
          )
        }

        for (i in 1:ntraits) {
          tmp <- mapply(
            function(v, w, x, y, z, xo, yo) {
              c(fill_matrix(t(interp::interp(x = x, y = y, z = z[, i], xo = xo, yo = yo, linear = F, extrap = T, duplicate = "mean")$z)[1:v, 1:w]))
            },
            v = nrows, w = ncols, x = xInterp_list, y = yInterp_list, z = zInterp_list, xo = col_centres_lst, yo = row_centres_lst, SIMPLIFY = FALSE
          )
          if (i == 1) {
            plot_error_lst1 <- tmp
          }
          if (i > 1) {
            plot_error_lst1 <- Map("cbind", plot_error_lst1, tmp)
          }
        }
        if (y == 100) {
          stop("Appropriate bivariate spatial error not obtained in 100 iterations. Consider reseting 'complexity'")
        }
      }
    }
  }

  plot_error_lst2 <- lapply(n_plots, function(x) matrix(0, nrow = x, ncol = ntraits))
  if (any(prop_rand > 0)) {
    x <- T
    y <- 0
    while (x | y == 100) {
      y <- y + 1
      plot_error_lst2 <- mapply(function(x) scale(matrix(c(stats::rnorm(x)), ncol = ntraits)),
        x = ncols * nrows * ntraits, SIMPLIFY = FALSE
      )
      x <- any(unlist(lapply(plot_error_lst2, function(x) eigen(stats::var(x))$values < 1e-7)))
      if (y == 100) {
        stop("Appropriate random error not obtained in 100 iterations")
      }
    }

    if (ntraits > 1) {
      plot_error_lst2 <- mapply(function(x) scale(x %*% solve(chol(stats::var(x))) %*% chol(RcorR)),
        x = plot_error_lst2, SIMPLIFY = FALSE
      )
    }
    if (ntraits == 1) {
      plot_error_lst2 <- mapply(function(x) scale(x),
        x = plot_error_lst2, SIMPLIFY = FALSE
      )
    }
  }

  n_plots <- mapply(function(x, y) x * y, x = ncols, y = nrows, SIMPLIFY = FALSE)
  plot_error_lst3c <- lapply(n_plots, function(x) matrix(0, nrow = x, ncol = ntraits))
  plot_error_lst3r <- lapply(n_plots, function(x) matrix(0, nrow = x, ncol = ntraits))

  if (any(prop.ext_col > 0)) {
    plot_error_lst3c <- mapply(function(x, y) scale(matrix(c(stats::rnorm(x)), ncol = ntraits)),
      x = ncols * ntraits, SIMPLIFY = FALSE
    )

    x <- any(unlist(lapply(plot_error_lst3c, function(x) eigen(stats::var(x))$values < 1e-7)))
    y <- 0
    while (x & all(ncols > ntraits) | y == 100 & all(ncols > ntraits)) {
      plot_error_lst3c <- mapply(function(x, y) scale(matrix(c(stats::rnorm(x)), ncol = ntraits)),
        x = ncols * ntraits, SIMPLIFY = FALSE
      )
      x <- any(unlist(lapply(plot_error_lst3c, function(x) eigen(stats::var(x))$values < 1e-7)))
      y <- y + 1
      if (y == 100) {
        stop("Appropriate column extraneous error not obtained in 100 iterations")
      }
    }

    if (ext.ord == "zig-zag") {
      x <- any(rep(ceiling(ncols / 2), each = ntraits) != unlist(lapply(plot_error_lst3c, function(x) colSums(x <= 0))))
      y <- 0
      while (x) {
        tmp1 <- mapply(function(x, y) scale(matrix(c(stats::rnorm(x)), ncol = ntraits)),
          x = ncols * ntraits, SIMPLIFY = FALSE
        )
        tmp2 <- mapply(function(x, y) x - y, x = lapply(tmp1, function(x) colSums(x <= 0)), y = lapply(tmp1, function(x) rep(ceiling(nrow(x) / 2), ntraits)), SIMPLIFY = FALSE)

        for (i in 1:nenvs) {
          plot_error_lst3c[[i]][, tmp2[[i]] == 0] <- tmp1[[i]][, tmp2[[i]] == 0]
        }
        x <- any(unlist(mapply(function(x, y) x - y, x = lapply(plot_error_lst3c, function(x) colSums(x <= 0)), y = lapply(plot_error_lst3c, function(x) rep(ceiling(nrow(x) / 2), ntraits)), SIMPLIFY = FALSE)) != 0)
        y <- y + 1
        if (y == 1000) {
          stop("Appropriate column extraneous error not obtained in 1000 iterations. Try random extraneous ordering")
        }
      }
      for (i in 1:ntraits) { # i <- 1
        tmp <- suppressWarnings(mapply(function(x, y) c(t(cbind(sample(sort(x[, i])[1:ceiling(nrow(x) / 2)]), sample(rev(sort(x[, i]))[1:floor(nrow(x) / 2)]))))[1:y], x = plot_error_lst3c, y = ncols, SIMPLIFY = FALSE))
        if (i == 1) {
          plot_error_lst3c1 <- tmp
        }
        if (i > 1) {
          plot_error_lst3c1 <- Map("cbind", plot_error_lst3c1, tmp)
        }
      }
      plot_error_lst3c <- plot_error_lst3c1
    }

    if (ntraits > 1 & all(ncols > ntraits) & ext.ord == "random") {
      plot_error_lst3c <- mapply(function(x) scale(x %*% solve(chol(stats::var(x))) %*% chol(EcorR)),
        x = plot_error_lst3c, SIMPLIFY = FALSE
      )
    }
    if (ntraits == 1 | any(ncols <= ntraits) | ext.ord == "zig-zag") {
      plot_error_lst3c <- mapply(function(x) scale(x),
        x = plot_error_lst3c, SIMPLIFY = FALSE
      )
    }
    Zc <- lapply(seq_len(nenvs), function(i) stats::model.matrix(~ col - 1, droplevels(error_df[error_df$env == i, ])))
    plot_error_lst3c <- mapply(function(w, x) scale(w %*% x), w = Zc, x = plot_error_lst3c, SIMPLIFY = FALSE)
  }

  if (any(prop.ext_row > 0)) {
    plot_error_lst3r <- mapply(function(x, y) scale(matrix(c(stats::rnorm(x)), ncol = ntraits)),
      x = nrows * ntraits, SIMPLIFY = FALSE
    )

    x <- any(unlist(lapply(plot_error_lst3r, function(x) eigen(stats::var(x))$values < 1e-7)))
    y <- 0
    while (x & all(nrows > ntraits) | y == 100 & all(nrows > ntraits)) {
      plot_error_lst3r <- mapply(function(x, y) scale(matrix(c(stats::rnorm(x)), ncol = ntraits)),
        x = nrows * ntraits, SIMPLIFY = FALSE
      )
      x <- any(unlist(lapply(plot_error_lst3r, function(x) eigen(stats::var(x))$values < 1e-7)))
      y <- y + 1
      if (y == 100) {
        stop("Appropriate row extraneous error not obtained in 100 iterations")
      }
    }

    if (ext.ord == "zig-zag") {
      x <- any(rep(ceiling(nrows / 2), each = ntraits) != unlist(lapply(plot_error_lst3r, function(x) colSums(x <= 0))))
      y <- 0
      while (x) {
        tmp1 <- mapply(function(x, y) scale(matrix(c(stats::rnorm(x)), ncol = ntraits)),
          x = nrows * ntraits, SIMPLIFY = FALSE
        )
        tmp2 <- mapply(function(x, y) x - y, x = lapply(tmp1, function(x) colSums(x <= 0)), y = lapply(tmp1, function(x) rep(ceiling(nrow(x) / 2), ntraits)), SIMPLIFY = FALSE)

        for (i in 1:nenvs) {
          plot_error_lst3r[[i]][, tmp2[[i]] == 0] <- tmp1[[i]][, tmp2[[i]] == 0]
        }
        x <- any(unlist(mapply(function(x, y) x - y, x = lapply(plot_error_lst3r, function(x) colSums(x <= 0)), y = lapply(plot_error_lst3r, function(x) rep(ceiling(nrow(x) / 2), ntraits)), SIMPLIFY = FALSE)) != 0)
        y <- y + 1
        if (y == 1000) {
          stop("Appropriate row extraneous error not obtained in 1000 iterations. Try random extraneous ordering")
        }
      }
      for (i in 1:ntraits) {
        tmp <- suppressWarnings(mapply(function(x, y) c(t(cbind(sample(sort(x[, i])[1:ceiling(nrow(x) / 2)]), sample(rev(sort(x[, i]))[1:floor(nrow(x) / 2)]))))[1:y], x = plot_error_lst3r, y = nrows, SIMPLIFY = FALSE))
        if (i == 1) {
          plot_error_lst3r1 <- tmp
        }
        if (i > 1) {
          plot_error_lst3r1 <- Map("cbind", plot_error_lst3r1, tmp)
        }
      }
      plot_error_lst3r <- plot_error_lst3r1
    }

    if (ntraits > 1 & all(nrows > ntraits) & ext.ord == "random") {
      plot_error_lst3r <- mapply(function(x) scale(x %*% solve(chol(stats::var(x))) %*% chol(EcorR)),
        x = plot_error_lst3r, SIMPLIFY = FALSE
      )
    }
    if (ntraits == 1 | any(nrows <= ntraits) | ext.ord == "zig-zag") {
      plot_error_lst3r <- mapply(function(x) scale(x),
        x = plot_error_lst3r, SIMPLIFY = FALSE
      )
    }
    Zr <- lapply(seq_len(nenvs), function(i) stats::model.matrix(~ row - 1, droplevels(error_df[error_df$env == i, ])))
    plot_error_lst3r <- mapply(function(w, x) scale(w %*% x), w = Zr, x = plot_error_lst3r, SIMPLIFY = FALSE)
  }

  varR <- as.data.frame(t(matrix(varR, ncol = ntraits)))
  varR <- lapply(X = varR, FUN = c)
  prop.spatial <- as.data.frame(t(matrix(prop.spatial, ncol = ntraits)))
  prop.ext_col <- as.data.frame(t(matrix(prop.ext_col, ncol = ntraits)))
  prop.ext_row <- as.data.frame(t(matrix(prop.ext_row, ncol = ntraits)))
  if (ntraits > 1) {
    prop.spatial <- lapply(X = prop.spatial, FUN = diag)
    prop.ext_col <- lapply(X = prop.ext_col, FUN = diag)
    prop.ext_row <- lapply(X = prop.ext_row, FUN = diag)
  }
  if (ntraits == 1) {
    prop.spatial <- lapply(X = prop.spatial, FUN = diag, nrow = 1)
    prop.ext_col <- lapply(X = prop.ext_col, FUN = diag, nrow = 1)
    prop.ext_row <- lapply(X = prop.ext_row, FUN = diag, nrow = 1)
  }
  e_spat <- mapply(function(x, y) (scale(x) %*% sqrt(y)), x = plot_error_lst1, y = prop.spatial, SIMPLIFY = FALSE)
  e_rand <- mapply(function(w, x, y, z) (w %*% sqrt(round(diag(1, nrow = ntraits) - x - y - z, 8))), w = plot_error_lst2, x = prop.spatial, y = prop.ext_col, z = prop.ext_row, SIMPLIFY = FALSE)
  e_ext_c <- mapply(function(x, y) (x %*% sqrt(y)), x = plot_error_lst3c, y = prop.ext_col, SIMPLIFY = FALSE)
  e_ext_r <- mapply(function(x, y) (x %*% sqrt(y)), x = plot_error_lst3r, y = prop.ext_row, SIMPLIFY = FALSE)

  if (ntraits > 1) {
    e_scale <- mapply(function(w, x, y, z) sqrt(diag(1 / diag(as.matrix(stats::var(w + x + y + z))))), w = e_spat, x = e_rand, y = e_ext_c, z = e_ext_r, SIMPLIFY = FALSE)
    e_spat <- mapply(function(x, y, z) x %*% y %*% diag(sqrt(z)), x = e_spat, y = e_scale, z = varR, SIMPLIFY = FALSE)
    e_rand <- mapply(function(x, y, z) x %*% y %*% diag(sqrt(z)), x = e_rand, y = e_scale, z = varR, SIMPLIFY = FALSE)
    e_ext_c <- mapply(function(x, y, z) x %*% y %*% diag(sqrt(z)), x = e_ext_c, y = e_scale, z = varR, SIMPLIFY = FALSE)
    e_ext_r <- mapply(function(x, y, z) x %*% y %*% diag(sqrt(z)), x = e_ext_r, y = e_scale, z = varR, SIMPLIFY = FALSE)
  }

  if (ntraits == 1) {
    e_scale <- mapply(function(w, x, y, z) sqrt(1 / stats::var(w + x + y + z)), w = e_spat, x = e_rand, y = e_ext_c, z = e_ext_r, SIMPLIFY = FALSE)
    e_spat <- mapply(function(x, y, z) x %*% y %*% sqrt(z), x = e_spat, y = e_scale, z = varR, SIMPLIFY = FALSE)
    e_rand <- mapply(function(x, y, z) x %*% y %*% sqrt(z), x = e_rand, y = e_scale, z = varR, SIMPLIFY = FALSE)
    e_ext_c <- mapply(function(x, y, z) x %*% y %*% sqrt(z), x = e_ext_c, y = e_scale, z = varR, SIMPLIFY = FALSE)
    e_ext_r <- mapply(function(x, y, z) x %*% y %*% sqrt(z), x = e_ext_r, y = e_scale, z = varR, SIMPLIFY = FALSE)
  }

  plot_error_lst <- mapply(function(w, x, y, z) w + x + y + z, w = e_spat, x = e_rand, y = e_ext_c, z = e_ext_r, SIMPLIFY = FALSE)
  plot_error <- do.call(what = "rbind", plot_error_lst)
  colnames(plot_error) <- paste0("e.Trait", 1:ntraits)
  error_df <- cbind(error_df, plot_error)

  if (return.effects) {
    e_spat <- do.call("rbind", e_spat)
    e_rand <- do.call("rbind", e_rand)
    e_ext_c <- do.call("rbind", e_ext_c)
    e_ext_r <- do.call("rbind", e_ext_r)
    e_all <- lapply(seq_len(ncol(e_spat)), function(i) cbind(e_spat[, i], e_rand[, i], e_ext_c[, i], e_ext_r[, i]))
    effects_df <- lapply(e_all, function(x) {
      data.frame(error_df[, 1:4],
        e.spat = x[, 1],
        e.rand = x[, 2],
        e.ext.col = x[, 3],
        e.ext.row = x[, 4]
      )
    })

    list_names <- c("error.df", paste0("Trait", 1:ntraits))
    error_df <- list(error_df)
    error_df <- c(error_df, effects_df)
    names(error_df) <- list_names
  }
  error_df
  return(error_df)
}
