#' Simulate plot errors in plant breeding trials
#'
#' Creates a data frame with simulated plot errors for one or more traits in plant breeding
#' trials across one or more environments. The simulated error consists of a spatial error term,
#' a random error term and an extraneous error term. The spatial error term is constructed
#' according to either 1) bivariate interpolation using the \link[interp]{interp} function of
#' the package 'interp', or 2) a separable first-order autoregressive process (AR1:AR1). The
#' random error term is constructed using an independent process. The extraneous error term is
#' constructed as the sum of column and/or row terms.
#' The spatial, random and extraneous error terms are combined according to a user-defined ratio.
#' \cr
#' For multiple traits, correlated error terms can be generated assuming 1) correlated spatial
#' error between traits, 2) correlated random error between traits, 3) correlated extraneous
#' error between traits, or 4) some combination of 1-3. \cr
#' A separable covariance structure is assumed between traits and environments.
#'
#' @param n_envs Number of environments to be simulated (same as for \code{compsym_asr_input}
#'   or \code{unstr_asr_output}, where applicable).
#' @param n_traits Number of traits to be simulated.
#' @param n_reps A vector specifying the number of complete replicates in each environment. If only
#'   one value is provided and \code{n_envs > 1}, all environments will be assigned the same
#'   number of replicates.
#' @param n_cols A vector specifying the total number of columns in each environment. If only one
#'   value is provided and \code{n_envs > 1}, all environments will be assigned the same number
#'   of columns.
#' @param n_rows A vector specifying the total number of rows in each environment. If only one
#'   value is provided and \code{n_envs > 1}, all environments will be assigned the same number
#'   of rows.
#' @param rep_dir A character string specifying the direction of replicate blocks. One of either
#'   "column" (side-by-side, the default) or "row" (above-and-below). \code{rep_dir} is ignored
#'   when \code{n_reps = 1}.
#' @param var_R A vector of error variances for each trait by environment combination (ordered
#'   as environments within traits). If the length of \code{var_R} is equal to \code{n_traits},
#'   all environments will be assigned the same error variance for each trait.
#' @param S_cor_R A matrix of spatial error correlations between more than one trait. If not
#'   defined and \code{n_traits > 1}, a diagonal matrix is constructed.
#' @param R_cor_R A matrix of random error correlations between more than one trait. If not
#'   defined and \code{n_traits > 1}, a diagonal matrix is constructed.
#' @param E_cor_R A matrix of extraneous error correlations between more than one trait. If not
#'   defined and \code{n_traits > 1}, a diagonal matrix is constructed. The same correlation between 
#'   traits is used for the column and row errors (where applicable).
#' @param spatial_model A character string specifying the model used to simulate the two-dimensional
#'   spatial error term. One of either "Bivariate" (bivariate interpolation, the default) or "AR1:AR1"
#'   (separable first-order autoregressive process).
#' @param complexity A vector specifying the complexity of the bivariate interpolation in each
#'   environment. If only one value is provided and \code{n_envs > 1}, all environments will be
#'   assigned the same complexity. If not specified and \code{spatial_model = "Bivariate"}, the 
#'   complexity is set to the maximum number of columns and rows for each environment.
#'   This generally provides good results. See \link[interp]{interp} for further details.
#' @param plot_length A vector of plot lengths (column direction, usually longer side) for
#'   each environment. If only one value is provided and \code{n_envs > 1}, all environments will be 
#'   assigned the same plot length. Only required when \code{spatial_model = "Bivariate"}.
#' @param plot_width A vector of plot widths (row direction, usually shorter side) for
#'   each environment. If only one value is provided and \code{n_envs > 1}, all environments will be 
#'   assigned the same plot width. Only required when \code{spatial_model = "Bivariate"}.
#' @param col_cor A vector of column autocorrelations for each environment used in the AR1:AR1
#'   spatial error model. If only one value is provided, all environments will be assigned the
#'   same column autocorrelation. Only required when \code{spatial_model = "AR1:AR1"}.
#' @param row_cor A vector of row autocorrelations for each environment used in the AR1:AR1
#'   spatial error model. If only one value is provided, all environments will be assigned the
#'   same row autocorrelation. Only required when \code{spatial_model = "AR1:AR1"}.
#' @param prop_spatial A vector specifying the proportion of spatial error variance to total error
#'   variance (spatial + random + extraneous) for each trait by environment combination. If the
#'   length of \code{prop_spatial} is equal to \code{n_traits}, all environments will be assigned the same
#'   proportion for each trait By default, \code{prop_spatial = 0.5}.
#' @param ext_col_cor A vector of column autocorrelations for each environment used in the extraneous
#'   error model. If only one value is provided, all environments will be assigned the same column
#'   autocorrelation. By default, \code{ext_col_cor = 0}.
#' @param ext_row_cor A vector of row autocorrelations for each environment used in the extraneous
#'   error model. If only one value is provided, all environments will be assigned the same row
#'   autocorrelation. By default, \code{ext_row_cor = 0}.
#' @param prop_ext A vector specifying the proportion of extraneous error variance to total error
#'   variance (spatial + random + extraneous) for each trait by environment combination. If the
#'   length of \code{prop_ext} is equal to \code{n_traits}, all environments will be assigned the same
#'   same proportion for each trait. By default, \code{prop_ext = 0}.
#' @param ext_dir A character string specifying the direction of extraneous variation. One of either
#'   "column", "row" or "both". When "both", half the variance is assigned to the columns and half
#'   is assigned to the rows.
#' @param return_effects When TRUE, a list is returned with additional entries for each trait
#'   containing the spatial, random and extraneous errors. By default, return_effects = FALSE.
#'
#' @return A data frame containing the environment, block, column and row identifiers, as well as the
#'   simulated error for each trait. When \code{return_effects = TRUE}, a list is returned with
#'   additional entries for each trait containing the spatial and random error values.
#'
#' @examples
#' # Simulation of plot errors for two traits in three environments using a bivariate
#' # interpolation model for spatial variation.
#'
#' n_envs <- 3 # Number of simulated environments.
#' n_traits <- 2 # Number of simulated traits.
#'
#' # Field layout
#' n_cols <- 10 # Total number of columns in each environment.
#' n_rows <- c(30, 30, 20) # Total number of rows in each environment.
#' plot_length <- 5 # Plot length set to 5 meters in each environment.
#' plot_width <- 2 # Plot width set to 2 meters in each environment.
#' n_reps <- c(3, 3, 2) # Number of complete replicates (blocks) per environment.
#'
#' # Error variances for traits 1 and 2.
#' var_R <- c(0.4, 15)
#'
#' # Spatial error correlations between traits 1 and 2.
#' S_cor_R <- matrix(
#'   c(
#'     1.0, 0.2,
#'     0.2, 1.0
#'   ),
#'   ncol = 2
#' )
#'
#' error_df <- field_trial_error(
#'   n_envs = n_envs,
#'   n_traits = n_traits,
#'   n_reps = n_reps,
#'   n_cols = n_cols,
#'   n_rows = n_rows,
#'   rep_dir = "row",
#'   var_R = var_R,
#'   S_cor_R = S_cor_R,
#'   spatial_model = "Bivariate",
#'   plot_length = plot_length,
#'   plot_width = plot_width,
#'   prop_spatial = 0.4,
#'   ext_row_cor = -0.6,
#'   prop_ext = 0.2,
#'   ext_dir = "row",
#'   return_effects = TRUE
#' )
#' @export
field_trial_error <- function(n_envs,
                              n_traits,
                              n_reps,
                              n_cols,
                              n_rows,
                              rep_dir = "column",
                              var_R,
                              S_cor_R = NULL,
                              R_cor_R = NULL,
                              E_cor_R = NULL,
                              spatial_model = "bivariate",
                              complexity = NULL,
                              plot_length,
                              plot_width,
                              col_cor = NULL,
                              row_cor = NULL,
                              prop_spatial = 0.5,
                              ext_col_cor = 0,
                              ext_row_cor = 0,
                              prop_ext = 0,
                              ext_dir = NULL,
                              return_effects = FALSE) {
  if (n_envs < 1 | n_envs %% 1 != 0) stop("'n_envs' must be an integer > 0")
  if (n_traits < 1 | n_traits %% 1 != 0) stop("'n_traits' must be an integer > 0")

  if (min(n_cols) < 1 | any(n_cols %% 1 != 0)) {
    stop("'n_cols' must contain integers > 0")
  }
  if (length(n_cols) == 1) n_cols <- rep(n_cols, n_envs)
  if (length(n_cols) != n_envs) {
    stop("Length of 'n_cols' does not match the number of environments")
  }
  if (min(n_rows) < 1 | any(n_rows %% 1 != 0)) {
    stop("'n_rows' must contain integers > 0")
  }
  if (length(n_rows) == 1) n_rows <- rep(n_rows, n_envs)
  if (length(n_rows) != n_envs) {
    stop("Length of 'n_rows' does not match the number of environments")
  }

  if (min(plot_length) <= 0) stop("'plot_length' must contain values > 0")
  if (length(plot_length) == 1) plot_length <- rep(plot_length, n_envs)
  if (length(plot_length) != n_envs) {
    stop("Length of 'plot_length' does not match the number of environments")
  }
  if (min(plot_width) <= 0) stop("'plot_width' must be > 0")
  if (length(plot_width) == 1) plot_width <- rep(plot_width, n_envs)
  if (length(plot_width) != n_envs) {
    stop("Length of 'plot_width' does not match the number of environments")
  }

  if (min(n_reps) < 1 | any(n_reps %% 1 != 0)) {
    stop("'n_reps' must contain integers > 0")
  }

  rep_dir <- tolower(rep_dir)
  if (rep_dir == "column" & any((n_cols / n_reps) %% 1 != 0)) {
    stop("Number of columns not divisible by number of reps in at least one
         environment. Review your trial design!")
  }
  if (rep_dir == "row" & any((n_rows / n_reps) %% 1 != 0)) {
    stop("Number of rows not divisible by number of reps in at least one
           environment. Review your trial design!")
  }

  if (!rep_dir %in% c("column", "row")) {
    stop("'rep_dir' must be either 'column' or 'row'")
  }

  if (length(n_reps) == 1) n_reps <- rep(n_reps, n_envs)
  if (length(n_reps) != n_envs) {
    stop("Length of 'n_reps' does not match the number of environments")
  }

  if (length(var_R) == 1) var_R <- rep(var_R, n_traits)
  if (length(var_R) == n_traits) var_R <- rep(var_R, each = n_envs)
  if (any(var_R <= 0)) {
    stop("'var_R' must contain values greater than 0")
  }
  if (length(var_R) != (n_envs * n_traits)) {
    stop("Length of 'var_R' does not match the number of traits
         or the number of trait by environment combinations")
  }

  if (is.null(S_cor_R)) S_cor_R <- diag(n_traits)
  if (is.null(R_cor_R)) R_cor_R <- diag(n_traits)
  if (is.null(E_cor_R)) E_cor_R <- diag(n_traits)

  if (any(eigen(S_cor_R)$values < 1e-7)) {
    stop("'S_cor_R' is not postiive definite")
  }
  if (any(eigen(R_cor_R)$values < 1e-7)) {
    stop("'R_cor_R' is not postiive definite")
  }
  if (any(eigen(E_cor_R)$values < 1e-7)) {
    stop("'E_cor_R' is not postiive definite")
  }

  spatial_model <- tolower(spatial_model)
  if (spatial_model != "bivariate" & spatial_model != "ar1:ar1") {
    stop("'spatial_model' must be 'Bivariate' or 'AR1:AR1'")
  }

  if (any(prop_spatial < 0) | any(prop_spatial > 1)) {
    stop("'prop_spatial' must contain values between 0 and 1")
  }
  if (length(prop_spatial) == 1) prop_spatial <- rep(prop_spatial, n_traits)
  if (length(prop_spatial) == n_traits) prop_spatial <- rep(prop_spatial, each = n_envs)
  if (length(prop_spatial) != n_traits * n_envs) {
    stop("Length of 'prop_spatial' does not match the number of environments or the number of trait by environment combinations")
  }
  
  if(is.null(ext_dir) | prop_ext == 0 ){ext_dir <- "both"}
  if (any(prop_ext < 0) | any(prop_ext > 1)) {
    stop("'prop_ext' must contain values between 0 and 1")
  }
  if (length(prop_ext) == 1) prop_ext <- rep(prop_ext, n_traits)
  if (length(prop_ext) == n_traits) prop_ext <- rep(prop_ext, each = n_envs)
  if (length(prop_ext) != n_traits * n_envs) {
    stop("Length of 'prop_ext' does not match the number of environments or the number of trait by
         environment combinations")
  }
  if (any(prop_spatial + prop_ext > 1)) {
    stop("The sum of 'prop_spatial' and 'prop_ext' must be between 0 and 1")
  }
  if (any(prop_spatial + prop_ext == 1)) {
    warning("The random error is zero for some environments")
  }
  ext_dir <- tolower(ext_dir)
  if (any(prop_ext > 0) & !any(ext_dir %in% c("column", "row", "both"))) {
    stop("'ext_dir' must be either 'column', 'row' or 'both'")
  }

  if (ext_dir %in% c("column", "both") & any(n_cols <= 3 & any(prop_ext > 0))) {
    stop("'n_cols' must be greater than 3 when simulating column extraneous variation.")
  }

  if (ext_dir %in% c("row", "both") & any(n_rows <= 3 & prop_ext > 0)) {
    stop("'n_rows' must be greater than 3 when simulating row extraneous variation.")
  }

  envs <- rep(1:n_envs, times = n_cols * n_rows)
  reps <- c(unlist(mapply(function(x, y, z) rep(1:z, each = c(x * y / z)), x = n_cols, y = n_rows, z = n_reps)))
  cols <- c(unlist(mapply(function(x, y) rep(1:x, each = y), x = n_cols, y = n_rows)))
  rows <- c(unlist(mapply(function(x, y) rep(1:x, times = y), y = n_cols, x = n_rows)))

  if (rep_dir == "row") {
    cols <- c(unlist(mapply(function(x, y) rep(1:x, times = y), x = n_cols, y = n_rows)))
    rows <- c(unlist(mapply(function(x, y) rep(1:x, each = y), y = n_cols, x = n_rows)))
  }

  plot_df <- data.frame(
    env = factor(envs),
    block = factor(reps),
    col = factor(cols),
    row = factor(rows)
  )
  plot_df <- plot_df[order(plot_df$env, plot_df$col, plot_df$row), ]
  rownames(plot_df) <- NULL

  if (spatial_model == "ar1:ar1") {
        if (length(col_cor) == 1) col_cor <- rep(col_cor, n_envs)
    if (length(col_cor) != n_envs) {
      stop("Length of vector 'col_cor' does not match total number of environments")
    }
    if (any(col_cor < -1) | any(col_cor > 1)) {
      stop("'col_cor' must contain values between -1 and 1'")
    }
    if(any(abs(col_cor) == 1)) {col_cor[abs(col_cor) == 1] <- sign(col_cor[abs(col_cor) == 1])*(1-1e-7)}

    if (length(row_cor) == 1) row_cor <- rep(row_cor, n_envs)
    if (length(row_cor) != n_envs) {
      stop("Length of vector 'row_cor' does not match total number of environments")
    }
    if (any(row_cor < -1) | any(row_cor > 1)) {
      stop("row_cor must be between -1 and 1")
    }
    if(any(abs(row_cor) == 1)) {row_cor[abs(row_cor) == 1] <- sign(row_cor[abs(row_cor) == 1])*(1-1e-7)}

    power_lst <- lapply(n_cols, function(x) abs(outer(1:x, 1:x, "-")))
    col_ar1 <- mapply(function(x, y) x^y, x = col_cor, y = power_lst, SIMPLIFY = FALSE)

    power_lst <- lapply(n_rows, function(x) abs(outer(1:x, 1:x, "-")))
    row_ar1 <- mapply(function(x, y) x^y, x = row_cor, y = power_lst, SIMPLIFY = FALSE)

    cor_mat_lst <- mapply(function(x, y) kronecker(t(chol(x)), t(chol(y))), x = col_ar1, y = row_ar1, SIMPLIFY = FALSE)
    x <- T
    y <- 0
    while (x | y == 100) {
      y <- y + 1
      plot_error_lst1 <- mapply(function(x, y) y %*% scale(matrix(c(stats::rnorm(x)), ncol = n_traits)),
        x = n_cols * n_rows * n_traits, y = cor_mat_lst, SIMPLIFY = FALSE
      )

      x <- any(unlist(lapply(plot_error_lst1, function(x) eigen(stats::var(x))$values < 1e-7)))
      if (y == 100) {
        stop("Appropriate AR1:AR1 spatial error not obtained in 100 iterations.")
      }
    }
    if (n_traits > 1) {
      plot_error_lst1 <- mapply(function(x) scale(x %*% solve(chol(stats::var(x))) %*% chol(S_cor_R)),
        x = plot_error_lst1, SIMPLIFY = FALSE
      )
    }
    if (n_traits == 1) {
      plot_error_lst1 <- mapply(function(x) scale(x),
        x = plot_error_lst1, SIMPLIFY = FALSE
      )
    }
  }

  if (spatial_model == "bivariate") {
    if (is.null(complexity)) complexity <- apply(cbind(n_cols, n_rows), 1, max)
    if (length(complexity) == 1) complexity <- rep(complexity, n_envs)
    if (length(complexity) != n_envs) {
      stop("Length of 'complexity' does not match the number of environments")
    }
    if (any(complexity <= 0)) stop("'complexity' values must be integers > 0")

    cols_lst <- with(plot_df, tapply(col, env, function(x) c(1:max(as.numeric(trimws(x))), max(as.numeric(trimws(x))) + 1)))
    rows_lst <- with(plot_df, tapply(row, env, function(x) c(1:max(as.numeric(trimws(x))), max(as.numeric(trimws(x))) + 1)))
    max_lst <- mapply(function(x, y) 1:max(x, y), x = cols_lst, y = rows_lst, SIMPLIFY = FALSE)

    col_centres_lst <- mapply(function(x, y) round(y * (x - 0.5), 8),
      x = max_lst, y = plot_length, SIMPLIFY = FALSE
    )

    row_centres_lst <- mapply(function(x, y) round(y * (x - 0.5), 8),
      x = max_lst, y = plot_width, SIMPLIFY = FALSE
    )

    col_gap <- plot_length / 4
    row_gap <- plot_width / 4

    plot_error_lst1 <- NA
    y <- 0
    while (sum(is.na(unlist(plot_error_lst1))) > 0 | y == 100) {
      y <- y + 1
      xInterp_list <- mapply(function(w, x, y, z) c(0 - z, x * y + z, 0 - z, x * y + z, sample(stats::runif(n = w, min = 0, max = (x * y)))),
        w = complexity, x = n_cols, y = plot_length, z = col_gap, SIMPLIFY = FALSE
      )
      yInterp_list <- mapply(function(w, x, y, z) c(0 - z, 0 - z, x * y + z, x * y + z, sample(stats::runif(n = w, min = 0, max = (x * y)))),
        w = complexity, x = n_rows, y = plot_width, z = row_gap, SIMPLIFY = FALSE
      )
      x <- T
      w <- 0
      while (x | w == 100) {
        w <- w + 1
        zInterp_list <- mapply(function(w) scale(matrix(stats::rnorm((4 + w) * n_traits), ncol = n_traits)),
          w = complexity, SIMPLIFY = FALSE
        )

        x <- any(unlist(lapply(zInterp_list, function(x) eigen(stats::var(x))$values < 1e-7)))
        if (y == 100) {
          stop("Appropriate bivariate spatial error not obtained in 100 iterations.")
        }
      }
      if (n_traits > 1) {
        zInterp_list <- mapply(function(x) scale(x %*% solve(chol(stats::var(x))) %*% chol(S_cor_R)),
          x = zInterp_list, SIMPLIFY = FALSE
        )
      }

      if (n_traits == 1) {
        zInterp_list <- mapply(function(x) scale(x),
          x = zInterp_list, SIMPLIFY = FALSE
        )
      }

      for (i in 1:n_traits) {
        tmp <- mapply(
          function(v, w, x, y, z, xo, yo) {
            c(fill_matrix(t(interp::interp(x = x, y = y, z = z[, i], xo = xo, yo = yo, linear = F, extrap = T, duplicate = "mean")$z)[1:v, 1:w]))
          },
          v = n_rows, w = n_cols, x = xInterp_list, y = yInterp_list, z = zInterp_list, xo = col_centres_lst, yo = row_centres_lst, SIMPLIFY = FALSE
        )
        if (i == 1) {
          plot_error_lst1 <- tmp
        }
        if (i > 1) {
          plot_error_lst1 <- Map("cbind", plot_error_lst1, tmp)
        }
      }
      if (y == 100) {
        stop("Appropriate bivariate spatial error not obtained in 100 iterations. Consider reseting 'complexity' argument.")
      }
    }
  }

  x <- T
  y <- 0
  while (x | y == 100) {
    y <- y + 1
    plot_error_lst2 <- mapply(function(x) scale(matrix(c(stats::rnorm(x)), ncol = n_traits)),
      x = n_cols * n_rows * n_traits, SIMPLIFY = FALSE
    )
    x <- any(unlist(lapply(plot_error_lst2, function(x) eigen(stats::var(x))$values < 1e-7)))
    if (y == 100) {
      stop("Appropriate random error not obtained in 100 iterations.")
    }
  }

  if (n_traits > 1) {
    plot_error_lst2 <- mapply(function(x) scale(x %*% solve(chol(stats::var(x))) %*% chol(R_cor_R)),
      x = plot_error_lst2, SIMPLIFY = FALSE
    )
  }
  if (n_traits == 1) {
    plot_error_lst2 <- mapply(function(x) scale(x),
      x = plot_error_lst2, SIMPLIFY = FALSE
    )
  }

  n_plots <- mapply(function(x, y) x * y, x = n_cols, y = n_rows)
  plot_error_lst3c <- lapply(n_plots, function(x) matrix(0, nrow = x, ncol = n_traits))
  if (any(ext_dir %in% c("column", "both"))) {

    if (is.null(ext_col_cor)) ext_col_cor <- 0
    if (length(ext_col_cor) == 1) ext_col_cor <- rep(ext_col_cor, n_envs)
    if (length(ext_col_cor) != n_envs) {
      stop("Length of vector 'ext_col_cor' does not match total number of environments")
    }
    if (any(ext_col_cor < -1) | any(ext_col_cor > 1)) {
      stop("'ext_col_cor' must contain values between -1 and 1'")
    }
    if(any(abs(ext_col_cor) == 1)) {ext_col_cor[abs(ext_col_cor) == 1] <- sign(ext_col_cor[abs(ext_col_cor) == 1])*(1-1e-7)}

    power_lst <- lapply(n_cols, function(x) abs(outer(1:x, 1:x, "-")))
    ext_col_ar1 <- mapply(function(x, y) x^y, x = ext_col_cor, y = power_lst, SIMPLIFY = FALSE)
    ext_col_cor_mat_lst <- mapply(function(x) t(chol(x)), x = ext_col_ar1, SIMPLIFY = FALSE)

    plot_error_lst3c <- mapply(function(x,y) y %*% scale(matrix(c(stats::rnorm(x)), ncol = n_traits)),
                               x = n_cols * n_traits, y = ext_col_cor_mat_lst, SIMPLIFY = FALSE
    )
    x <- any(unlist(lapply(plot_error_lst3c, function(x) eigen(stats::var(x))$values < 1e-7)))
    y <- 0
    while (x & all(n_cols > n_traits) | y == 100 & all(n_cols > n_traits)) {
      plot_error_lst3c <- mapply(function(x,y) y %*% scale(matrix(c(stats::rnorm(x)), ncol = n_traits)),
                                 x = n_cols * n_traits, y = ext_col_cor_mat_lst, SIMPLIFY = FALSE
      )
      x <- any(unlist(lapply(plot_error_lst3c, function(x) eigen(stats::var(x))$values < 1e-7)))
      y <- y + 1
      if (y == 100) {
        stop("Appropriate column extraneous error not obtained in 100 iterations.")
      }
    }
    if (n_traits > 1 & all(n_cols > n_traits)) {
      plot_error_lst3c <- mapply(function(x) scale(x %*% solve(chol(stats::var(x))) %*% chol(E_cor_R)),
                                 x = plot_error_lst3c, SIMPLIFY = FALSE
      )
    }
    if (n_traits == 1 | any(n_cols <= n_traits)) {
      plot_error_lst3c <- mapply(function(x) scale(x),
                                 x = plot_error_lst3c, SIMPLIFY = FALSE
      )
    }
    Zc <- lapply(seq_len(n_envs), function(i) stats::model.matrix(~ col - 1, droplevels(plot_df[plot_df$env == i, ])))
    plot_error_lst3c <- mapply(function(w, x) scale(w %*% x), w = Zc, x = plot_error_lst3c, SIMPLIFY = F)
  }

  plot_error_lst3r <- lapply(n_plots, function(x) matrix(0, nrow = x, ncol = n_traits))
  if (any(ext_dir %in% c("row", "both"))) {

    if (is.null(ext_row_cor)) ext_row_cor <- 0
    if (length(ext_row_cor) == 1) ext_row_cor <- rep(ext_row_cor, n_envs)
    if (length(ext_row_cor) != n_envs) {
      stop("Length of vector 'ext_row_cor' does not match total number of environments")
    }
    if (any(ext_row_cor < -1) | any(ext_row_cor > 1)) {
      stop("'ext_row_cor' must contain values between -1 and 1'")
    }
    if(any(abs(ext_row_cor) == 1)) {ext_row_cor[abs(ext_row_cor) == 1] <- sign(ext_row_cor[abs(ext_row_cor) == 1])*(1-1e-7)}

    power_lst <- lapply(n_rows, function(x) abs(outer(1:x, 1:x, "-")))
    ext_row_ar1 <- mapply(function(x, y) x^y, x = ext_row_cor, y = power_lst, SIMPLIFY = FALSE)
    ext_row_cor_mat_lst <- mapply(function(x) t(chol(x)), x = ext_row_ar1, SIMPLIFY = FALSE)

    plot_error_lst3r <- mapply(function(x,y) y %*% scale(matrix(c(stats::rnorm(x)), ncol = n_traits)),
                               x = n_rows * n_traits, y = ext_row_cor_mat_lst, SIMPLIFY = FALSE
    )
    x <- any(unlist(lapply(plot_error_lst3r, function(x) eigen(stats::var(x))$values < 1e-7)))
    y <- 0
    while (x & all(n_rows > n_traits) | y == 100 & all(n_rows > n_traits)) {
      plot_error_lst3r <- mapply(function(x,y) y %*% scale(matrix(c(stats::rnorm(x)), ncol = n_traits)),
                                 x = n_rows * n_traits, y = ext_row_cor_mat_lst, SIMPLIFY = FALSE
      )
      x <- any(unlist(lapply(plot_error_lst3r, function(x) eigen(stats::var(x))$values < 1e-7)))
      y <- y + 1
      if (y == 100) {
        stop("Appropriate row extraneous error not obtained in 100 iterations.")
      }
    }
    if (n_traits > 1 & all(n_rows > n_traits)) {
      plot_error_lst3r <- mapply(function(x) scale(x %*% solve(chol(stats::var(x))) %*% chol(E_cor_R)),
                                 x = plot_error_lst3r, SIMPLIFY = FALSE
      )
    }
    if (n_traits == 1 | any(n_rows <= n_traits)) {
      plot_error_lst3r <- mapply(function(x) scale(x),
                                 x = plot_error_lst3r, SIMPLIFY = FALSE
      )
    }
    Zr <- lapply(seq_len(n_envs), function(i) stats::model.matrix(~ row - 1, droplevels(plot_df[plot_df$env == i, ])))
    plot_error_lst3r <- mapply(function(w, x) scale(w %*% x), w = Zr, x = plot_error_lst3r, SIMPLIFY = F)
  }

  var_R <- as.data.frame(t(matrix(var_R, ncol = n_traits)))
  var_R <- lapply(X = var_R, FUN = c)
  prop_spatial <- as.data.frame(t(matrix(prop_spatial, ncol = n_traits)))
  if (n_traits > 1) prop_spatial <- lapply(X = prop_spatial, FUN = diag)
  if (n_traits == 1) prop_spatial <- lapply(X = prop_spatial, FUN = diag, nrow = 1)
  prop_ext <- as.data.frame(t(matrix(prop_ext, ncol = n_traits)))
  if (n_traits > 1) prop_ext <- lapply(X = prop_ext, FUN = diag)
  if (n_traits == 1) prop_ext <- lapply(X = prop_ext, FUN = diag, nrow = 1)
  e_spat <- mapply(function(x, y) (scale(x) %*% sqrt(y)), x = plot_error_lst1, y = prop_spatial, SIMPLIFY = F)
  e_rand <- mapply(function(x, y, z) (x %*% sqrt(diag(1, nrow = n_traits) - y - z)), x = plot_error_lst2, y = prop_spatial, z = prop_ext, SIMPLIFY = F)
  e_ext_c <- mapply(function(x, y) (x %*% sqrt(y)), x = plot_error_lst3c, y = prop_ext, SIMPLIFY = F)
  e_ext_r <- mapply(function(x, y) (x %*% sqrt(y)), x = plot_error_lst3r, y = prop_ext, SIMPLIFY = F)

  if (any(ext_dir == "both")) {
    e_ext_c <- lapply(e_ext_c, function(x) sqrt(0.5) * x)
    e_ext_r <- lapply(e_ext_r, function(x) sqrt(0.5) * x)
  }

  if (n_traits > 1) {
    e_scale <- mapply(function(w, x, y, z) sqrt(diag(1 / diag(as.matrix(stats::var(w + x + y + z))))), w = e_spat, x = e_rand, y = e_ext_c, z = e_ext_r, SIMPLIFY = F)
    e_spat <- mapply(function(x, y, z) x %*% y %*% diag(sqrt(z)), x = e_spat, y = e_scale, z = var_R, SIMPLIFY = F)
    e_rand <- mapply(function(x, y, z) x %*% y %*% diag(sqrt(z)), x = e_rand, y = e_scale, z = var_R, SIMPLIFY = F)
    e_ext_c <- mapply(function(x, y, z) x %*% y %*% diag(sqrt(z)), x = e_ext_c, y = e_scale, z = var_R, SIMPLIFY = F)
    e_ext_r <- mapply(function(x, y, z) x %*% y %*% diag(sqrt(z)), x = e_ext_r, y = e_scale, z = var_R, SIMPLIFY = F)
  }

  if (n_traits == 1) {
    e_scale <- mapply(function(x, y) sqrt(1 / stats::var(x + y)), x = e_spat, y = e_rand, SIMPLIFY = F)
    e_spat <- mapply(function(x, y, z) x %*% y %*% sqrt(z), x = e_spat, y = e_scale, z = var_R, SIMPLIFY = F)
    e_rand <- mapply(function(x, y, z) x %*% y %*% sqrt(z), x = e_rand, y = e_scale, z = var_R, SIMPLIFY = F)
    e_ext_c <- mapply(function(x, y, z) x %*% y %*% sqrt(z), x = e_ext_c, y = e_scale, z = var_R, SIMPLIFY = F)
    e_ext_r <- mapply(function(x, y, z) x %*% y %*% sqrt(z), x = e_ext_r, y = e_scale, z = var_R, SIMPLIFY = F)
  }

  plot_error_lst <- mapply(function(w, x, y, z) w + x + y + z, w = e_spat, x = e_rand, y = e_ext_c, z = e_ext_r, SIMPLIFY = F)
  plot_error <- do.call(what = "rbind", plot_error_lst)
  colnames(plot_error) <- paste0("e.Trait.", 1:n_traits)
  plot_df <- cbind(plot_df, plot_error)

  if (return_effects) {
    e_spat <- do.call("rbind", e_spat)
    e_rand <- do.call("rbind", e_rand)
    e_ext_c <- do.call("rbind", e_ext_c)
    e_ext_r <- do.call("rbind", e_ext_r)
    e_all <- lapply(seq_len(ncol(e_spat)), function(i) cbind(e_spat[, i], e_rand[, i], e_ext_c[, i], e_ext_r[, i]))
    resids <- lapply(e_all, function(x) {
      data.frame(plot_df[, 1:4],
        e_spat = x[, 1],
        e_rand = x[, 2],
        e_ext_col = x[, 3],
        e_ext_row = x[, 4]
      )
    })

    list_names <- c("plot_df", paste0("Trait.", 1:n_traits))
    plot_df <- list(plot_df)
    plot_df <- c(plot_df, resids)
    names(plot_df) <- list_names
  }
  plot_df
  return(plot_df)
}
